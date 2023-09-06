module impute_mod

   use geostat
   use kdtree2_module
   use mtmod, only: grnd, gaussrnd
   use network_mod, only: network_forward
   use types_mod, only: variogram, network, kdtrees
   use vario_mod, only: set_sill, set_rotmatrix
   use covasubs, only: get_cov
   use constants, only: MINCOV, IMPEPS, EPSLON, SMALLDBLE
   use subs, only: shuffle, gauinv, sortem, nscore, &
                   locate, powint, dblecumsum

   implicit none

contains

   subroutine impute()

      real(8) :: start, finish

      call cpu_time(start)

      write (*, *) " "
      write (*, *) "Starting imputation..."
      write (*, *) " "

      call nmr_imputer(pool, nreals, nsearch, imputed, iter1, iter2, tol1, tol2)

      call cpu_time(finish)

      write (*, *) " "
      print '("Imputation took ", f5.2, " minutes")', (finish - start)/60

   end subroutine impute

   subroutine nmr_imputer(pool, nreals, nsearch, imputed, iter1, iter2, tol1, tol2)

      ! sequential Gaussian rejection imputation

      type(variogram), intent(inout) :: pool(:)
      integer, intent(in) :: nreals, nsearch
      integer, intent(in) :: iter1, iter2
      real(8), intent(in) :: tol1, tol2
      real(8), allocatable, intent(out) :: imputed(:, :, :) ! (ndata, nfact, nreals)

      ! kdtree required variables
      !   type(kdtree2), pointer :: tree
      type(kdtrees) :: trees(ngvarg) ! array of pointers
      type(kdtree2_result), allocatable :: results(:)
      integer :: nfound

      ! normal equations
      real(8), allocatable :: rhs(:), lhs(:, :), kwts(:)
      real(8), allocatable :: cmeans(:), cstdevs(:)
      integer, allocatable :: nuse(:, :), useidx(:, :)

      ! imputation
      real(8), allocatable :: sim(:, :)
      real(8), allocatable :: anisxyz(:, :)
      integer, allocatable :: isim(:, :), randpath(:)
      integer, allocatable :: factpath(:)
      real(8) :: zimp1(1), yimp1(1, ngvarg + 1)
      real(8) :: zimp2(1), yimp2(1, ngvarg + 1)
      real(8) :: ztry(1), ytry(1, ngvarg + 1)
      real(8) :: p, xp, gmin, gmax
      real(8) :: axyz(3)
      integer :: simidx, ierr
      real(8) :: diff1, diff2, pert
      integer :: nfact

      ! reference distribution
      integer, parameter :: nsamp = 1000000
      integer :: nlook
      real(8) :: zref(nsamp)
      real(8) :: yref(nsamp, ngvarg + 1)
      real(8), allocatable :: nsref(:)
      real(8) :: usnsref(nsamp)
      real(8) :: minref, maxref

      ! dynamic tolerances
      real(8), allocatable :: vsort(:), wsort(:), vord(:)
      real(8), allocatable :: vcdf(:)
      real(8) :: qv, qdiff
      integer :: qidx, qj

      ! indexes
      integer :: i, j, k1, k2, igv, ix, iy, fi, nst, ireal, luidx

      ! allocate arrays based on number of data and max search
      allocate (rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch))
      allocate (nuse(ndata, ngvarg), useidx(ndata, ndata))
      allocate (sim(ndata, ngvarg + 1), isim(ndata, ngvarg), randpath(ndata))
      allocate (results(ndata), anisxyz(3, ndata))
      allocate (cmeans(ngvarg + 1), cstdevs(ngvarg + 1)) ! +1 for nugget
      allocate (imputed(ndata, ngvarg + 2, nreals)) ! +2 for nugget and zval
      allocate (zinit(ndata, nreals)) ! +2 for nugget and zval

      ! build reference distribution for transformations
      call build_refcdf(nsamp, zref, nsref, yref, usnsref) ! this also calc norm moments
      minref = minval(zref)
      maxref = maxval(zref)
      write (*, *) "Reference min and max:", minref, maxref
      write (*, *) "NScore min and max:", minval(nsref), maxval(nsref)
      write (*, *) " "

      ! write out the reference dist and cumulative probs.
      do i = 1, nsamp
         write (ldbg, "(*(g14.8,1x))") dble(i)/(dble(nsamp) + 1.d0), zref(i), nsref(i)
      end do

      ! build declustered data cdf for quantile lookup
      vsort = var
      wsort = wts/sum(wts)
      vord = [(i, i=1, ndata)]
      call sortem(1, ndata, vsort, 2, wsort, vord, vsort, vsort, vsort, vsort, &
                  vsort, vsort)
      vcdf = dblecumsum(wsort)
      vcdf = vcdf(2:)
      vcdf = vcdf - vcdf(1)/2.0

      ! min and max Gaussian values ~ [-4, 4]
      call gauinv(0.0001d0, gmin, ierr)
      call gauinv(0.9999d0, gmax, ierr)

      ! total number of factors, plus one for nugget
      nfact = ngvarg + 1

      ! hard code number of structs to be 1 for each Gaussian
      nst = 1

      ! set rotmats if not already done
      do igv = 1, ngvarg
         if (.not. allocated(pool(igv)%rm)) then
            call set_rotmatrix(pool)
         end if
      end do

      ! setup anisotropic kdtree for each factor
      do igv = 1, ngvarg
         do i = 1, ndata
            anisxyz(:, i) = matmul(pool(igv)%rm(:, :, nst), xyz(:, i))
         end do
         trees(igv)%tree => kdtree2_create(input_data=anisxyz, dim=3, &
                                           sort=.true., rearrange=.true.)
      end do

      !
      ! main loop over realizations
      !
      REALLOOP: do ireal = 1, nreals

         write (*, *) "imputing realization", ireal

         ! define random path through nodes
         randpath = [(i, i=1, ndata)]
         call shuffle(randpath)

         ! all location are initially unsimulated
         isim = 0
         useidx = 0
         nuse = 0
         nlook = 0

         !
         ! loop over data locations
         !
         DATALOOP: do i = 1, ndata

            ! get random simulation index
            simidx = randpath(i)

            !
            ! loop over factors at this simidx
            !
            FACTLOOP: do igv = 1, ngvarg

               ! anisotropic coord for this simidx
               axyz = matmul(pool(igv)%rm(:, :, nst), xyz(:, simidx))

               ! query the tree for this location
               call kdtree2_r_nearest(tp=trees(igv)%tree, qv=axyz, &
                                      r2=pool(igv)%aa(nst)**2, nfound=nfound, &
                                      nalloc=ndata, results=results)
               !
               ! loop over samples found in search
               !
               nuse(simidx, igv) = 0
               do j = 1, nfound

                  ! check if this data index is the simulation index
                  if (results(j)%idx .eq. simidx) cycle

                  ! check if this data index is simulated or not
                  if (isim(results(j)%idx, igv) .eq. 0) cycle ! no conditioning value here

                  ! meet minimum covariance? (not collocated)
                  rhs(1) = get_cov(pool(igv), xyz(:, simidx), xyz(:, results(j)%idx))
                  if (rhs(1) .lt. MINCOV) cycle

                  ! if we got this far increment number found at ith data location
                  nuse(simidx, igv) = nuse(simidx, igv) + 1

                  ! track conditioning indices found at ith data location
                  useidx(nuse(simidx, igv), simidx) = results(j)%idx

                  ! have we met the max search?
                  if (nuse(simidx, igv) .ge. nsearch) exit

               end do

               ! build and solve normal equations
               if (nuse(simidx, igv) .gt. 0) then
                  call krige(pool(igv), xyz, rhs, lhs, kwts, nuse(:, igv), useidx, &
                             sim(:, igv), simidx, cmeans(igv), cstdevs(igv))
               else
                  ! if no data the distribution is N(0,1)
                  cmeans(igv) = 0.d0
                  cstdevs(igv) = 1.d0
               end if

               ! min stdev to prevent getting stuck?
               if (cstdevs(igv) .lt. 0.1) cstdevs(igv) = 0.1

               ! update this location and factor as simulated
               isim(simidx, igv) = 1

            end do FACTLOOP

            ! cmean and cstdev for nugget
            cmeans(ngvarg + 1) = 0.d0
            cstdevs(ngvarg + 1) = 1.d0

            !
            ! dyanmic tolerances
            !
            ! get the quantile of the variable we are trying to match
            ! so we can find the corresponding activation in the refdist
            !
            call locate(vsort, ndata, 1, ndata, var(simidx), qj)
            qj = min(max(1, qj), (ndata - 1))
            qv = powint(vsort(qj), vsort(qj + 1), vcdf(qj), vcdf(qj + 1), &
                        var(simidx), 1.d0)

            ! find the index in nsamp that corresponds to qv
            qidx = int(nsamp*qv)

            ! get the difference between adjacent samples at this
            ! position in the cdf
            qdiff = abs(zref(qidx + 1) - zref(qidx))
            qdiff = (qdiff + abs(zref(qidx - 1) - zref(qidx)))/2.d0

            !
            ! now that we have the conditional moments for each factor,
            ! repeatedly simulate values and check against the first tolerance
            !
            ! coarse search
            !
            diff1 = 999.0
            k1 = 0

            COARSE: do while (diff1 .gt. tol1)

               k1 = k1 + 1

               ! simulate each factor at this simidx
               do igv = 1, ngvarg + 1

                  p = grnd()
                  call gauinv(p, xp, ierr)
                  sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)

                  ! enforce reasonable min/max Gaussian values
                  do while ((sim(simidx, igv) .lt. gmin) .or. &
                            (sim(simidx, igv) .gt. gmax))
                     ! resimulate if necessary
                     p = grnd()
                     call gauinv(p, xp, ierr)
                     sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)
                  end do

               end do

               ! calculate imputed value
               yimp1 = reshape(sim(simidx, :), shape=[1, nfact]) ! reshape to preserve first dim
               call network_forward(nnet, yimp1, zimp1, nstrans=.false., &
                                    norm=nnet%norm, calc_mom=.false.)

               zinit(simidx, ireal) = zimp1(1)

               call transform_to_refcdf(zimp1(1), zref, nsref, zimp1(1))

               ! get difference with true data value
               diff1 = abs(zimp1(1) - var(simidx))

               ! break if we need to
               if (k1 .ge. iter1) then
                  write (*, *) "coarse search did not converge after", iter1, &
                     "iterations at data index", simidx, "diff=", diff1

                  do j = 1, nfact
                     write (*, *) cmeans(j), cstdevs(j)
                  end do

                  write (*, *) " "
                  write (*, *) yimp1
                  write (*, *) " "

                  exit
               end if

            end do COARSE

            !
            ! yimp1 contains imputed values which meet the first tolerance; now
            ! polish the previous solution and check against the second tolerance
            !
            ! current diff b/w imputed and true values
            diff2 = diff1

            ! current imputation values
            zimp2 = zimp1
            yimp2 = yimp1
            ytry = yimp1

            !
            ! solution polishing
            !
            k2 = 0
            POLISH: do while (diff2 .gt. tol2)

               k2 = k2 + 1

               ! random path though the factors
               factpath = [(fi, fi=1, nfact - 1)]
               call shuffle(factpath)

               do j = 1, nfact - 1

                  iy = factpath(j)

                  ! perturb a factor
                  pert = -qdiff + grnd()*(qdiff - (-qdiff))
                  ytry(1, iy) = yimp2(1, iy) + pert

                  ! revert to previous value if the new one is extreme
                  if (ytry(1, iy) .lt. gmin) then
                     ! write (*, *) "extreme value in polish", ytry(1, iy)
                     ytry(1, iy) = yimp2(1, iy)
                  end if

                  if (ytry(1, iy) .gt. gmax) then
                     ! write (*, *) "extreme value in polish", ytry(1, iy)
                     ytry(1, iy) = yimp2(1, iy)
                  end if

                  ! calculate the new imputed value
                  call network_forward(nnet, ytry, ztry, nstrans=.false., &
                                       norm=nnet%norm, calc_mom=.false.)
                  zinit(simidx, ireal) = ztry(1)
                  call transform_to_refcdf(ztry(1), zref, nsref, ztry(1))

                  ! did this pertubation improve the solution?
                  if (abs(ztry(1) - var(simidx)) .lt. diff2) then
                     ! update imputed values
                     yimp2(1, iy) = ytry(1, iy)
                     zimp2 = ztry
                     ! update the new difference
                     diff2 = abs(ztry(1) - var(simidx))
                  else
                     ! revert back to the previous value
                     ytry(1, iy) = yimp2(1, iy)
                  end if

               end do

               ! break if we need to
               if (k2 .ge. iter2) then
                  write (*, *) " "
                  write (*, *) "solution polishing did not converge after", iter2, &
                     "iterations at data index", simidx, "diff=", diff2

                  write (*, *) yimp2
                  write (*, *) " "

                  exit
               end if

            end do POLISH

            ! !
            ! ! If the polish doesn't converge, grab the correct zval and corresponding
            ! ! factors from the reference look up table. If there are only a few
            ! ! samples that dont converge, variogram reproduction shound't be affected
            ! ! too much.
            ! !
            ! ! this tolerance is arbitrary (could be tol2?)
            ! if (diff2 .gt. tol2) then

            !    ! grab NS z value and corresponding factors from table
            !    call lookup_table(var(simidx), usnsref, yref, zimp2, yimp2, luidx)
            !    zinit(simidx, ireal) = zref(luidx)

            !    ! write (*, *) "conditional means:", cmeans
            !    ! write (*, *) "conditional stdevs:", cstdevs

            !    ! track lookups
            !    nlook = nlook + 1

            ! end if

            !
            ! update the conditioning values for this realization
            !
            sim(simidx, :) = yimp2(1, :)

            !
            ! store the final values
            !
            imputed(simidx, 1:nfact, ireal) = yimp2(1, :)
            imputed(simidx, nfact + 1, ireal) = zimp2(1)

         end do DATALOOP

         if (nlook .gt. 0) write (*, *) nlook, "values taken from lookup table"

      end do REALLOOP

   end subroutine nmr_imputer

   subroutine krige(vm, coords, rhs, lhs, kwts, nuse, useidx, &
                    sim, simidx, cmean, cstdev)

      ! simple kriging conditional mean and variance

      integer, intent(in) :: nuse(:), useidx(:, :), simidx
      real(8), intent(in) :: coords(:, :), sim(:)
      type(variogram) :: vm
      real(8), intent(inout) :: rhs(:), lhs(:, :), kwts(:)
      real(8), intent(inout) :: cmean, cstdev
      integer :: j, k, test

      ! calculate matrices for normal equations
      do j = 1, nuse(simidx)
         ! build rhs vector
         rhs(j) = get_cov(vm, coords(:, simidx), coords(:, useidx(j, simidx)))
         do k = j, nuse(simidx)
            ! diagonal
            if (j .eq. k) then
               lhs(j, j) = 1.d0
            else
               ! build lhs matrix
               lhs(j, k) = get_cov(vm, coords(:, useidx(k, simidx)), &
                                   coords(:, useidx(j, simidx)))
               if (lhs(j, k) .lt. 0) stop 'ERROR: Negative covariance.'
               lhs(k, j) = lhs(j, k)
            end if
         end do
      end do

      ! solve the kriging system - external call to LAPACK
      call solve(lhs(1:nuse(simidx), 1:nuse(simidx)), kwts(1:nuse(simidx)), &
                 rhs(1:nuse(simidx)), nuse(simidx), 1, test)

      ! calcualte conditional mean and standard devaition
      cmean = 0.d0
      cstdev = get_cov(vm, coords(:, simidx), coords(:, simidx)) ! variance
      do j = 1, nuse(simidx)
         ! cmean considers previously simulated values
         cmean = cmean + kwts(j)*sim(useidx(j, simidx))
         cstdev = cstdev - kwts(j)*rhs(j)
      end do
      cstdev = sqrt(max(cstdev, 0.0))

   end subroutine krige

   subroutine build_refcdf(nsamp, zref, nsref, yref, usnsref)

      ! build reference CDF for the normal score transform
      ! of imputed values

      integer, intent(in) :: nsamp ! num samples in ref CDF
      real(8), intent(inout) :: zref(nsamp) ! sorted zref
      real(8), allocatable, intent(inout) :: nsref(:) ! sorted NS values
      real(8), intent(inout) :: usnsref(nsamp) ! unsorted NS values
      real(8), intent(inout) :: yref(nsamp, ngvarg + 1) ! Gauss lookup table
      real(8), allocatable :: tmp(:)
      real(8) :: p, xp, wt(nsamp)
      integer :: i, j, ierr, nfact

      nfact = ngvarg + 1
      wt = 1.d0

      ! initialize independent N(0,1) y realizations
      do j = 1, nfact
         do i = 1, nsamp
            p = grnd()
            call gauinv(p, xp, ierr)
            yref(i, j) = xp
         end do
         ! write (*, *) minval(yref(:, j)), maxval(yref(:, j))
      end do

      ! calculate the corresponding z values
      call network_forward(nnet, yref, zref, nstrans=.false., &
                           norm=nnet%norm, calc_mom=.true.)

      ! normal score to build transform table
      do i = 1, nsamp
         zref(i) = zref(i) + grnd()*SMALLDBLE ! random despike
      end do
      call nscore(nsamp, zref, dble(-1e21), dble(1e21), 1, wt, &
                  tmp, nsref, ierr)

      ! keep an unsorted NS copy for lookup table
      usnsref = nsref

      ! sort the reference distribution (for interp.) and return
      call sortem(1, nsamp, zref, 1, nsref, zref, zref, zref, &
                  zref, zref, zref, zref)

   end subroutine build_refcdf

   subroutine transform_to_refcdf(z, zref, nsref, zimp)

      ! normal score transform network output according to
      ! specified reference distribution

      real(8), intent(in) :: z ! network output
      real(8), intent(in) :: zref(:) ! sorted zvalues
      real(8), intent(in) :: nsref(:) ! corresponding sorted NS values
      real(8), intent(out) :: zimp ! N(0,1) imputed value
      integer :: j, nsamp

      nsamp = size(zref, dim=1)
      call locate(zref, nsamp, 1, nsamp, z, j)
      j = min(max(1, j), (nsamp - 1))
      zimp = powint(zref(j), zref(j + 1), nsref(j), nsref(j + 1), z, 1.d0)

   end subroutine transform_to_refcdf

   subroutine lookup_table(z, zref, yref, zimp, yimp, idx)

      real(8), intent(in) :: z ! target data value
      real(8), intent(in) :: zref(:) ! unsorted reference z values
      real(8), intent(in) :: yref(:, :) ! unsorted Gauss factors corresponding to zref
      real(8), intent(out) :: zimp(:) ! new imputed value from lookup table
      real(8), intent(out) :: yimp(:, :) ! Gauss factors corresponding to zimp
      real(8) :: diff(size(zref))
      integer, intent(out) :: idx

      ! find the minimum diff between target and table
      diff = abs(zref - z)
      idx = minloc(diff, dim=1)

      ! grab the values from that index
      zimp = zref(idx)
      yimp(1, :) = yref(idx, :)

   end subroutine lookup_table

end module impute_mod
