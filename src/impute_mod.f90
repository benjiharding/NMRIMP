module impute_mod

   use geostat
   use kdtree2_module
   use mtmod, only: grnd, gaussrnd
   use network_mod, only: network_forward
   use types_mod, only: variogram, network, kdtrees
   use vario_mod, only: set_sill, set_rotmatrix
   use covasubs, only: get_cov
   use constants, only: MINCOV, IMPEPS, EPSLON
   use subs, only: shuffle, gauinv, sortem, nscore, &
                   locate, powint

   implicit none

contains

   subroutine impute()

      real(8) :: start, finish

      call cpu_time(start)

      call nmr_imputer(pool, nreals, nsearch, imputed, iter1, iter2, tol1, tol2)

      call cpu_time(finish)

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
      real(8) :: zimp1(1), yimp1(1, ngvarg + 1)
      real(8) :: zimp2(1), yimp2(1, ngvarg + 1)
      real(8) :: ztry(1), ytry(1, ngvarg + 1)
      real(8) :: p, xp, axyz(3)
      integer :: simidx, ierr
      real(8) :: diff1, diff2, pert
      integer :: nfact

      ! reference distribution
      integer, parameter :: nsamp = 100000
      real(8) :: zref(nsamp)
      real(8), allocatable :: nsref(:)

      ! indexes
      integer :: i, j, k1, k2, igv, iy, nst, ireal

      ! allocate arrays based on number of data and max search
      allocate (rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch))
      allocate (nuse(ndata, ngvarg), useidx(ndata, ndata))
      allocate (sim(ndata, ngvarg + 1), isim(ndata, ngvarg), randpath(ndata))
      allocate (results(ndata), anisxyz(3, ndata))
      allocate (cmeans(ngvarg), cstdevs(ngvarg))
      allocate (imputed(ndata, ngvarg + 2, nreals)) ! +2 for nugget and zval

      ! build reference distribution for transformations
      call build_refcdf(nsamp, zref, nsref)
      write (*, *) minval(zref), maxval(zref)
      write (*, *) minval(nsref), maxval(nsref)

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

         ! define random path through nodes
         randpath = [(i, i=1, ndata)]
         call shuffle(randpath)

         ! all location are initially unsimulated
         isim = 0
         useidx = 0
         nuse = 0

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

               ! update this location and factor as simulated
               isim(simidx, igv) = 1

            end do FACTLOOP

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
               do igv = 1, ngvarg
                  p = grnd()
                  call gauinv(p, xp, ierr)
                  sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)
               end do

               ! simulate nugget
               p = grnd()
               call gauinv(p, xp, ierr)
               sim(simidx, ngvarg + 1) = xp

               ! calculate imputed value
               yimp1 = reshape(sim(simidx, :), shape=[1, nfact]) ! reshape to preserve first dim
               call network_forward(nnet, yimp1, zimp1, .false.)
               call transform_to_refcdf(zimp1(1), zref, nsref, zimp1(1))

               ! get difference with true data value
               diff1 = abs(zimp1(1) - var(simidx))

               ! break if we need to
               if (k1 .ge. iter1) then
                  write (*, *) "coarse search did not converge after", iter1, &
                     "iterations at data index", simidx
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

            k2 = 0
            POLISH: do while (diff2 .gt. tol2)

               k2 = k2 + 1

               do iy = 1, nfact

                  ! perturb a factor and calculate new imputed value
                  pert = -IMPEPS + grnd()*(2*IMPEPS)
                  ytry(1, iy) = yimp2(1, iy) + pert
                  call network_forward(nnet, ytry, ztry, .false.)
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
                  write (*, *) "solution polishing did not converge after", iter2, &
                     "iterations at data index", simidx
                  exit
               end if

            end do POLISH

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

   subroutine build_refcdf(nsamp, zref, nsref)

      ! build reference CDF for the normal score transform
      ! of imputed values

      integer, intent(in) :: nsamp
      real(8), intent(inout) :: zref(nsamp)
      real(8), allocatable, intent(inout) :: nsref(:)
      real(8), allocatable :: tmp(:)
      real(8) :: yrand(nsamp, ngvarg + 1), wt(nsamp)
      real(8) :: p, xp
      integer :: i, j, ierr, nfact

      nfact = ngvarg + 1
      wt = 1.d0

      ! initialize independent N(0,1) realizations
      do i = 1, nsamp
         do j = 1, nfact
            p = grnd()
            call gauinv(p, xp, ierr)
            yrand(i, j) = xp
         end do
      end do

      ! calculate the corresponding z values
      call network_forward(nnet, yrand, zref, nstrans=.false.)

      ! normal score to build transform table
      do i = 1, nsamp
         zref(i) = zref(i) + grnd()*EPSLON ! random despike
      end do
      call nscore(nsamp, zref, dble(-1e21), dble(1e21), 1, wt, &
                  tmp, nsref, ierr)

      ! sort the reference distribution and return
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

end module impute_mod
