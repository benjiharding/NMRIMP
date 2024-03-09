module impute_mod

   use geostat
   use kdtree2_module
   use mtmod, only: grnd
   use network_mod, only: network_forward, network_forward2
   use types_mod, only: variogram, network, kdtrees
   use vario_mod, only: set_sill, set_rotmatrix, gammabar
   use covasubs, only: get_cov
   use constants, only: MINCOV, EPSLON, SMALLDBLE, MAXRESIM
   use subs, only: shuffle, gauinv, sortem, powint, locate, nscore, reverse

   implicit none

contains

   subroutine impute()

      !
      ! convience subroutine to call imputer and handle any
      ! resimulation if required
      !

      real(8) :: start, finish
      integer :: ireal, nresim, i, j, luidx, nlook
      integer, allocatable :: rsidxs(:), idxs(:) ! resim indices

      call cpu_time(start)

      write (*, *) " "
      write (*, *) "Initializing imputer..."
      write (*, *) " "

      !
      ! intitalize the imputer
      !
      call init_imputer()

      write (*, *) " "
      write (*, *) "Starting imputation..."
      write (*, *) " "

      !
      ! main loop over realizations
      !
      REALLOOP: do ireal = 1, nreals

         write (*, *) "imputing realization", ireal

         j = 0
         nresim = 0
         call nmr_imputer(pool, ireal, nsearch, imputed, iter1, iter2, &
                          tol1, tol2, nresim, rsidxs)

         ! do we need to resim any locations?
         ! nresim is updated by call to imputer
         do while (nresim .gt. 0)

            j = j + 1
            if (j .gt. MAXRESIM) exit

            if (idbg .gt. 1) write (*, *) "resimulating at", nresim, "locations"

            ! reset the sim indicators for these locations
            idxs = rsidxs
            isim(rsidxs, :) = 0

            ! call the imputer again
            call nmr_imputer(pool, ireal, nsearch, imputed, iter1, iter2, &
                             tol1, tol2, nresim, rsidxs, idxs, j)

         end do

         ! use lookup table here?
         nlook = 0
         if (any(isim(:, 1) .lt. 1)) then
            do i = 1, ndata
               if (isim(i, 1) .lt. 1) then
                  call lookup_table(var(i), usnsref, yref, &
                                    imputed(i, nfact + 1, ireal), &
                                    imputed(i, 1:nfact, ireal), luidx)
                  zinit(i, ireal) = zref(luidx)
                  nlook = nlook + 1
                  if (idbg .ge. 1) then
                     write (*, *) "value taken from lookup table at data index ", i
                  end if
               end if
            end do
         end if

         if (nlook .gt. int(ndata*0.01)) then
            write (*, *) "Warning: more than 1% of values taken from lookup table"
         end if

      end do REALLOOP

      call cpu_time(finish)

      write (*, *) " "
      print '("Imputation took ", f5.2, " minutes")', (finish - start)/60

   end subroutine impute

   subroutine init_imputer()

      !
      ! initialize the NMR imputer by building reference distributions
      ! and allocating the required arrays
      !

      integer :: i, igv, n

      ! allocate arrays based on number of data and max search
      allocate (rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch))
      allocate (nuse(ndata, ngvarg), useidx(ndata, ndata))
      allocate (sim(ndata, ngvarg + 1), isim(ndata, ngvarg), randpath(ndata))
      allocate (results(ndata), anisxyz(3, ndata))
      allocate (cmeans(ngvarg + 1), cstdevs(ngvarg + 1)) ! +1 for nugget
      allocate (imputed(ndata, ngvarg + 2, nreals)) ! +2 for nugget and zval
      allocate (zinit(ndata, nreals))
      allocate (yref(nsamp, ngvarg + 1))
      allocate (trees(ngvarg))
      allocate (gbar(ngvarg), gbaridx(ngvarg))
      allocate (abratio(ndata, nreals), rsc(ndata, nreals))
      allocate (seeded(ndata, nreals))
      allocate (lb(size(var, dim=1)), ub(size(var, dim=1))) ! constraints

      ! intization simulation indicator
      isim = 0
      rsc = 0
      seeded = 0

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

      ! min and max Gaussian values
      gmin = minval(nsref)
      gmax = maxval(nsref)

      ! imputation bounds
      lb = gmin
      ub = gmax

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

      ! sort precedence array for indexing
      if (ifp) then
         fidx = [(i, i=1, ngvarg + 1)]
         fprec_r = fprec
         call sortem(1, ngvarg + 1, fprec_r, 1, fidx, fidx, fidx, fidx, &
                     fidx, fidx, fidx, fidx)
      end if

   end subroutine init_imputer

   subroutine nmr_imputer(pool, ireal, nsearch, imputed, iter1, iter2, &
                          tol1, tol2, nresim, rsidxs, idxs, irs)
      !
      ! sequential Gaussian rejection imputation
      !

      ! inputs
      type(variogram), intent(inout) :: pool(:)
      integer, intent(in) :: ireal, nsearch
      integer, intent(in) :: iter1, iter2
      real(8), intent(in) :: tol1, tol2
      real(8), intent(out) :: imputed(:, :, :) ! (ndata, nfact, nreals)
      integer, intent(out) :: nresim
      integer, allocatable, intent(out) :: rsidxs(:)
      integer, optional, intent(in) :: idxs(:) ! resim data indices
      integer, optional, intent(in) :: irs ! resim realization index

      ! imputation
      real(8) :: zimp1, yimp1(1, ngvarg + 1)
      real(8) :: zimp2, yimp2(1, ngvarg + 1)
      real(8) :: ytry(1, ngvarg + 1)
      real(8) :: p, xp, simt
      real(8) :: axyz(3)
      integer :: simidx, ierr
      real(8) :: diff1, diff2
      integer :: nfound

      ! factor sensitivity
      real(8), allocatable :: dp(:), dn(:)
      integer, allocatable :: didx(:)
      real(8) :: dy, dz, a, ab
      logical :: inv(ngvarg + 1)

      ! resimulation
      ! integer, parameter :: nreset = 6 ! simidx + 5 neighbours
      integer, parameter :: nreset = 11 ! simidx + 10 neighbours
      integer :: nr ! number of resims
      integer :: reset_idx(nreset)
      integer :: que(ndata*nreset) ! max possible size
      logical :: iresim

      ! indexes
      integer :: i, j, jj, k1, k2, igv, idx

      ! latent sensitivity arrays
      allocate (dp(ngvarg + 1), dn(ngvarg + 1), didx(ngvarg + 1))

      ! resimulaiton flag
      iresim = .false.
      if (nresim .gt. 0) iresim = .true.

      ! all location are initially unsimulated
      if (.not. iresim) isim = 0 ! dont reset if resimulating
      useidx = 0
      nuse = 0

      ! reset resim counter every realization
      nr = 0

      ! get the semi-random or random path for this iteration
      if (isd) then
         ! deallocate the previous path as size may change
         if (allocated(rp1s)) deallocate (rp1, rp1s, rp2)
         call semirandom_path(ireal)
      else
         ! true random path
         randpath = [(i, i=1, ndata)]
         call shuffle(randpath)
      end if

      !
      ! loop over data locations
      !
      DATALOOP: do i = 1, ndata

         ! get random simulation index
         simidx = randpath(i)

         ! consider a subset of indices?
         if (present(idxs)) then
            if (.not. any(idxs .eq. simidx)) cycle
         end if

         !
         ! loop over factors at this simidx
         !
         cmeans = 0.d0
         cstdevs = 1.d0
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
               ! if no conditioning the distribution is N(0,1)
               cmeans(igv) = 0.d0
               cstdevs(igv) = 1.d0
            end if

            ! min stdev to prevent getting stuck?
            if (cstdevs(igv) .lt. 0.1) cstdevs(igv) = 0.1

            ! crazy means?
            if (cmeans(igv) .lt. gmin) cmeans(igv) = gmin
            if (cmeans(igv) .gt. gmax) cmeans(igv) = gmax

         end do FACTLOOP

         ! cmean and cstdev for nugget
         cmeans(ngvarg + 1) = 0.d0
         cstdevs(ngvarg + 1) = 1.d0

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
            do j = 1, ngvarg + 1

               if (ifp .and. isd) then
                  ! impute in order of precedence
                  igv = fidx(j)
                  ! accept/reject based on threshold
                  ! this only applies to factor with precedence
                  if (j .eq. 1) then
                     p = grnd()
                     call gauinv(p, xp, ierr)
                     simt = xp*cstdevs(igv) + cmeans(igv)
                     jj = 0
                     ! enforce the constraints
                     do while ((simt .lt. lb(simidx)) .or. &
                               (simt .gt. ub(simidx)))
                        jj = jj + 1
                        p = grnd()
                        call gauinv(p, xp, ierr)
                        simt = xp*cstdevs(igv) + cmeans(igv)
                        if (jj .gt. 1000) then
                           ! write (*, *) "Stuck in coarse resim!"
                           exit
                        end if
                     end do
                     sim(simidx, igv) = simt
                  else
                     p = grnd()
                     call gauinv(p, xp, ierr)
                     sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)
                  end if
               else
                  ! no precedence
                  igv = j
                  p = grnd()
                  call gauinv(p, xp, ierr)
                  sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)
               end if

               ! enforce reasonable min/max Gaussian values
               do while ((sim(simidx, igv) .lt. gmin) .or. &
                         (sim(simidx, igv) .gt. gmax))
                  ! resimulate if necessary
                  p = grnd()
                  call gauinv(p, xp, ierr)
                  sim(simidx, igv) = xp*cstdevs(igv) + cmeans(igv)
               end do

            end do

            ! calculate observed value for this y vector
            yimp1 = reshape(sim(simidx, :), shape=[1, nfact]) ! reshape to preserve first dim
            zimp1 = f(yimp1)
            diff1 = abs(zimp1 - var(simidx))

            ! check for funky behaviour that will cause stress
            ! in the polishing step

            ! break if we need to
            if (k1 .ge. iter1) then
               if (idbg .gt. 1) then
                  write (*, *) "coarse search did not converge after", iter1, &
                     "iterations at data index", simidx, "diff=", diff1
               end if
               exit
            end if

         end do COARSE

         ! update conditioning array with intial coarse values
         zinit(simidx, ireal) = nmr(yimp1)
         sim(simidx, :) = yimp1(1, :)

         !
         ! if the coarse search doesn't converge track the locations
         ! for resimulation
         !
         if (diff1 .ge. tol1) then

            ! dont reset the seeded locations
            if (isd) then
               if (any(rp2 .eq. simidx)) then
                  ! track number of reset locations
                  nr = nr + 1
                  rsc(simidx, ireal) = rsc(simidx, ireal) + 1
                  ! get neighbouring indices to reset
                  do j = 1, nreset
                     reset_idx(j) = results(j)%idx
                  end do
                  ! track data indices so we can revisit them
                  que((nr - 1)*nreset + 1:nr*nreset) = reset_idx
               end if

            else
               ! track number of reset locations
               nr = nr + 1
               rsc(simidx, ireal) = rsc(simidx, ireal) + 1
               ! get neighbouring indices to reset
               do j = 1, nreset
                  reset_idx(j) = results(j)%idx
               end do
               ! track data indices so we can revisit them
               que((nr - 1)*nreset + 1:nr*nreset) = reset_idx
            end if

            ! cycle to start a new iteration of DATALOOP
            if (.not. present(irs)) cycle

            ! are we resimulating?
            if (present(irs)) then
               ! cycle to start a new iteration of DATALOOP unless
               ! we are on the last resim iteration
               if (irs .lt. MAXRESIM) cycle
            end if
            rsc(simidx, ireal) = rsc(simidx, ireal) - 1

         end if

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
         ab = 0.d0
         POLISH: do while (diff2 .gt. tol2)

            k2 = k2 + 1

            ! get the direction we need to move and factor sensitivities
            dz = var(simidx) - zimp2
            call latent_sensitivity(zimp2, yimp2, dz, dp, dn, simidx, &
                                    didx, a, dy, inv)
            ! track the final ratio
            ab = a/abs(dz)

            ! zimp is too high, we need to reduce
            if (dz .lt. 0.d0) then

               ! are we inside the sensitivity tolerances for negative dz?
               ! if so we can simply find a solution
               if (a .lt. 0.d0) then

                  ! find the factor with the smallest delta that
                  ! is bigger than dz
                  call findidx(dp, dn, dz, 0, idx)

                  ! do some binary search based on dp(idx)
                  ! this updates diff2 which should break do while condition
                  call binary_search(var(simidx), yimp2, dy, inv(didx(idx)), &
                                     didx(idx), tol2, zimp2, diff2)

               else if (a .gt. 0.d0) then
                  ! if we are outside the senesitvity tolerances for negative
                  ! dz, we need to perturb and reasses the sensitivities

                  ! sort by gamma bar?
                  if (isort .gt. 0) then
                     do j = 1, nsort
                        call set_to_bound(dy, dp, dn, yimp2, int(gbaridx(j)), 0, zimp2)
                     end do
                  else ! set the most sensitive factor to bound and update z;
                     zimp2 = dynamic_set_to_bound(dy, dp, dn, yimp2, didx, 0)
                  end if

                  ! check the new difference
                  diff2 = abs(zimp2 - var(simidx))
               end if

            else if (dz .gt. 0.d0) then
               ! zimp is too low, we need to increase

               ! are we inside the sensitivity tolerances for postive dz?
               ! if so we can simply find a solution
               if (a .lt. 0.d0) then

                  ! find the factor with the smallest delta that
                  ! is bigger than dz
                  call findidx(dp, dn, dz, 1, idx)

                  ! do some binary search based on dp(idx)
                  ! this updates diff2 which should break do while condition
                  call binary_search(var(simidx), yimp2, dy, inv(didx(idx)), &
                                     didx(idx), tol2, zimp2, diff2)

               else if (a .gt. 0.d0) then
                  ! if we are outside the senesitvity tolerances for
                  ! positive dz, we need to perturb and reasses the sensitivities

                  ! sort by gamma bar?
                  if (isort .gt. 0) then
                     do j = 1, nsort
                        call set_to_bound(dy, dp, dn, yimp2, int(gbaridx(j)), 1, zimp2)
                     end do
                  else ! set the most sensitive factor to bound and update z;
                     zimp2 = dynamic_set_to_bound(dy, dp, dn, yimp2, didx, 1)
                  end if

                  ! check the new difference
                  diff2 = abs(zimp2 - var(simidx))

               end if

            end if

            ! break if we need to
            if (k2 .ge. iter2) then
               if (idbg .gt. 1) then
                  write (*, *) " "
                  write (*, *) "solution polishing did not converge after", iter2, &
                     "iterations at data index", simidx, "diff=", diff2
                  ! write (*, *) yimp2
                  ! write (*, *) " "
               end if
               exit
            end if

         end do POLISH

         ! polishing difficulty
         abratio(simidx, ireal) = ab

         !
         ! update the conditioning values for this realization
         !
         zinit(simidx, ireal) = nmr(yimp2)
         sim(simidx, :) = yimp2(1, :)

         !
         ! store the final values
         !
         imputed(simidx, 1:nfact, ireal) = yimp2(1, :)
         imputed(simidx, nfact + 1, ireal) = zimp2

         ! update this location and factor as simulated
         do igv = 1, ngvarg
            isim(simidx, igv) = 1
            ! unless we didnt converge
            if (diff2 .gt. tol2) isim(simidx, igv) = 0 ! flag for lookup
         end do

      end do DATALOOP

      ! get the resim indices from the que array if required
      nresim = nr
      if (nr .gt. 0) then
         allocate (rsidxs(nr*nreset))
         rsidxs = que(1:nr*nreset)
      else
         allocate (rsidxs(1))
         rsidxs = 0
      end if

   end subroutine nmr_imputer

   subroutine semirandom_path(ireal)

      integer, intent(in) :: ireal

      real(8), allocatable :: tvar(:)
      type(kdtree2), pointer :: nntree ! euclidean nearest neighbour tree
      type(kdtree2_result), allocatable :: nnresults(:)
      real(8), allocatable :: txyz(:, :)
      integer, allocatable :: tseed(:), exclude(:, :)
      integer :: nfound
      real(8) :: r2

      integer :: i, j, n, m

      ! define semi-random path through nodes to seed values
      n = 0
      do i = 1, ndata
         if (var(i) .gt. tfp) n = n + 1
      end do
      allocate (rp1(n), rp1_r(n), tvar(n))
      ! reset counters and track indices
      n = 0
      m = 0
      do i = 1, ndata
         if (var(i) .gt. tfp) then
            n = n + 1
            rp1(n) = i
            tvar(n) = var(i) + grnd()*0.025 ! small random component for sorting
         end if
      end do
      ! sort rp1 (locs > threshold) by tvar (grade > threshold)
      rp1_r = real(rp1)
      call sortem(1, n, tvar, 1, rp1_r, rp1_r, rp1_r, rp1_r, rp1_r, rp1_r, rp1_r, rp1_r)
      call reverse(tvar) ! descending
      call reverse(rp1_r) ! descending
      rp1 = int(rp1_r)

      ! build a euclidean kdtree for neighbour search
      allocate (txyz(3, n))
      do i = 1, n
         txyz(:, i) = xyz(:, rp1(i))
      end do
      nntree => kdtree2_create(input_data=txyz, dim=3, &
                               sort=.true., rearrange=.true.)

      ! iterate over rp1 and reject samples within a distance threshold
      allocate (tseed(n), exclude(n, n), nnresults(n))
      tseed = 0
      exclude = 0
      r2 = sr**2
      do i = 1, n
         if (any(exclude .eq. rp1(i))) cycle
         tseed(i) = rp1(i)
         call kdtree2_r_nearest(tp=nntree, qv=txyz(:, i), &
                                r2=r2, nfound=nfound, &
                                nalloc=n, results=nnresults)
         do j = 1, nfound
            if (rp1(nnresults(j)%idx) .eq. rp1(i)) cycle
            if (nnresults(j)%dis .lt. r2) then
               exclude(j - 1, i) = rp1(nnresults(j)%idx)
            else
               cycle
            end if
         end do
      end do

      ! get the final seed indices above the threshold
      allocate (rp1s(count(tseed > 0)))
      n = 0
      do i = 1, size(tseed)
         if (tseed(i) .gt. 0) then
            n = n + 1
            rp1s(n) = tseed(i)
         end if
      end do
      allocate (rp2(ndata - n))
      ! get all the others indices for the complete path
      n = 0
      m = 0
      do i = 1, ndata
         if (any(rp1s .eq. i)) then
            n = n + 1
            seeded(rp1s(n), ireal) = 1
            cycle
         else
            m = m + 1
            rp2(m) = i
         end if
      end do

      if (idbg .gt. 1) then
         write (*, *) "seeding", n, "locations above the", qfp, "quantile"
      end if

      ! data above threshold are simulated first
      call shuffle(rp1s) ! sparse seed locations > threhshold
      call shuffle(rp2) ! every other data location
      randpath(1:n) = rp1s
      randpath(n + 1:ndata) = rp2

      ! setup upper/lower thresholds for constrained imputation
      ! constraints only apply to the factor with precedence at seed locations
      do i = 1, n
         if (var(rp1s(i)) .gt. tfp) lb(rp1s(i)) = tfp
         if (var(i) .lt. tfp) ub(i) = tfp
      end do

      ! deallocate arrays as size may change across realizations
      deallocate (rp1_r, tvar, txyz)
      deallocate (tseed, exclude, nnresults)

   end subroutine semirandom_path

   subroutine krige(vm, coords, rhs, lhs, kwts, nuse, useidx, &
                    sim, simidx, cmean, cstdev)
      !
      ! simple kriging conditional mean and variance
      !

      type(variogram), intent(in) :: vm
      integer, intent(in) :: nuse(:), useidx(:, :), simidx
      real(8), intent(in) :: coords(:, :), sim(:)
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

   subroutine latent_sensitivity(z, y, dz, dp, dn, sidx, didx, a, dy, inv)

      !
      ! determine what latent factors have the greatest
      ! influence on zimp
      !

      ! paramters
      real(8), intent(in) :: z ! current imputed value
      real(8), intent(in) :: y(:, :) ! latent vector
      real(8), intent(in) :: dz ! current delta with observed
      real(8), intent(inout) :: dp(:), dn(:)
      integer, intent(in) :: sidx ! simidx
      integer, intent(out) :: didx(:) ! sorted delta z
      real(8), intent(out) :: a ! measure of in/out tols
      real(8), intent(out) :: dy ! delta y
      logical, intent(out) :: inv(:) ! perturbations inverse to sign?

      ! locals
      real(8) :: d(ngvarg + 1)
      real(8) :: ytmp(size(y, dim=1), size(y, dim=2))
      real(8) :: yy(size(y, dim=1), size(y, dim=2))
      real(8) :: y1, y2, zt
      real(8) :: b, c
      real(8) :: ffidx(ngvarg + 1)
      integer :: i, ierr

      ! factor indices
      ffidx = [(i, i=1, ngvarg + 1)]

      ! maximum step size in probability space
      call gauinv(0.5d000, y1, ierr)
      call gauinv(0.525d0, y2, ierr)
      dy = y2 - y1

      ! If inv is false postive perts increase z, negative
      ! decrease z. If true, sign of pert and the direction
      ! z moves are inverse
      inv = .false.

      ! calculate +/- sensitivity
      dn = 0.d0
      dp = 0.d0
      do i = 1, ngvarg + 1 ! exclude the nugget here? should always be least sensitive.
         ytmp = y
         ytmp(1, i) = ytmp(1, i) - dy
         zt = f(ytmp)
         dn(i) = z - zt
         ytmp(1, i) = ytmp(1, i) + dy*2.d0 ! revert back to zero plus dy
         zt = f(ytmp)
         dp(i) = z - zt
         if (dn(i) .lt. 0.d0) inv(i) = .true.
      end do

      ! ! set sensitivity of seeded values to zero to prevent polishing
      ! if (any(rp1 .eq. sidx)) then
      !    dn(int(fidx(1))) = 0.d0
      !    dp(int(fidx(1))) = 0.d0
      ! end if

      ! catch edge cases where dp/dn are all the same sign
      ! use a larger/smaller perturbation?
      if (dz .lt. 0) then
         if (all(dp < 0) .and. all(dn < 0)) then ! all neg
            dp = 0.d0
            dn = 0.d0
            do i = 1, ngvarg + 1
               ytmp = y
               ytmp(1, i) = ytmp(1, i) + dy
               dp(i) = z - f(ytmp)
               ytmp(1, i) = ytmp(1, i) - dy*2.d0 ! reset
               dn(i) = z - f(ytmp)
            end do
         end if
      else if (dz .gt. 0) then
         if (all(dp > 0) .and. all(dn > 0)) then ! all pos
            dp = 0.d0
            dn = 0.d0
            do i = 1, ngvarg + 1
               ytmp = y
               ytmp(1, i) = ytmp(1, i) - dy
               dn(i) = z - f(ytmp)
               if (dn(i) .lt. 0.d0) inv(i) = .true.
               ytmp(1, i) = ytmp(1, i) + dy*2.d0 ! reset
               dp(i) = z - f(ytmp)
            end do
         end if
      end if

      ! calc metrics of "difficulty"
      if (dz .lt. 0.d0) c = max(maxval(dp), maxval(dn))
      if (dz .gt. 0.d0) c = min(minval(dp), minval(dn))
      ! these absolute values make a < 0 inside regardless of the
      ! sign of dz
      b = abs(dz)
      a = b - abs(c)

      ! sort indices by delta - least to most sensitive
      d = abs(dp - dn)
      call sortem(1, ngvarg, d, 3, ffidx, dp, dn, d, d, d, d, d)
      didx = ffidx

   end subroutine latent_sensitivity

   subroutine binary_search(z, y, dy, inv, idx, tol, zhat, delta)

      !
      ! search for the value of y(idx) such that f(y) matches z;
      ! the correct value lies within y(idx) +/- dy
      ! https://codereview.stackexchange.com/questions/168730/binary-search-on-real-space
      !

      ! parameters
      real(8), intent(in) :: z, dy, tol
      integer, intent(in) :: idx
      logical, intent(in) :: inv
      real(8), intent(inout) :: y(:, :)
      real(8), intent(out) :: zhat, delta

      ! locals
      real(8) :: ytmp(size(y, dim=1), size(y, dim=2))
      real(8) :: k, m, n
      integer :: i

      ! set search boundaries
      m = y(1, idx) - dy
      n = y(1, idx) + dy

      if (m .lt. gmin) m = gmin
      if (n .gt. gmax) n = gmax

      ! check lower bound
      ytmp = y
      ytmp(1, idx) = m
      delta = z - f(ytmp)
      if (abs(delta) .le. tol) then
         y(1, idx) = m
         zhat = f(ytmp)
         delta = abs(delta)
         return
      end if

      ! check upper bound
      ytmp = y
      ytmp(1, idx) = n
      delta = z - f(ytmp)
      if (abs(delta) .le. tol) then
         y(1, idx) = n
         zhat = f(ytmp)
         delta = abs(delta)
         return
      end if

      ! binary search within [m, n]
      ytmp = y
      do while (m < n)
         i = i + 1
         k = (n + m)/2
         ytmp(1, idx) = k
         delta = z - f(ytmp)
         if (abs(delta) .le. tol) then
            y(1, idx) = k
            zhat = f(ytmp)
            delta = abs(delta)
            return
         end if
         if (delta .gt. 0) then
            if (inv) n = k ! decrease upper bound
            if (.not. inv) m = k ! increase lower bound
         else if (delta .lt. 0) then
            if (inv) m = k ! increase lower bound
            if (.not. inv) n = k ! decrease upper bound
         end if

         if (i .gt. 10000) then
            ! write (*, *) "infinte loop in binary search"
            return
         end if
      end do

   end subroutine binary_search

   subroutine build_refcdf(nsamp, zref, nsref, yref, usnsref)

      !
      ! build reference CDF for the normal score transform
      ! of imputed values
      !

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
      if (ifp) then
         call network_forward2(nnet, yref, zref, nstrans=.false., fprec=fprec, &
                               sigwt=sigwt)
      else
         call network_forward(nnet, yref, zref, nstrans=.false.)
      end if

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

      !
      ! normal score transform network output according to
      ! specified reference distribution
      !

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
      real(8), intent(out) :: zimp ! new imputed value from lookup table
      real(8), intent(out) :: yimp(:) ! Gauss factors corresponding to zimp
      real(8) :: diff(size(zref))
      integer, intent(out) :: idx

      ! find the minimum diff between target and table
      diff = abs(zref - z)
      idx = minloc(diff, dim=1)

      ! grab the values from that index
      zimp = zref(idx)
      yimp = yref(idx, :)

   end subroutine lookup_table

   subroutine findidx(dp, dn, dz, updn, idx)

      real(8), intent(in) :: dp(:), dn(:), dz
      integer, intent(in) :: updn
      integer, intent(out) :: idx

      ! locals
      integer :: i, j, k
      logical :: jj(size(dp)), kk(size(dp))

      jj = .false.
      kk = .false.

      if (updn .eq. 0) then
         ! if we are decreasing z the two scenarios are:
         ! 1. dp is pos and inv is true
         ! 2. dn is pos and inv is false

         do i = 1, size(dp)
            if (dp(i) .gt. 0.d0) then
               if (dp(i) .gt. abs(dz)) then
                  jj(i) = .true.
               end if
            end if
            if (dn(i) .gt. 0.d0) then
               if (dn(i) .gt. abs(dz)) then
                  kk(i) = .true.
               end if
            end if
         end do

         ! minimum index of values greater than dz
         j = minloc(dp, dim=1, mask=jj)
         k = minloc(dn, dim=1, mask=kk)

      else if (updn .eq. 1) then
         ! if we are increasing z the two scenarios are:
         ! 1. dp is neg and inv is false
         ! 2. dn is neg and inv is true

         do i = 1, size(dp)
            if (dp(i) .lt. 0.d0) then
               if (abs(dp(i)) .gt. dz) then
                  jj(i) = .true.
               end if
            end if
            if (dn(i) .lt. 0.d0) then
               if (abs(dn(i)) .gt. dz) then
                  kk(i) = .true.
               end if
            end if
         end do

         ! minimum index of values greater than dz
         j = maxloc(dp, dim=1, mask=jj)
         k = maxloc(dn, dim=1, mask=kk)
      end if

      ! these should never both be zero as a condition of
      ! entring this subroutine is at least one factor is
      ! sensitive enough
      if (j .eq. 0 .and. k .eq. 0) then
         stop "STOP! something is funky"
      end if

      ! kk is all false
      if (any(jj) .and. .not. any(kk)) then
         idx = j
         ! jj is all false
      else if (any(kk) .and. .not. any(jj)) then
         idx = k
      else
         idx = min(j, k)
      end if

   end subroutine findidx

   function findloc(arr, x, ltgt) result(idx)

      ! find the first occurence where the
      ! condition is true

      real(8) :: arr(:), x
      integer :: ltgt ! .lt. (0) or .gt. (1)
      integer :: j, idx

      select case (ltgt)
      case (0) ! find index where arr(j) < x
         do j = 1, size(arr)
            if (arr(j) .lt. x) then
               idx = j
               return
            end if
         end do
         idx = size(arr, dim=1) ! default is max
      case (1) ! find index where arr(j) > x
         do j = 1, size(arr)
            if (arr(j) .gt. x) then
               idx = j
               return
            end if
         end do
         idx = size(arr, dim=1) ! default is max
      end select

   end function findloc

   function dynamic_set_to_bound(dy, dp, dn, y, fidx, updn) result(z)

      ! determine the most sensitive factor to and set it
      ! to its sensitivity upper/lower bound

      real(8) :: dy, dp(:), dn(:), y(:, :)
      integer :: fidx(:)
      integer :: updn ! increase (1) or decrease (0) z
      real(8) :: z
      integer :: i, j
      logical :: ii(size(dp)), jj(size(dp))

      ii = .false.
      jj = .false.

      if (updn .eq. 0) then
         ! if we are decreasing z the two scenarios are:
         ! 1. dp is pos and inv is true
         ! 2. dn is pos and inv is false

         where (dp .gt. 0.d0) ii = .true.
         where (dn .gt. 0.d0) jj = .true.

         ! edge case where both masks are false?
         if (.not. any(ii) .and. .not. any(jj)) then
            i = 1
            j = 1
            ! write (*, *) "both masks are false"
         else
            i = maxloc(dp, dim=1, mask=ii)
            j = maxloc(dn, dim=1, mask=jj)
         end if

         if (dp(i) .gt. dn(j)) then
            y(1, fidx(i)) = y(1, fidx(i)) + dy
            if (y(1, fidx(i)) .gt. gmax) y(1, fidx(i)) = gmax
         else
            y(1, fidx(j)) = y(1, fidx(j)) - dy
            if (y(1, fidx(j)) .lt. gmin) y(1, fidx(j)) = gmin
         end if

      else if (updn .eq. 1) then
         ! if we are increasing z the two scenarios are:
         ! 1. dp is neg and inv is false
         ! 2. dn is neg and inv is true

         where (dp .lt. 0.d0) ii = .true.
         where (dn .lt. 0.d0) jj = .true.

         ! edge case where both masks are false?
         if (.not. any(ii) .and. .not. any(jj)) then
            i = 1
            j = 1
            ! write (*, *) "both masks are false"
         else
            i = minloc(dp, dim=1, mask=ii)
            j = minloc(dn, dim=1, mask=jj)
         end if

         if (dp(i) .lt. dn(j)) then
            y(1, fidx(i)) = y(1, fidx(i)) + dy
            if (y(1, fidx(i)) .gt. gmax) y(1, fidx(i)) = gmax
         else
            y(1, fidx(j)) = y(1, fidx(j)) - dy
            if (y(1, fidx(j)) .lt. gmin) y(1, fidx(j)) = gmin
         end if

      end if

      ! calculate the new imputed value
      z = f(y)

   end function dynamic_set_to_bound

   subroutine set_to_bound(dy, dp, dn, y, idx, updn, z)

      !
      ! Set the specified factor index to its upper/lower sensitivity bound.
      ! This subroutine modifies the y vector so multiple sucessive
      ! adjustments can be made.
      !

      real(8) :: dy, dp(:), dn(:), y(:, :)
      integer :: idx
      integer :: updn ! increase (1) or decrease (0) z
      real(8) :: z

      if (updn .eq. 0) then
         ! if we are decreasing z the two scenarios are:
         ! 1. dp is pos and inv is true
         ! 2. dn is pos and inv is false

         if (dp(idx) .gt. dn(idx)) then
            y(1, idx) = y(1, idx) + dy
         else
            y(1, idx) = y(1, idx) - dy
         end if

      else if (updn .eq. 1) then
         ! if we are increasing z the two scenarios are:
         ! 1. dp is neg and inv is false
         ! 2. dn is neg and inv is true

         if (dp(idx) .lt. dn(idx)) then
            y(1, idx) = y(1, idx) + dy
         else
            y(1, idx) = y(1, idx) - dy
         end if

      end if

      ! calculate the new imputed value
      z = f(y)

   end subroutine set_to_bound

   function f(y) result(z)

      ! convenience nmr function - latent to observed mapping

      real(8), intent(in) :: y(:, :)
      real(8) :: zz(1), z

      if (ifp) then
         call network_forward2(nnet, y, zz, nstrans=.false., fprec=fprec, &
                               sigwt=sigwt)
      else
         call network_forward(nnet, y, zz, nstrans=.false.)
      end if
      call transform_to_refcdf(zz(1), zref, nsref, zz(1))
      z = zz(1)

   end function f

   function nmr(y) result(a)

      ! conveience function for activation values

      real(8), intent(in) :: y(:, :)
      real(8) :: aa(1), a

      if (ifp) then
         call network_forward2(nnet, y, aa, nstrans=.false., fprec=fprec, &
                               sigwt=sigwt)
      else
         call network_forward(nnet, y, aa, nstrans=.false.)
      end if
      a = aa(1)

   end function nmr

end module impute_mod
