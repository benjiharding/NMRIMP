module impute_mod

   use geostat
   use mtmod
   use kdtree2_module
   use network_mod, only: network_forward
   use types_mod, only: variogram, network, kdtrees
   use vario_mod, only: set_sill, set_rotmatrix
   use covasubs, only: get_cov
   use constants, only: MINCOV
   use subs, only: shuffle

   implicit none

contains

   subroutine nmr_imputer(pool, nreals, nsearch, imputed, iter1, iter2, tol1, tol2)

      ! sequential Gaussian rejection imputation

      type(variogram), intent(inout) :: pool(:)
      integer, intent(in) :: nreals, nsearch
      integer, intent(in) :: iter1, iter2
      real(8), intent(in) :: tol1, tol2
      real(8), allocatable, intent(out) :: imputed(:, :)

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
      real(8) :: zimp(1)
      real(8) :: p, xp, axyz(3)
      integer :: simidx, ierr
      real(8) :: diff1, diff2

      ! indexes
      integer :: i, j, igv, nst, ireal

      ! allocate arrays based on number of data and max search
      allocate (rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch))
      allocate (nuse(ngvarg, ndata), useidx(ndata, ndata))
      allocate (sim(ngvarg + 1, ndata), isim(ngvarg, ndata), randpath(ndata))
      allocate (results(ndata), anisxyz(3, ndata))
      allocate (cmeans(ngvarg), cstdevs(ngvarg))

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
      do ireal = 1, nreals

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
         do i = 1, ndata

            ! get random simulation index
            simidx = randpath(i)

            !
            ! loop over factors at this simidx
            !
            do igv = 1, ngvarg

               ! anisotropic coord for this simidx
               axyz = matmul(pool(igv)%rm(:, :, nst), xyz(:, simidx))

               ! query the tree for this location
               call kdtree2_r_nearest(tp=trees(igv)%tree, qv=axyz, &
                                      r2=pool(igv)%aa(nst)**2, nfound=nfound, &
                                      nalloc=ndata, results=results)
               !
               ! loop over samples found in search
               !
               nuse(igv, simidx) = 0
               do j = 1, nfound

                  ! check if this data index is the simulation index
                  if (results(j)%idx .eq. simidx) cycle

                  ! check if this data index is simulated or not
                  if (isim(igv, results(j)%idx) .eq. 0) cycle ! no conditioning value here

                  ! meet minimum covariance? (not collocated)
                  rhs(1) = get_cov(pool(igv), xyz(:, simidx), xyz(:, results(j)%idx))
                  if (rhs(1) .lt. MINCOV) cycle

                  ! if we got this far increment number found at ith data location
                  nuse(igv, simidx) = nuse(igv, simidx) + 1

                  ! track conditioning indices found at ith data location
                  useidx(nuse(igv, simidx), simidx) = results(j)%idx

                  ! have we met the max search?
                  if (nuse(igv, simidx) .ge. nsearch) exit

               end do

               ! build and solve normal equations
               if (nuse(igv, simidx) .gt. 0) then
                  call krige(pool(igv), xyz, rhs, lhs, kwts, nuse(igv, :), useidx, &
                             sim(igv, :), simidx, cmeans(igv), cstdevs(igv))
               else
                  ! if no data the distribution is N(0,1)
                  cmeans(igv) = 0.d0
                  cstdevs(igv) = 1.d0
               end if

               ! update this location and factor as simulated
               isim(igv, simidx) = 1

            end do

            ! now that we have the conditional moments for each factor
            ! simulate values and check against the first tolerance
            !
            ! coarse search
            !
            diff1 = 999.0
            do while (diff1 .gt. tol1)

               ! simulate each factor
               do igv = 1, ngvarg
                  p = grnd()
                  call gauinv(p, xp, ierr)
                  sim(igv, simidx) = xp*cstdevs(igv) + cmeans(igv)
               end do

               ! simulate nugget
               p = grnd()
               call gauinv(p, xp, ierr)
               sim(ngvarg + 1, simidx) = xp

               ! calculate imputed value
               call network_forward(nnet, sim(simidx, :), zimp, .false.)

            end do

            !
            !
            !

         end do

      end do

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

end module impute_mod
