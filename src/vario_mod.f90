module vario_mod

   use covasubs, only: get_cov
   use rotationmatrix, only: set_rmat
   use types_mod, only: variogram

   implicit none

contains

   subroutine indicator_transform(zval, zc, nd, ncut, iz, ivars)

      ! indicator transform of zval based on cutoffs zc

      ! parameters
      real(8), intent(in) :: zval(:), zc(:)
      integer, intent(in) :: nd, ncut

      ! result
      integer, intent(inout) :: iz(:, :)
      real(8), intent(inout) :: ivars(:)

      ! local variables
      integer :: i, j
      real(8) :: prop

      ! initalize indicators
      iz = 0

      ! iterate over cutoffs and set indicators
      do i = 1, nd
         do j = 1, ncut
            if (zval(i) .le. zc(j)) then
               iz(i, j) = 1
            end if
         end do
      end do

      ! indicator variance to scale sill
      do j = 1, ncut
         prop = 0.d0
         do i = 1, nd
            prop = prop + iz(i, j)
         end do
         prop = prop/nd
         ivars(j) = prop*(1 - prop)
      end do

   end subroutine indicator_transform

   subroutine calc_expsill(var, sill)

      ! sill of traditional variogram (variance)

      real(8), intent(in) :: var(:)
      real(8), intent(out) :: sill
      real(8) :: mean, sumsqs
      integer :: i, nd

      nd = size(var)
      mean = 0.d0
      sumsqs = 0.d0

      do i = 1, nd
         mean = mean + var(i)
         sumsqs = sumsqs + var(i)*var(i)
      end do

      mean = mean/nd
      sumsqs = sumsqs/nd
      sill = sumsqs - mean*mean

   end subroutine calc_expsill

   subroutine gammabar(vm, xsiz, ysiz, zsiz, nx, ny, nz, gamb)

      !
      ! average variogram value of model vm in volume (xsiz, ysiz, zsiz)
      ! discretized by (nx, ny, nz)
      !

      ! paramerers
      type(variogram), intent(in) :: vm
      real(8), intent(in) :: xsiz, ysiz, zsiz
      integer, intent(in) :: nx, ny, nz
      real(8), intent(out) :: gamb

      ! locals
      real(8) :: maxcov, cov
      real(8) :: xmn, ymn, zmn
      real(8) :: xsz, ysz, zsz
      real(8) :: xzero, yzero, zzero
      integer :: ix, iy, iz
      integer :: jx, jy, jz
      real(8) :: xxi, yyi, zzi
      real(8) :: xxj, yyj, zzj

      ! discretization setup
      xsz = xsiz/real(nx)
      ysz = ysiz/real(ny)
      zsz = zsiz/real(nz)
      xmn = xsz/2.d0
      ymn = ysz/2.d0
      zmn = zsz/2.d0
      xzero = xsz*0.0001
      yzero = ysz*0.0001
      zzero = zsz*0.0001

      ! maximum covariance (variance)
      maxcov = get_cov(vm, [0.d0, 0.d0, 0.d0], [0.d0, 0.d0, 0.d0])
      gamb = 0.d0

      ! main loop over the volume
      do ix = 1, nx
         xxi = xmn + real(ix - 1)*xsz + xzero
         do iy = 1, ny
            yyi = ymn + real(iy - 1)*ysz + yzero
            do iz = 1, nz
               zzi = zmn + real(iz - 1)*zsz + zzero
               ! loop over pairs
               do jx = 1, nx
                  xxj = xmn + real(jx - 1)*xsz
                  do jy = 1, ny
                     yyj = ymn + real(jy - 1)*ysz
                     do jz = 1, nz
                        zzj = zmn + real(jz - 1)*zsz
                        ! calculate covariance for this pair
                        cov = get_cov(vm, [xxi, yyi, zzi], [xxj, yyj, zzj])
                        gamb = gamb + (maxcov - cov)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! final measure
      gamb = gamb/(real(nx*ny*nz)**2)

   end subroutine gammabar

   subroutine set_sill(vm)

      implicit none

      type(variogram), intent(inout) :: vm(:)
      integer :: n, i, j

      n = ubound(vm(:), 1)

      do i = 1, n
         vm(i)%sill = vm(i)%c0
         do j = 1, vm(i)%nst
            if (vm(i)%it(j) .eq. 4) then
               vm(i)%sill = vm(i)%sill + 999.D+00
            else
               vm(i)%sill = vm(i)%sill + vm(i)%cc(j)
            end if
         end do
      end do

   end subroutine set_sill

   subroutine set_rotmatrix(vm)

      implicit none

      type(variogram), intent(inout) :: vm(:)

      integer :: n, i, j, test

      n = ubound(vm(:), 1)

      do i = 1, n
         if (.not. (allocated(vm(i)%rm))) then
            allocate (vm(i)%rm(3, 3, vm(i)%nst), stat=test)
            if (test .ne. 0) stop 'ERROR: Allocation failed due to insufficient memory.'
         end if
         do j = 1, vm(i)%nst
            vm(i)%rm(:, :, j) = set_rmat([vm(i)%ang1(j), vm(i)%ang2(j), vm(i)%ang3(j)], &
                                         [1.0D+00, vm(i)%anis1(j), vm(i)%anis2(j)])
         end do
      end do

   end subroutine set_rotmatrix

end module vario_mod
