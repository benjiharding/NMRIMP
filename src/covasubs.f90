module covasubs

   use types_mod, only: variogram
   use constants, only: PI

   implicit none

contains

!Calculate the covariance between two points
   real*8 function get_cov(vm, a, b) result(cov)
      real*8, parameter :: EPSLON = 1.0e-5
      type(variogram), intent(in) :: vm
      real*8, dimension(3), intent(in) :: a, b
      integer :: i
      real*8 :: h, hr, hsqd

! Check for "zero" distance, return with cmax if so:
      if (dsqrd(a, b) < EPSLON) then
         cov = vm%sill; return
      end if

! Loop over all the structures:
      cov = 0.0D+00
      do i = 1, vm%nst
         hsqd = dsqrd(a, b, vm%rm(:, :, i))
         h = dsqrt(hsqd)
         hr = h/vm%aa(i)
         VARIOTYPE:select case(vm%it(i))
         case (1)   ! Spherical Variogram Model?
         if (hr < 1.0D+00) cov = cov + vm%cc(i)*(1.0D+00 - hr*(1.5D+00 - 0.5D+00*hr**2))
         case (2) ! Exponential Variogram Model?
         cov = cov + vm%cc(i)*exp(-3.0D+00*hr)
         case (3) ! Gaussian Variogram Model?
         cov = cov + vm%cc(i)*exp(-3.0D+00*hr**2)
         case (4) ! Power Variogram Model?
         cov = cov + vm%sill - vm%cc(i)*(h**vm%aa(i))
         case (5) ! Hole Effect Model?
         cov = cov + vm%cc(i)*dcos(hr*PI)
         end select VARIOTYPE
      end do
   end function get_cov

!Get the anisotropic distance between two points
   real*8 function dsqrd(a, b, rm)
      real*8, intent(in) :: a(3), b(3)
      real*8, optional :: rm(3, 3)
      real*8 :: ba(3)
      ba = b - a
      if (present(rm)) then
         dsqrd = (rm(1, 1)*ba(1) + rm(1, 2)*ba(2) + rm(1, 3)*ba(3))**2 + &
                 (rm(2, 1)*ba(1) + rm(2, 2)*ba(2) + rm(2, 3)*ba(3))**2 + &
                 (rm(3, 1)*ba(1) + rm(3, 2)*ba(2) + rm(3, 3)*ba(3))**2
      else
         dsqrd = ba(1)**2 + ba(2)**2 + ba(3)**2
      end if
   end function dsqrd

   subroutine cova3(x1, y1, z1, x2, y2, z2, ivarg, nst, MAXNST, c0, it, cc, aa, &
                    irot, MAXROT, rotmat, cmax, cova)
!-----------------------------------------------------------------------
!
!                    Covariance Between Two Points
!                    *****************************
!
! This subroutine calculated the covariance associated with a variogram
! model specified by a nugget effect and nested varigoram structures.
! The anisotropy definition can be different for each nested structure.
!
!
!
! INPUT VARIABLES:
!
!   x1,y1,z1         coordinates of first point
!   x2,y2,z2         coordinates of second point
!   nst(ivarg)       number of nested structures (maximum of 4)
!   ivarg            variogram number (set to 1 unless doing cokriging
!                       or indicator kriging)
!   MAXNST           size of variogram parameter arrays
!   c0(ivarg)        isotropic nugget constant
!   it(i)            type of each nested structure:
!                      1. spherical model of range a;
!                      2. exponential model of parameter a;
!                           i.e. practical range is 3a
!                      3. gaussian model of parameter a;
!                           i.e. practical range is a*sqrt(3)
!                      4. power model of power a (a must be gt. 0  and
!                           lt. 2).  if linear model, a=1,!=slope.
!                      5. hole effect model
!   cc(i)            multiplicative factor of each nested structure.
!                      (sill-c0) for spherical, exponential,and gaussian
!                      slope for linear model.
!   aa(i)            parameter "a" of each nested structure.
!   irot             index of the rotation matrix for the first nested
!                    structure (the second nested structure will use
!                    irot+1, the third irot+2, and so on)
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!
!
! OUTPUT VARIABLES:
!
!   cmax             maximum covariance
!   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
!
!
!
! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                      rotmat    computes rotation matrix for distance
!-----------------------------------------------------------------------
      real(8), parameter :: PI = 3.14159265D0, PMX = 999.D0, EPSLON = 1.d-5
      integer :: nst, it(*), istart, ivarg, MAXNST, MAXROT, irot, is, ist, ir
      real(8) :: c0, cc(*), aa(*), cova, cmax, hr, h
      real(8) :: rotmat(MAXROT, 3, 3), hsqd, x1, y1, z1, x2, y2, z2

! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):

      istart = 1 + (ivarg - 1)*MAXNST
      cmax = c0
      do is = 1, nst
         ist = istart + is - 1
         if (it(ist) .eq. 4) then
            cmax = cmax + PMX
         else
            cmax = cmax + cc(ist)
         end if
      end do

! Check for "zero" distance, return with cmax if so:

      hsqd = sqdist(x1, y1, z1, x2, y2, z2, irot, MAXROT, rotmat)
      if (dble(hsqd) .lt. EPSLON) then
         cova = cmax
         return
      end if

! Loop over all the structures:

      cova = 0.D0
      do is = 1, nst
         ist = istart + is - 1

         ! Compute the appropriate distance:

         if (ist .ne. 1) then
            ir = min((irot + is - 1), MAXROT)
            hsqd = sqdist(x1, y1, z1, x2, y2, z2, ir, MAXROT, rotmat)
         end if
         h = dble(dsqrt(hsqd))

         ! Spherical Variogram Model?

         if (it(ist) .eq. 1) then
            hr = h/aa(ist)
            if (hr .lt. 1.D0) cova = cova + cc(ist)*(1.D0 - hr*(1.5 - .5*hr*hr))

            ! Exponential Variogram Model?

         else if (it(ist) .eq. 2) then
            cova = cova + cc(ist)*exp(-3.D0*h/aa(ist))

            ! Gaussian Variogram Model?

         else if (it(ist) .eq. 3) then
            cova = cova + cc(ist)*exp(-3.D0*(h/aa(ist))*(h/aa(ist)))

            ! Power Variogram Model?

         else if (it(ist) .eq. 4) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))

            ! Hole Effect Model?

         else if (it(ist) .eq. 5) then
            !                 d = 10.0 * aa(ist)
            !                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
            cova = cova + cc(ist)*cos(h/aa(ist)*PI)
         end if
      end do

   end subroutine cova3

   subroutine setrot(ang1, ang2, ang3, anis1, anis2, ind, MAXROT, rotmat)
      !-----------------------------------------------------------------------
      !
      !              Sets up an Anisotropic Rotation Matrix
      !              **************************************
      !
      ! Sets up the matrix to transform cartesian coordinates to coordinates
      ! accounting for angles and anisotropy (see manual for a detailed
      ! definition):
      !
      !
      ! INPUT PARAMETERS:
      !
      !   ang1             Azimuth angle for principal direction
      !   ang2             Dip angle for principal direction
      !   ang3             Third rotation angle
      !   anis1            First anisotropy ratio
      !   anis2            Second anisotropy ratio
      !   ind              matrix indicator to initialize
      !   MAXROT           maximum number of rotation matrices dimensioned
      !   rotmat           rotation matrices
      !
      !
      ! NO EXTERNAL REFERENCES
      !
      !
      !-----------------------------------------------------------------------
      integer MAXROT, ind
      real(8), parameter :: DEG2RAD = 3.141592654/180.D0, EPSLON = 1.d-20
      real(8) :: rotmat(MAXROT, 3, 3), afac1, afac2, sina, sinb, sint, &
                 cosa, cosb, cost, ang1, ang2, ang3, anis1, anis2, &
                 alpha, beta, theta, sin, cos

      !
      ! Converts the input angles to three angles which make more
      !  mathematical sense:
      !
      !         alpha   angle between the major axis of anisotropy and the
      !                 E-W axis. Note: Counter clockwise is positive.
      !         beta    angle between major axis and the horizontal plane.
      !                 (The dip of the ellipsoid measured positive down)
      !         theta   Angle of rotation of minor axis about the major axis
      !                 of the ellipsoid.
      !
      if (ang1 .ge. 0.D0 .and. ang1 .lt. 270.D0) then
         alpha = (90.D0 - ang1)*DEG2RAD
      else
         alpha = (450.D0 - ang1)*DEG2RAD
      end if
      beta = -1.D0*ang2*DEG2RAD
      theta = ang3*DEG2RAD
      !
      ! Get the required sines and cosines:
      !
      sina = dble(sin(alpha))
      sinb = dble(sin(beta))
      sint = dble(sin(theta))
      cosa = dble(cos(alpha))
      cosb = dble(cos(beta))
      cost = dble(cos(theta))
      !
      ! Construct the rotation matrix in the required memory:
      !
      afac1 = 1.D0/dble(max(anis1, EPSLON))
      afac2 = 1.D0/dble(max(anis2, EPSLON))
      rotmat(ind, 1, 1) = (cosb*cosa)
      rotmat(ind, 1, 2) = (cosb*sina)
      rotmat(ind, 1, 3) = (-sinb)
      rotmat(ind, 2, 1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind, 2, 2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind, 2, 3) = afac1*(sint*cosb)
      rotmat(ind, 3, 1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind, 3, 2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind, 3, 3) = afac2*(cost*cosb)

   end subroutine setrot

   subroutine setrotmat(azm, dip, tilt, forwardrotmat, reverserotmat)
      !-----------------------------------------------------------------------
      !
      !        Sets up the forward rotation matrix and reverse matrix
      !
      ! (c) Jared Deutsch, 2014
      !-----------------------------------------------------------------------
      implicit none
      real(kind=8), parameter :: PI = 4*atan(1.0d0)
      real(kind=8), intent(in) :: azm, dip, tilt
      real(kind=8), dimension(3, 3) :: azmrotmat, diprotmat, tiltrotmat
      real(kind=8), dimension(3, 3), intent(out) :: forwardrotmat, reverserotmat
      real(kind=8) :: angle

      ! Get the angle to rotate about the Z axis in radians
      angle = (90d0 - azm)*PI/180d0
      ! Rotation matrix for azimuth correction
      azmrotmat = 0d0
      azmrotmat(3, 3) = 1d0
      azmrotmat(1, 1) = cos(angle)
      azmrotmat(1, 2) = -1d0*sin(angle)
      azmrotmat(2, 1) = sin(angle)
      azmrotmat(2, 2) = cos(angle)

      ! Get the angle to rotate about the new X' axis
      angle = -1d0*(dip)*PI/180d0
      ! Rotation matrix for dip correction
      diprotmat = 0d0
      diprotmat(2, 2) = 1d0
      diprotmat(1, 1) = cos(angle)
      diprotmat(1, 3) = sin(angle)
      diprotmat(3, 1) = -1d0*sin(angle)
      diprotmat(3, 3) = cos(angle)

      ! Get the angle to rotate about the new Y' axis
      angle = -1d0*(tilt)*PI/180d0
      ! Rotation matrix for tilt correction
      tiltrotmat = 0d0
      tiltrotmat(1, 1) = 1d0
      tiltrotmat(2, 2) = cos(angle)
      tiltrotmat(2, 3) = sin(angle)
      tiltrotmat(3, 2) = -1d0*sin(angle)
      tiltrotmat(3, 3) = cos(angle)

      ! Complete forward rotation matrix is the product of these 2 matrices
      forwardrotmat = matmul(matmul(azmrotmat, diprotmat), tiltrotmat)
      ! Reverse rotation matrix is the transpose of the forward matrix
      ! as these matrices are orthogonal
      reverserotmat = transpose(forwardrotmat)
   end subroutine setrotmat

   real(kind=8) function distance(xyz_a, xyz_b)
      !-----------------------------------------------------------------------
      !
      !                  Distance given 2 x,y,z vectors
      !
      ! (c) Jared Deutsch, 2014
      !-----------------------------------------------------------------------
      real(kind=8), dimension(3), intent(in) :: xyz_a, xyz_b
      distance = sqrt((xyz_a(1) - xyz_b(1))**2 + &
                      (xyz_a(2) - xyz_b(2))**2 + &
                      (xyz_a(3) - xyz_b(3))**2)
   end function distance

   function sqdist(x1, y1, z1, x2, y2, z2, ind, MAXROT, rotmat) result(hsqd)
      !-----------------------------------------------------------------------
      !
      !    Squared Anisotropic Distance Calculation Given Matrix Indicator
      !    ***************************************************************
      !
      ! This routine calculates the anisotropic distance between two points
      !  given the coordinates of each point and a definition of the
      !  anisotropy.
      !
      !
      ! INPUT VARIABLES:
      !
      !   x1,y1,z1         Coordinates of first point
      !   x2,y2,z2         Coordinates of second point
      !   ind              The rotation matrix to use
      !   MAXROT           The maximum number of rotation matrices dimensioned
      !   rotmat           The rotation matrices
      !
      !
      !
      ! OUTPUT VARIABLES:
      !
      !   sqdist           The squared distance accounting for the anisotropy
      !                      and the rotation of coordinates (if any).
      !
      !
      ! NO EXTERNAL REFERENCES
      !
      !
      !-----------------------------------------------------------------------
      integer :: i, ind, MAXROT
      real(8) :: rotmat(MAXROT, 3, 3), cont, dx, dy, dz, x1, y1, z1, x2, y2, z2, &
                 hsqd

      !
      ! Compute component distance vectors and the squared distance:
      !
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      hsqd = 0.D0
      do i = 1, 3
         cont = rotmat(ind, i, 1)*dx &
                + rotmat(ind, i, 2)*dy &
                + rotmat(ind, i, 3)*dz
         hsqd = hsqd + cont*cont
      end do

   end function sqdist

end module covasubs
