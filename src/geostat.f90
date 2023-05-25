module geostat

   use types_mod

   implicit none

   !
   ! GLOBAL VARIABLE DECLARATIONS
   !

   ! flags
   integer :: itrans
   integer :: idbg

   ! network variables
   type(network) :: nnet

   ! simulation variables
   integer :: nsearch
   integer :: nreals
   integer :: rseed
   real(8) :: radius, radius1, radius2
   real(8) :: sang1, sang2, sang3
   real(8) :: sanis1, sanis2
   real(8) :: srotmat(3, 3)

   ! data variables
   real(8), allocatable :: var(:)
   real(8), allocatable :: nsvar(:)
   real(8), allocatable :: wts(:)

   ! data parameters
   real(8), allocatable :: xyz(:, :) ! data coordinates
   integer, allocatable :: dhids(:)
   integer :: ndata

   ! Gaussian pool parameters
   integer :: ngvarg
   type(variogram), allocatable :: pool(:)

contains

end module geostat
