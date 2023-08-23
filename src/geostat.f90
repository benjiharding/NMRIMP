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
   real(8), allocatable :: imputed(:, :, :) ! (ndata, nfact, nreals)
   real(8), allocatable :: zinit(:, :)
   integer :: nsearch
   integer :: nreals
   integer :: rseed
   integer :: iter1, iter2
   real(8) :: tol1, tol2

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

   ! output file
   integer :: lout, ldbg, lmom

contains

end module geostat
