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
   real(8), allocatable :: zinit(:, :) ! (ndata, nfact)
   integer, allocatable :: isim(:, :) ! (ndata, nfact)
   real(8), allocatable :: sim(:, :) ! (ndata, nfact)
   real(8), allocatable :: anisxyz(:, :)
   integer, allocatable :: randpath(:)
   integer :: nsearch
   integer :: nreals
   integer :: nfact
   integer :: nst
   integer :: rseed
   integer :: iter1, iter2
   real(8) :: tol1, tol2

   ! normal equations
   real(8), allocatable :: rhs(:), lhs(:, :), kwts(:)
   real(8), allocatable :: cmeans(:), cstdevs(:)
   integer, allocatable :: nuse(:, :), useidx(:, :)

   ! reference distribution
   integer, parameter :: nsamp = 1000000
   integer :: nlook
   real(8) :: zref(nsamp)
   real(8), allocatable :: nsref(:), yref(:, :)
   real(8) :: usnsref(nsamp)
   real(8) :: minref, maxref
   real(8) :: gmin, gmax

   ! dynamic tolerances
   real(8), allocatable :: vsort(:), wsort(:), vord(:)
   real(8), allocatable :: vcdf(:)

   ! kdtree required variables
   type(kdtrees), allocatable :: trees(:) ! array of pointers
   type(kdtree2_result), allocatable :: results(:)

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
