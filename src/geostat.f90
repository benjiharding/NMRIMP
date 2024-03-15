module geostat

   use types_mod

   implicit none

   !
   ! GLOBAL VARIABLE DECLARATIONS
   !

   ! flags
   integer :: itrans
   integer :: idbg
   integer :: isort = 0 ! sort factors by gammabar?

   ! network variables
   type(network) :: nnet

   ! simulation variables
   real(8), allocatable :: imputed(:, :, :) ! (ndata, nfact, nreals)
   real(8), allocatable :: zinit(:, :) ! (ndata, nfact)
   integer, allocatable :: isim(:, :) ! (ndata, nfact)
   real(8), allocatable :: sim(:, :) ! (ndata, nfact)
   real(8), allocatable :: anisxyz(:, :)
   real(8), allocatable :: abratio(:, :) ! (ndata, nreals)
   integer, allocatable :: rp1(:), rp2(:), rp1s(:), randpath(:)
   real(8), allocatable :: rp1_r(:) ! rand path 1 as reals for sorting
   integer, allocatable :: rsc(:, :) ! resim count
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
   real(8) :: zref(nsamp)
   real(8), allocatable :: nsref(:), yref(:, :)
   real(8) :: usnsref(nsamp)
   real(8) :: minref, maxref
   real(8) :: gmin, gmax

   ! kdtree required variables
   type(kdtrees), allocatable :: trees(:) ! array of pointers
   type(kdtree2_result), allocatable :: results(:)

   ! gammabar for sorting by continuity
   real(8) :: gb
   real(8), allocatable :: gbar(:), gbaridx(:)
   integer :: nsort

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

   ! factor precedence
   integer, allocatable :: fprec(:)
   real(8), allocatable :: fprec_r(:), fidx(:)
   real(8), allocatable :: sigwt(:) ! sigmoid weighting factor
   real(8) :: qfp, tfp ! quantile and threshold for factor with preference
   real(8), allocatable :: ub(:), lb(:)
   real(8) :: sr ! exclution radius for seeding
   integer, allocatable :: seeded(:, :) ! (ndata, nreals)
   logical :: isd, ifp

   integer :: fpid, sdid
   real(8) :: fpwt

   ! output file
   integer :: lout, ldbg, lmom

contains

end module geostat
