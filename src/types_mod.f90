module types_mod

   use kdtree2_module

   implicit none

   ! variogram model
   type variogram
      integer :: nst ! number of nested structure
      real(8) :: c0 ! nugget
      real(8) :: sill
      integer, allocatable :: it(:) ! structure types
      real(8), allocatable :: cc(:) ! structure contributions
      real(8), allocatable :: ang1(:), ang2(:), ang3(:) ! structure angles
      real(8), allocatable :: aa(:), anis1(:), anis2(:) ! structure anisotropy
      real(8), allocatable :: ahmin(:), avert(:) ! minor and vert ranges
      real(8), allocatable :: rm(:, :, :) ! rotation matrix
   end type variogram

   ! neural network parameters
   type weights
      real(8), allocatable :: nnwts(:, :) ! (n, n)
      real(8), allocatable :: nnbias(:, :) ! (n, 1)
      real(8), allocatable :: nnmu(:) ! normalization mean (1, n)
      real(8), allocatable :: nnsig(:) ! normalization stdev (1, n)
      real(8), allocatable :: gmma(:) ! normalization loc (1, n)
      real(8), allocatable :: beta(:) ! normalization scale (1, n)
      integer :: sw(2), sb(2) ! shape of weights and bias matrices
   end type weights

   type network
      integer :: nl ! number of layers
      integer, allocatable :: ld(:) ! layer dimensions
      type(weights), allocatable :: layer(:) ! layer matrices (n, n)
      integer, allocatable :: iwts(:), ibias(:) ! vector indices
      integer, allocatable :: igmma(:), ibeta(:) ! vector indices
      integer :: af ! activation function
      integer :: dims ! total number of weights + biases
      integer :: ireg ! regularizer (0=none, 1=L1, 2=L2)
      real(8) :: regconst ! regularization constant
      logical :: norm ! batch normalization?
   end type network

   ! kdtree pointers
   type kdtrees
      type(kdtree2), pointer :: tree
   end type kdtrees

end module types_mod
