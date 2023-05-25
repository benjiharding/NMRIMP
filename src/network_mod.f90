module network_mod

   use geostat, only: nnet, wts
   use types_mod, only: network
   use subs, only: nscore

   implicit none

   ! interface for procedural pointer
   abstract interface
      function afunc(y) result(a)
         real(8), intent(in) :: y(:, :)
         real(8) :: a(size(y, 1), size(y, 2))
      end function afunc
   end interface

contains

   subroutine init_network(net)

      ! initialize network layer parameters and weight/bias matrices

      ! parameters
      type(network), intent(inout) :: net

      ! local variables
      integer, allocatable :: nwts(:), nbias(:)
      integer :: i, j

      ! allocate some counters matrix indices
      allocate (nwts(net%nl - 1), nbias(net%nl - 1))
      allocate (net%iwts(net%nl), net%ibias(net%nl))

      ! number of weights
      do i = 1, net%nl - 1
         nwts(i) = net%ld(i)*net%ld(i + 1)
      end do

      ! number of bias terms
      do j = 2, net%nl
         nbias(j - 1) = net%ld(j)
      end do

      ! weight matrix indices from cumulative sums
      net%iwts(1) = 0
      net%iwts(2:) = nwts
      do i = 2, net%nl
         net%iwts(i) = net%iwts(i - 1) + nwts(i - 1)
      end do

      ! bias matrix indices from cumulative sums
      net%ibias(1) = sum(nwts)
      net%ibias(2:) = nbias
      do j = 2, net%nl
         net%ibias(j) = net%ibias(j - 1) + nbias(j - 1)
      end do

      ! allocate network weight and bias matrices
      allocate (net%layer(net%nl - 1)) ! excludes input layer
      do i = 1, net%nl - 1
         ! get matrix shapes
         net%layer(i)%sw = [net%ld(i + 1), net%ld(i)] ! weight matrix shape
         net%layer(i)%sb = [net%ld(i + 1), 1] ! bias vector shape
         ! allocate
         allocate (net%layer(i)%nnwts(net%layer(i)%sw(1), net%layer(i)%sw(2)))
         allocate (net%layer(i)%nnbias(net%layer(i)%sb(1), net%layer(i)%sb(1)))
      end do

      ! total number of dimensions
      net%dims = sum(nwts) + sum(nbias)

   end subroutine init_network

   subroutine network_forward(net, Ymat, AL, nstrans)

      ! forward pass through network

      ! iterate over layers(2:)
      ! get number of connection in each layer (n^l * n^l-1)
      ! get cumulative sum of connecitons to get 1D ids, ie for each matrix
      ! do the same for the bias terms
      ! on each iteration reshape the array slice to (n^l * n^l-1)
      ! do matrix math and apply activations
      ! linear activation on final layer followed by nscore

      ! parameters
      type(network), intent(inout) :: net ! neural network object
      real(8), intent(in) :: Ymat(:, :) ! simulated factors
      logical, intent(in) :: nstrans ! nscore transform flag

      ! return
      real(8), intent(inout) :: AL(:) ! output mixture vector

      ! internal variables
      procedure(afunc), pointer :: f_ptr => null()
      real(8), allocatable :: Amat(:, :), A_prev(:, :)
      real(8), allocatable :: W(:, :), WL(:, :), b(:, :), bL(:, :), &
                              Zmat(:, :), ZL(:, :)
      real(8), allocatable :: vrg(:), tmp(:)
      integer :: i, ierr

      ! initialize the activation matrix
      Amat = Ymat

      ! pointer to activation function
      select case (net%af)
      case (1)
         f_ptr => sigmoid
      case (2)
         f_ptr => hyptan
      case (3)
         f_ptr => relu
      case (4)
         f_ptr => linear
      case (5)
         f_ptr => silu
      case (6)
         f_ptr => gelu
      end select

      ! hidden layers
      do i = 1, net%nl - 2 ! excludes input and output

         A_prev = Amat

         ! transpose prior to forward pass
         W = transpose(net%layer(i)%nnwts)
         b = transpose(net%layer(i)%nnbias)

         ! forward pass
         b = spread(b(1, :), 1, size(A_prev, dim=1))
         Zmat = matmul(A_prev, W) + b
         Amat = f_ptr(Zmat)

      end do

      ! output layer
      WL = transpose(net%layer(net%nl - 1)%nnwts)
      bL = transpose(net%layer(net%nl - 1)%nnbias)
      bL = spread(bL(1, :), 1, size(Amat, dim=1))
      ZL = matmul(Amat, WL) + bL

      ! linear activation and reduce dims
      AL = ZL(:, 1)

      ! normal score transform if required
      if (nstrans) then
         call nscore(size(AL), AL, dble(-1.0e21), dble(1.0e21), 1, &
                     wts, tmp, vrg, ierr)
         if (ierr .gt. 0) stop "Error in normal score transform"
         AL = vrg
      end if

   end subroutine network_forward

   subroutine vector_to_matrices(vector, net)

      ! reshape trial vector (DE output) to neural network weight matrices
      ! this subroutine updates the inupt type(network) object

      type(network), intent(inout) :: net
      real(8), intent(in) :: vector(:)
      integer :: i

      do i = 1, net%nl - 1

         ! reshape weights and biases
         net%layer(i)%nnwts = reshape(vector(net%iwts(i) + 1:net%iwts(i + 1)), &
                                      shape=(net%layer(i)%sw), order=[2, 1])
         net%layer(i)%nnbias = reshape(vector(net%ibias(i) + 1:net%ibias(i + 1)), &
                                       shape=(net%layer(i)%sb), order=[2, 1])
      end do

   end subroutine vector_to_matrices

   subroutine calc_regularization(net, reg)

      type(network), intent(inout) :: net
      real(8), intent(out) :: reg
      integer :: i

      reg = 0.d0

      do i = 1, net%nl - 1

         if (net%ireg .eq. 0) then
            reg = reg + 0.d0

         else if (net%ireg .eq. 1) then ! L1
            reg = reg + sum(abs(net%layer(i)%nnwts))*net%regconst

         else if (net%ireg .eq. 2) then ! L2
            reg = reg + sum(net%layer(i)%nnwts**2)*net%regconst
         end if

      end do

   end subroutine calc_regularization

   function relu(yval) result(a)

      ! rectified lienar unit activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = max(0.d0, yval)

   end function relu

   function silu(yval) result(a)

      ! sigmoid lienar unit activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = yval*sigmoid(yval)

   end function silu

   function gelu(yval) result(a)

      ! gaussian lienar unit activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = yval*sigmoid(1.702*yval)

   end function gelu

   function sigmoid(yval) result(a)

      ! sigmoid activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = 1.d0/(1.d0 + exp(-yval))

   end function sigmoid

   function hyptan(yval) result(a)

      ! hyperbolic tangent activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = dble(tanh(yval))

   end function hyptan

   function linear(yval) result(a)

      ! linear activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = yval

   end function linear

end module network_mod
