module network_mod

   use geostat
   use types_mod, only: network
   use subs, only: nscore, powint, locate, sortem
   use mtmod, only: grnd
   use constants

   implicit none

   ! interface for procedural pointer
   abstract interface
      function afunc(y) result(a)
         real(8), intent(in) :: y(:, :)
         real(8) :: a(size(y, 1), size(y, 2))
      end function afunc
   end interface

contains

   ! subroutine init_network(net)

   !    ! initialize network layer parameters and weight/bias matrices

   !    ! parameters
   !    type(network), intent(inout) :: net

   !    ! local variables
   !    integer, allocatable :: nwts(:), nbias(:)
   !    integer, allocatable :: ngmma(:), nbeta(:)
   !    integer :: i, j

   !    ! allocate some counters matrix indices
   !    allocate (nwts(net%nl - 1), nbias(net%nl - 1))
   !    allocate (ngmma(net%nl - 1), nbeta(net%nl - 1))
   !    allocate (net%iwts(net%nl), net%ibias(net%nl))
   !    allocate (net%igmma(net%nl), net%ibeta(net%nl))

   !    ! number of weights
   !    do i = 1, net%nl - 1
   !       nwts(i) = net%ld(i)*net%ld(i + 1)
   !    end do

   !    ! number of bias terms
   !    do j = 2, net%nl
   !       nbias(j - 1) = net%ld(j)
   !    end do

   !    ! trainable loc and scale for batch norm
   !    if (net%norm) then
   !       do j = 2, net%nl
   !          ngmma(j - 1) = net%ld(j)
   !          nbeta(j - 1) = net%ld(j)
   !       end do
   !    end if

   !    ! weight matrix indices from cumulative sums
   !    net%iwts(1) = 0
   !    net%iwts(2:) = nwts
   !    do i = 2, net%nl
   !       net%iwts(i) = net%iwts(i - 1) + nwts(i - 1)
   !    end do

   !    ! bias vector indices from cumulative sums
   !    net%ibias(1) = sum(nwts)
   !    net%ibias(2:) = nbias
   !    do j = 2, net%nl
   !       net%ibias(j) = net%ibias(j - 1) + nbias(j - 1)
   !    end do

   !    if (net%norm) then
   !       ! gamma vector indices from cumulative sums
   !       net%igmma(1) = sum(nbias) + sum(nwts)
   !       net%igmma(2:) = ngmma
   !       do j = 2, net%nl
   !          net%igmma(j) = net%igmma(j - 1) + ngmma(j - 1)
   !       end do

   !       ! beta vector indices from cumulative sums
   !       net%ibeta(1) = sum(nbias) + sum(nwts) + sum(ngmma)
   !       net%ibeta(2:) = nbeta
   !       do j = 2, net%nl
   !          net%ibeta(j) = net%ibeta(j - 1) + nbeta(j - 1)
   !       end do
   !    end if

   !    ! allocate network weight and bias matrices
   !    allocate (net%layer(net%nl - 1)) ! excludes input layer
   !    do i = 1, net%nl - 1

   !       ! get matrix shapes
   !       net%layer(i)%sw = [net%ld(i + 1), net%ld(i)] ! weight matrix shape
   !       net%layer(i)%sb = [net%ld(i + 1), 1] ! bias vector shape

   !       ! allocate matrices
   !       allocate (net%layer(i)%nnwts(net%layer(i)%sw(1), net%layer(i)%sw(2)))
   !       allocate (net%layer(i)%nnbias(net%layer(i)%sb(1), net%layer(i)%sb(2)))

   !       ! allocate moment vectors
   !       allocate (net%layer(i)%nnmu(net%layer(i)%sb(1)))
   !       allocate (net%layer(i)%nnsig(net%layer(i)%sb(1)))

   !       ! allocate loc and scale vectors
   !       allocate (net%layer(i)%gmma(net%layer(i)%sb(1)))
   !       allocate (net%layer(i)%beta(net%layer(i)%sb(1)))
   !    end do

   !    ! total number of dimensions
   !    if (net%norm) then
   !       net%dims = sum(nwts) + sum(nbias) + sum(ngmma) + sum(nbeta)
   !    else
   !       net%dims = sum(nwts) + sum(nbias)
   !    end if

   ! end subroutine init_network

   subroutine init_network(net)

      ! initialize network layer parameters and weight/bias matrices

      ! parameters
      type(network), intent(inout) :: net

      net%dims = 2*net%ld(1) ! 2L trainable parameters

      ! allocate weight annd exponent arrays
      allocate (net%awts(1, net%ld(1))) ! row vector
      allocate (net%omega(net%ld(1)))

   end subroutine init_network

   ! subroutine network_forward(net, Ymat, AL, nstrans, norm, calc_mom)

   !    ! forward pass through network with optional normalization

   !    ! iterate over layers(2:)
   !    ! get number of connection in each layer (n^l * n^l-1)
   !    ! get cumulative sum of connecitons to get 1D ids, ie for each matrix
   !    ! do the same for the bias terms
   !    ! on each iteration reshape the array slice to (n^l * n^l-1)
   !    ! do matrix math and apply activations
   !    ! linear activation on final layer followed by nscore

   !    ! parameters
   !    type(network), intent(inout) :: net ! neural network object
   !    real(8), intent(in) :: Ymat(:, :) ! simulated factors
   !    logical, intent(in) :: nstrans ! nscore transform flag
   !    logical, intent(in) :: norm ! normalize activations flag
   !    logical, intent(in) :: calc_mom ! calculate moments?

   !    ! return
   !    real(8), intent(inout) :: AL(:) ! output mixture vector

   !    ! internal variables
   !    procedure(afunc), pointer :: f_ptr => null()
   !    real(8), allocatable :: Amat(:, :), A_prev(:, :), Anorm(:, :)
   !    real(8), allocatable :: W(:, :), WL(:, :), b(:, :), bL(:, :), &
   !                            Zmat(:, :), Znorm(:, :), ZL(:, :), ZLnorm(:, :)
   !    real(8), allocatable :: vrg(:), tmp(:)
   !    integer :: i, ierr, idx

   !    ! initialize the activation matrix
   !    Amat = Ymat

   !    ! pointer to activation function
   !    select case (net%af)
   !    case (1)
   !       f_ptr => sigmoid
   !    case (2)
   !       f_ptr => hyptan
   !    case (3)
   !       f_ptr => relu
   !    case (4)
   !       f_ptr => linear
   !    case (5)
   !       f_ptr => silu
   !    case (6)
   !       f_ptr => gelu
   !    end select

   !    ! hidden layers
   !    do i = 1, net%nl - 2 ! excludes input and output

   !       A_prev = Amat

   !       ! transpose prior to forward pass
   !       W = transpose(net%layer(i)%nnwts)
   !       b = transpose(net%layer(i)%nnbias)

   !       ! forward pass - BN prior to activation
   !       b = spread(b(1, :), 1, size(A_prev, dim=1))
   !       Zmat = matmul(A_prev, W) + b
   !       if (norm) then
   !          if (calc_mom) then
   !             call normalize_input(Zmat, Znorm, net, i, calc_mom=.true., &
   !                                  gmma=net%layer(i)%gmma, beta=net%layer(i)%beta)
   !          else
   !             if (.not. allocated(net%layer(i)%nnmu)) stop "moment vectors must be allocated"
   !             call normalize_input(Zmat, Znorm, net, i, calc_mom=.false., &
   !                                  mu=net%layer(i)%nnmu, sig=net%layer(i)%nnsig, &
   !                                  gmma=net%layer(i)%gmma, beta=net%layer(i)%beta)
   !          end if
   !       else
   !          Znorm = Zmat
   !       end if
   !       Amat = f_ptr(Znorm)

   !       ! ! forward pass - activation prior to BN
   !       ! b = spread(b(1, :), 1, size(A_prev, dim=1))
   !       ! Zmat = matmul(A_prev, W) + b
   !       ! Amat = f_ptr(Zmat)
   !       ! if (norm) then
   !       !    if (calc_mom) then
   !       !       call normalize_input(Amat, Anorm, net, i, calc_mom=.true., &
   !       !                            gmma=net%layer(i)%gmma, beta=net%layer(i)%beta)
   !       !    else
   !       !       if (.not. allocated(net%layer(i)%nnmu)) stop "moment vectors must be allocated"
   !       !       call normalize_input(Amat, Anorm, net, i, calc_mom=.false., &
   !       !                            mu=net%layer(i)%nnmu, sig=net%layer(i)%nnsig, &
   !       !                            gmma=net%layer(i)%gmma, beta=net%layer(i)%beta)
   !       !    end if
   !       !    Amat = Anorm
   !       ! end if

   !    end do

   !    ! output layer
   !    WL = transpose(net%layer(net%nl - 1)%nnwts)
   !    bL = transpose(net%layer(net%nl - 1)%nnbias)
   !    bL = spread(bL(1, :), 1, size(Amat, dim=1))
   !    ZL = matmul(Amat, WL) + bL

   !    if (norm) then
   !       idx = net%nl - 1
   !       if (calc_mom) then
   !          call normalize_input(ZL, ZLnorm, net, idx, calc_mom=.true., &
   !                               gmma=net%layer(idx)%gmma, beta=net%layer(idx)%beta)
   !       else
   !          if (.not. allocated(net%layer(idx)%nnmu)) stop "moment vectors must be allocated"
   !          call normalize_input(ZL, ZLnorm, net, idx, calc_mom=.false., &
   !                               mu=net%layer(idx)%nnmu, sig=net%layer(idx)%nnsig, &
   !                               gmma=net%layer(idx)%gmma, beta=net%layer(idx)%beta)
   !       end if
   !       ZL = ZLnorm
   !    end if

   !    ! linear activation and reduce dims
   !    AL = ZL(:, 1)

   !    ! normal score transform if required
   !    if (nstrans) then
   !       do i = 1, size(AL)
   !          AL(i) = AL(i) + grnd()*SMALLDBLE ! random despike
   !       end do
   !       call nscore(size(AL), AL, dble(-1.0e21), dble(1.0e21), 1, &
   !                   wts, tmp, vrg, ierr)
   !       if (ierr .gt. 0) stop "Error in normal score transform"
   !       AL = vrg
   !    end if

   ! end subroutine network_forward

   subroutine network_forward(net, Ymat, AL, nstrans, ttable)

      ! parameters
      type(network), intent(inout) :: net ! neural network object
      real(8), intent(in) :: Ymat(:, :) ! simulated factors
      logical, intent(in) :: nstrans ! nscore transform flag
      real(8), optional :: ttable(:, :) ! transform table

      ! return
      real(8), intent(inout) :: AL(:) ! output mixture vector

      ! locals
      real(8), allocatable :: Amat(:, :), ZL(:, :)
      integer :: i

      ! forward pass; assumes single layer network
      Amat = Ymat
      do i = 1, net%ld(1)
         Amat(:, i) = power(Ymat(:, i), net%omega(i), 1.d0)
      end do
      ZL = matmul(Amat, transpose(net%awts))
      AL = ZL(:, 1)

      ! random despike and normal score transform
      if (nstrans) then
         do i = 1, ndata
            AL(i) = AL(i) + grnd()*SMALLDBLE
            call transform_to_refcdf(AL(i), ttable, AL(i), nsamp)
         end do
      end if

   end subroutine network_forward

   subroutine network_forward2(net, Ymat, AL, nstrans, fp, sigwt, ttable)

      ! parameters
      type(network), intent(inout) :: net ! neural network object
      real(8), intent(in) :: Ymat(:, :) ! simulated factors
      logical, intent(in) :: nstrans ! nscore transform flag
      integer, intent(in) :: fp ! factor precedence
      real(8), intent(in) :: sigwt ! sigmoid weight
      real(8), optional :: ttable(:, :) ! transform table

      ! return
      real(8), intent(inout) :: AL(:) ! output mixture vector

      ! locals
      real(8), allocatable :: fw(:)
      integer :: i

      ! weight factor with highest precedence
      fw = sigmoid1d(Ymat(:, fp)*sigwt)
      AL = net%awts(1, fp)*power(Ymat(:, fp), net%omega(fp), 1.d0)

      ! forward pass through remaining factors
      do i = 1, net%ld(1)
         if (i .eq. fp) cycle
         AL = AL + net%awts(1, i)*power(Ymat(:, i), net%omega(i), 1.d0)*fw
      end do

      ! random despike and normal score transform
      if (nstrans) then
         do i = 1, ndata
            AL(i) = AL(i) + grnd()*SMALLDBLE
            call transform_to_refcdf(AL(i), ttable, AL(i), nsamp)
         end do
      end if

   end subroutine network_forward2

   subroutine build_refcdf(nsamp, yref, net, ttable)

      !
      ! build reference CDF for normal score transform of
      ! activation values
      !
      integer, intent(in) :: nsamp ! num samples in ref CDF
      real(8), intent(in) :: yref(:, :) ! Gauss lookup table
      type(network), intent(inout) :: net ! network object with weights
      real(8), intent(out) :: ttable(nsamp, 2) ! transform table

      ! locals
      real(8) :: zref(nsamp) ! activations
      real(8), allocatable :: nsref(:) ! nscore activations
      real(8), allocatable :: tmp(:)
      real(8) :: wt(nsamp)
      integer :: i, ierr, nfact

      nfact = ngvarg + 1
      wt = 1.d0

      ! calculate the corresponding z values
      if (ifp) then
         call network_forward2(net, yref, zref, nstrans=.false., fp=fpid, &
                               sigwt=fpwt)
      else
         call network_forward(net, yref, zref, nstrans=.false.)
      end if

      ! normal score to build transform table
      do i = 1, nsamp
         zref(i) = zref(i) + grnd()*SMALLDBLE ! random despike
      end do
      call nscore(nsamp, zref, dble(-1e21), dble(1e21), 1, wt, &
                  tmp, nsref, ierr)

      ! sort the reference distribution (for interp.) and return
      call sortem(1, nsamp, zref, 1, nsref, zref, zref, zref, &
                  zref, zref, zref, zref)
      ttable(:, 1) = zref
      ttable(:, 2) = nsref

   end subroutine build_refcdf

   subroutine transform_to_refcdf(a, ttable, az, nsamp)

      real(8), intent(in) :: a ! activation value
      real(8), intent(in) :: ttable(:, :) ! transform table
      integer, intent(in) :: nsamp
      real(8), intent(out) :: az ! transformed normal score value

      ! locals
      integer :: j

      call locate(ttable(:, 1), nsamp, 1, nsamp, a, j)
      j = min(max(1, j), (nsamp - 1))
      az = powint(ttable(j, 1), ttable(j + 1, 1), ttable(j, 2), &
                  ttable(j + 1, 2), a, 1.d0)

   end subroutine transform_to_refcdf

   ! subroutine vector_to_matrices(vector, net)

   !    ! reshape trial vector (DE output) to neural network weight matrices
   !    ! this subroutine updates the inupt type(network) object

   !    type(network), intent(inout) :: net
   !    real(8), intent(in) :: vector(:)
   !    integer :: i

   !    do i = 1, net%nl - 1

   !       ! reshape weights and biases
   !       net%layer(i)%nnwts = reshape(vector(net%iwts(i) + 1:net%iwts(i + 1)), &
   !                                    shape=(net%layer(i)%sw), order=[2, 1])
   !       net%layer(i)%nnbias = reshape(vector(net%ibias(i) + 1:net%ibias(i + 1)), &
   !                                     shape=(net%layer(i)%sb), order=[2, 1])
   !       if (net%norm) then
   !          ! get the gamma and beta vectors
   !          net%layer(i)%gmma = vector(net%igmma(i) + 1:net%igmma(i + 1))
   !          net%layer(i)%beta = vector(net%ibeta(i) + 1:net%ibeta(i + 1))
   !          ! net%layer(i)%gmma = 1.d0
   !          ! net%layer(i)%beta = 0.d0
   !       end if

   !    end do

   ! end subroutine vector_to_matrices

   subroutine vector_to_matrices(vector, net)

      ! reshape trial vector (DE output) to neural network weight matrices
      ! this subroutine updates the inupt type(network) object

      type(network), intent(inout) :: net
      real(8), intent(in) :: vector(:)
      integer :: L

      L = net%ld(1)
      net%awts(1, :) = vector(1:L)
      net%omega = vector(L + 1:2*L)

   end subroutine vector_to_matrices

   ! subroutine calc_regularization(net, reg)

   !    type(network), intent(inout) :: net
   !    real(8), intent(out) :: reg
   !    integer :: i

   !    reg = 0.d0

   !    do i = 1, net%nl - 1

   !       if (net%ireg .eq. 0) then
   !          reg = reg + 0.d0

   !       else if (net%ireg .eq. 1) then ! L1
   !          reg = reg + sum(abs(net%layer(i)%nnwts))*net%regconst

   !       else if (net%ireg .eq. 2) then ! L2
   !          reg = reg + sum(net%layer(i)%nnwts**2)*net%regconst
   !       end if

   !    end do

   ! end subroutine calc_regularization

   subroutine calc_regularization(net, reg)

      type(network), intent(inout) :: net
      real(8), intent(out) :: reg

      reg = 0.d0

      if (net%ireg .eq. 0) then
         reg = reg + 0.d0
      else if (net%ireg .eq. 1) then ! L1
         reg = reg + sum(abs(net%awts))*net%regconst
      else if (net%ireg .eq. 2) then ! L2
         reg = reg + sum(net%awts**2)*net%regconst
      end if

   end subroutine calc_regularization

   subroutine normalize_input(x, xnorm, net, idx, calc_mom, mu, sig, gmma, beta)

      ! parameters
      real(8), intent(in) :: x(:, :)
      real(8), allocatable, intent(out) :: xnorm(:, :)
      type(network), intent(inout) :: net
      integer, intent(in) :: idx ! layer index
      logical, intent(in) :: calc_mom
      real(8), intent(in) :: gmma(:), beta(:)
      real(8), optional :: mu(:), sig(:)

      ! internal variables
      real(8), allocatable :: mean(:), sumsqs(:), sigma(:)
      integer :: i, j, nd, nf

      nd = size(x, dim=1) ! number of data
      nf = size(x, dim=2) ! number of features
      allocate (xnorm(nd, nf), mean(nf), sumsqs(nf), sigma(nf))

      ! calculate the moments if required and store in network object
      if (calc_mom) then
         mean = 0.d0
         sumsqs = 0.d0
         do j = 1, nf
            do i = 1, nd
               mean(j) = mean(j) + x(i, j)
               sumsqs(j) = sumsqs(j) + x(i, j)*x(i, j)
            end do
            mean(j) = mean(j)/nd
            sumsqs(j) = sumsqs(j)/nd
            sigma(j) = sqrt(max(sumsqs(j) - mean(j)*mean(j), 0.d0))
         end do
         net%layer(idx)%nnmu = mean
         net%layer(idx)%nnsig = sigma
      end if

      ! use optional moments if supplied
      if (present(mu) .and. present(sig)) then
         mean = mu
         sigma = sig
      end if

      ! now normalize the input matrix
      do j = 1, nf
         xnorm(:, j) = gmma(j)*((x(:, j) - mean(j))/(sigma(j) + EPSLON)) + beta(j)
      end do

   end subroutine normalize_input

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

   function sigmoid1d(yval) result(a)

      ! sigmoid activation

      real(8), intent(in) :: yval(:)
      real(8) :: a(size(yval))

      a = 1.d0/(1.d0 + dexp(-yval))

   end function sigmoid1d

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

   function power(yval, w, f) result(a)

      ! power activation function
      real(8), intent(in) :: yval(:)
      real(8), intent(in) :: w ! trainable exponent
      real(8), intent(in) :: f ! scale factor
      real(8) :: a(size(yval))

      a = yval/f
      where (yval .gt. 0.d0) a = f*(a**w)
      where (yval .lt. 0.d0) a = f*(abs(a)**(1/w)*-(1.d0))

   end function power

end module network_mod
