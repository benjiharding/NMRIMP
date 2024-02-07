module readpar_mod

   use geostat
   use makepar_mod
   use network_mod, only: init_network, vector_to_matrices
   use vario_mod, only: set_sill, set_rotmatrix
   use constants, only: EPSLON, MAXGNST
   use subs, only: nscore

   implicit none

contains

   subroutine readpar()

      character(256), parameter :: parfile = 'nmrimp.par'
      character(256) :: datafile
      character(256) :: poolfile
      character(256) :: dbgfile
      character(256) :: outfile
      character(256) :: nnwtsfile
      character(256) :: str
      integer :: lin
      logical :: testfl
      integer :: test, i, j, iv, tmp
      integer :: dhcol, xyzcols(3), varcol, wtcol, ncols, numwts
      real(8) :: tmin, tmax
      real(8), allocatable :: nnwts(:)
      real(8), allocatable :: tmpvar(:), tmpnsvar(:), tmpwts(:)

      ! unit numbers
      lin = 1
      lout = 2
      ldbg = 3
      lmom = 4

      ! attempt to open the parfile
      do i = 1, 256
         str(i:i) = ' '
      end do
      call getarg(1, str)

      ! if the input argument is empty request a parameter file name
      if (str .eq. "") then
         write (*, "('Which parameter file do you want to use?')")
         read (*, '(a)') str
      end if
      str = trim(adjustl(str))
      if (str .eq. "") str = parfile

      ! verify if the parameter file exists
      inquire (file=str, exist=testfl)

      ! create a blank parameter file if required
      if (.not. testfl) then
         print *, "ERROR - the parameter file does not exist"
         print *
         if (str .eq. parfile) then
            print *, "        creating a blank parameter file"
            call makepar(parfile)
            print *
         end if
         stop
      end if

      ! open the parfile
      open (lin, file=str, status="old")
      read (lin, '(a4)', iostat=test) str(1:4)
      if (test .ne. 0) stop "ERROR in parameter file"

      ! find the start
      do while (str(1:4) .ne. 'STAR')
         read (lin, '(a4)', iostat=test) str(1:4)
         if (test .ne. 0) stop "ERROR in parameter file"
      end do

      ! read input data file
      write (*, *) " "
      write (*, *) " reading parameter file..."
      write (*, *) " "

      read (lin, '(a256)', iostat=test) datafile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(datafile, 256)
      inquire (file=datafile, exist=testfl)
      if (.not. testfl) stop "ERROR - the data file does not exist"
      write (*, "(2a)") '  data file: ', trim(adjustl(datafile))

      ! read columns for coordinates and values
      read (lin, *, iostat=test) dhcol, xyzcols, varcol, wtcol
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,6(i0,x))") '  column for dhid, x, y, z, var and wt: ', dhcol, xyzcols, varcol, wtcol
      if (dhcol .le. 0) then
         write (*, *) "ERROR: Column for dhid must be > 0. Drill hole IDs are required for sequences."
         stop
      end if

      ! nscore flag
      read (lin, *, iostat=test) itrans
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' normal score transform flag: ', itrans

      ! trimming limits
      read (lin, *, iostat=test) tmin, tmax
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' trimming limits: ', tmin, tmax

      ! number of  reals
      read (lin, *, iostat=test) nreals
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of unconditional realizations: ', nreals

      ! random number seed
      read (lin, *, iostat=test) rseed
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' random number seed: ', rseed

      ! debugging level
      read (lin, *, iostat=test) idbg
      if (test .ne. 0) stop "ERROR in parameter file"

      ! debugging output
      read (lin, '(a256)', iostat=test) dbgfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(dbgfile, 256)
      write (*, "(2a)") '  debugging file: ', trim(adjustl(dbgfile))
      if (idbg .gt. 0) then
         open (ldbg, file=dbgfile, status="UNKNOWN")
         write (ldbg, "(A)") "Reference Distribution"
         write (ldbg, "(i1)") 3
         write (ldbg, "(A)") "CDF y"
         write (ldbg, "(A)") "Raw Activation"
         write (ldbg, "(A)") "NScore Activation"
      end if

      ! imputation output
      read (lin, '(a256)', iostat=test) outfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(outfile, 256)
      write (*, "(2a)") '  output file: ', trim(adjustl(outfile))

      ! network layers
      read (lin, *, iostat=test) nnet%nl
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of network layers: ', nnet%nl

      allocate (nnet%ld(nnet%nl), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! network dimensions
      read (lin, *, iostat=test) nnet%ld
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  network layer dimensions: ', nnet%ld

      ! activation function
      read (lin, *, iostat=test) nnet%af
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  activation function: ', nnet%af

      ! batch normalization
      nnet%norm = .false.
      read (lin, *, iostat=test) tmp
      if (test .ne. 0) stop "ERROR in parameter file"
      if (tmp .gt. 0) nnet%norm = .true.
      write (*, *) '  normalize layer inputs?: ', nnet%norm

      ! no regularization at test time
      nnet%ireg = 0

      ! initilaze the network
      call init_network(nnet)
      allocate (nnwts(nnet%dims), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! network wts input
      read (lin, '(a256)', iostat=test) nnwtsfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(nnwtsfile, 256)
      write (*, "(2a)") '  output file: ', trim(adjustl(nnwtsfile))

      ! open the output file for layer moments if required
      if (nnet%norm) then
         open (lmom, file="moments.out", status="UNKNOWN")
         write (lmom, "(A)") "Network Moments"
         write (lmom, "(i2)") 2
         write (lmom, "(A)") "mu"
         write (lmom, "(A)") "sigma"
      end if

      ! previously simulated nodes to consider
      read (lin, *, iostat=test) nsearch
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of previously simulated nodes: ', nsearch

      ! max iterations for imputation steps
      read (lin, *, iostat=test) iter1, iter2
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' max iterations for step 1 and 2: ', iter1, iter2

      ! tolerances for imputations steps
      read (lin, *, iostat=test) tol1, tol2
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' rejection tolerances for step 1 and 2: ', tol1, tol2

      ! Gaussian pool file
      read (lin, '(a256)', iostat=test) poolfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(poolfile, 256)
      inquire (file=poolfile, exist=testfl)
      if (.not. testfl) stop "ERROR - the data file does not exist"
      write (*, "(2a)") '  file with Gaussian pool cov. struct.: ', &
         trim(adjustl(poolfile))

      ! finished reading parameters
      close (lin)

      !
      ! start reading the data file
      !
      write (*, *) " "
      write (*, *) " reading data file..."
      write (*, *) " "

      open (lin, file=datafile, status='OLD')
      read (lin, *)
      read (lin, *, iostat=test) ncols
      if (test .ne. 0) stop "ERROR in data file"

      ! check column numbers
      if (dhcol .gt. ncols .or. &
          any(xyzcols .gt. ncols) .or. &
          (varcol .gt. ncols) .or. &
          (wtcol .gt. ncols)) then
         write (*, *) 'there are only ', ncols, ' columns in the data file'
         write (*, *) '  your specification is out of range'
         stop
      end if
      allocate (tmpvar(ncols), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! jump headers and get names
      do i = 1, ncols
         read (lin, *, iostat=test) str
         if (test .ne. 0) stop "ERROR in data file"
      end do

      ! get the number of data
      ndata = 0
      do
         read (lin, *, iostat=test) tmpvar(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         ndata = ndata + 1
      end do

      ! allocate arrays for input data
      allocate (dhids(ndata), xyz(3, ndata), var(ndata), wts(ndata), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! restart data file
      rewind (lin)
      do i = 1, ncols + 2
         read (lin, *, iostat=test)
      end do

      ! read file again, but storing variable and weight
      ndata = 0
      do
         read (lin, *, iostat=test) tmpvar(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         ndata = ndata + 1
         dhids(ndata) = tmpvar(dhcol)
         xyz(:, ndata) = tmpvar(xyzcols)
         var(ndata) = tmpvar(varcol)
         wts(ndata) = tmpvar(wtcol)
      end do

      ! assume equal weighting if not specified
      if (wtcol .eq. 0) then
         wts = 1.d0
      end if

      ! nscore input var if required
      if (itrans .eq. 1) then
         call nscore(ndata, var, tmin, tmax, 1, wts, tmpnsvar, nsvar, test)
         if (test .gt. 0) stop "ERROR in normal score transform"
         var = nsvar
      end if

      write (*, *) " "
      write (*, *) " reading optimized network weights..."
      write (*, *) " "

      open (lin, file=nnwtsfile, status='OLD')
      read (lin, *)
      read (lin, *, iostat=test) ncols
      if (test .ne. 0) stop "ERROR in data file"

      ! jump headers and get names
      do i = 1, ncols
         read (lin, *, iostat=test) str
         if (test .ne. 0) stop "ERROR in data file"
      end do

      ! get the number of weights
      numwts = 0
      do
         read (lin, *, iostat=test) tmpwts(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         numwts = numwts + 1
      end do

      ! check weights match specified network architecture
      if (numwts .ne. nnet%dims) stop "ERROR Network weight dimensions do not &
&      match the specified network architecture!"

      allocate (tmpwts(1), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! restart weights file
      rewind (lin)
      do i = 1, ncols + 2
         read (lin, *, iostat=test)
      end do

      ! read file again, but storing weights
      numwts = 0
      do
         read (lin, *, iostat=test) tmpwts(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         numwts = numwts + 1
         nnwts(numwts) = tmpwts(1)
      end do

      ! reshape the weight vector into matrices
      call vector_to_matrices(nnwts, nnet)

      ! start reading Gaussian pool file
      write (*, *) " "
      write (*, *) " reading covariance structure of Gaussian pool..."
      write (*, *) " "

      ! allocate arrays for the pool
      ngvarg = nnet%ld(1) - 1
      allocate (pool(ngvarg), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! open the output file and write headers now that we know ngvarg
      open (lout, file=outfile, status="UNKNOWN")
      write (lout, "(A)") "Imputed Data Values"
      write (lout, "(i2)") 7 + ngvarg + 1
      write (lout, "(A)") "dhid"
      write (lout, "(A)") "x"
      write (lout, "(A)") "y"
      write (lout, "(A)") "z"
      write (lout, "(A)") "Data value"
      write (lout, "(A)") "Imputed value"
      write (lout, "(A)") "Raw activation"
      do iv = 1, ngvarg + 1
         write (lout, "(a6, i3)") "Factor", iv
      end do

      ! open the pool file
      open (lin, file=poolfile, status='OLD')

      ! parse the Gaussian variogram models
      do iv = 1, ngvarg

         read (lin, *, iostat=test) pool(iv)%nst, pool(iv)%c0
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  gnst, gc0: ', pool(iv)%nst, pool(iv)%c0

         if (pool(iv)%nst .gt. MAXGNST) then
            write (*, *) 'gnst must be equal to ', MAXGNST
            stop
         end if

         allocate (pool(iv)%it(pool(iv)%nst), pool(iv)%cc(pool(iv)%nst), &
                   pool(iv)%ang1(pool(iv)%nst), pool(iv)%ang2(pool(iv)%nst), pool(iv)%ang3(pool(iv)%nst), &
                   pool(iv)%aa(pool(iv)%nst), pool(iv)%anis1(pool(iv)%nst), pool(iv)%anis2(pool(iv)%nst), &
                   pool(iv)%ahmin(pool(iv)%nst), pool(iv)%avert(pool(iv)%nst), stat=test)
         if (test .ne. 0) stop "allocation failed due to insufficient memory!"

         do j = 1, pool(iv)%nst
            read (lin, *, iostat=test) pool(iv)%it(j), pool(iv)%cc(j), pool(iv)%ang1(j), &
               pool(iv)%ang2(j), pool(iv)%ang3(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            read (lin, *, iostat=test) pool(iv)%aa(j), pool(iv)%ahmin(j), pool(iv)%avert(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            pool(iv)%anis1(j) = pool(iv)%ahmin(j)/max(pool(iv)%aa(j), EPSLON)
            pool(iv)%anis2(j) = pool(iv)%avert(j)/max(pool(iv)%aa(j), EPSLON)
            write (*, *) ' git, gcc, gang[1,2,3]; ', pool(iv)%it(j), pool(iv)%cc(j), &
               pool(iv)%ang1(j), pool(iv)%ang2(j), pool(iv)%ang3(j)
            write (*, *) ' a1 a2 a3: ', pool(iv)%aa(j), pool(iv)%ahmin(j), pool(iv)%avert(j)
         end do
      end do

      call set_sill(pool)
      call set_rotmatrix(pool)

      ! finished reading data
      close (lin)

   end subroutine readpar

end module readpar_mod
