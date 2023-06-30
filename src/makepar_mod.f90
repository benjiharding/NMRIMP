module makepar_mod

   implicit none

contains

   subroutine makepar(parfile)
      character(256), intent(in) :: parfile

      open (unit=99, file=parfile, action="write", status="replace")

      write (99, "(A)") '               PARAMETERS FOR NMRIMP'
      write (99, "(A)") '               *********************'
      write (99, "(A)") ''
      write (99, "(A)") 'START OF PARAMETERS:'
      write (99, "(A)") 'data.dat                          - file with data'
      write (99, "(A)") '1 4 5 6 7 11                      - columns for dh, x, y, z, var and wt'
      write (99, "(A)") '0                                 - normal score transform var? (0=no, 1=yes)'
      write (99, "(A)") '-1.0e21    1.0e21                 - trimming limits'
      write (99, "(A)") '100                               - number of realizations'
      write (99, "(A)") '5841044                           - random number seed'
      write (99, "(A)") '1                                 - debugging level'
      write (99, "(A)") 'nmrimp.dbg                        - file for debugging output'
      write (99, "(A)") 'nmrimp.out                        - output file with imputed realizations'
      write (99, "(A)") '5                                 - number of network layers (input to output layer)'
      write (99, "(A)") '16 5 5 5 1                        - network layer dimensions (input + nugget to output layer)'
      write (99, "(A)") '1                                 - network activation func. (1=sigmoid, 2=tanh, 3=relu, 4=linear)'
      write (99, "(A)") 'nmrwts.out                        - input file with optimzed network weights'
      write (99, "(A)") '40                                - maximum previously simulated nodes'
      write (99, "(A)") '50000 10000                       - maximum iterations for step 1 and step 2'
      write (99, "(A)") '0.1 0.01                          - rejection tolerances for step 1 and step 2'
      write (99, "(A)") 'pool.dat                          - file with covariance structs. of Gaussian pool'
      close (99)

   end subroutine makepar

   subroutine chknam(str, len)
      !-----------------------------------------------------------------------
      !
      !                   Check for a Valid File Name
      !                   ***************************
      !
      ! This subroutine takes the character string "str" of length "len" and
      ! removes all leading blanks and blanks out all characters after the
      ! first blank found in the string (leading blanks are removed first).
      !
      !
      !
      !-----------------------------------------------------------------------
      character(len=*), intent(inout) :: str
      integer itrim, len
      !
      ! Remove leading blanks:
      str = adjustl(str)
      !
      ! find first two blanks and blank out remaining characters:
      itrim = index(str, '   ')
      if (itrim > 0) str(itrim:) = ' '
      !
      ! Look for "-fi"
      itrim = index(str, '-fi')
      if (itrim > 0) str(itrim:) = ' '
      !
      ! Look for "\fi"
      itrim = index(str, '\fi')
      if (itrim > 0) str(itrim:) = ' '

      return
   end subroutine chknam

end module makepar_mod
