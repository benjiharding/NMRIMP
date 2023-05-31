program main

   use readpar_mod
   use mtmod
   use impute_mod, only: impute
   use output_mod, only: write_files

   implicit none

   ! read the parfile
   call readpar

   ! initilize random generator
   call sgrnd(rseed)

   ! impute the factors
   call impute

   ! write output files
   call write_files

   ! close all output files
   close (lout)

end program main
