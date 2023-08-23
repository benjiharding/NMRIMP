program main

   ! Network model of regionalization imputation. Latent
   ! Gaussian factors are imputed using Sequential Gaussian
   ! Rejection (SGR) Impuation. Imputed factors honour the
   ! variograms specified by the pool and reproduce the hard
   ! data values when passed through the optimized network
   ! generated with NMROPT.
   !
   ! Author: Ben Harding
   ! Date: June 2023
   ! Location: Centre for Computational Geostatistics,
   ! University of Alberta, Edmonton, Canada
   ! Contact: bharding@ualberta.ca

   use readpar_mod
   use mtmod
   use impute_mod, only: impute
   use output_mod, only: write_files

   implicit none

   real, parameter :: VERSION = 1.000

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
   close (ldbg)
   close (lmom)

   ! finished
   write (*, *) " "
   write (*, "(A,f5.3,A)") "NMRIMP version ", VERSION, " finished"

end program main
