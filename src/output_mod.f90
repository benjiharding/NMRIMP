module output_mod

   use readpar_mod
   use geostat
   use impute_mod, only: imputed

   implicit none

contains

   subroutine write_files()

      integer :: i, j

      do i = 1, nreals
         do j = 1, ndata
            write (lout, "((i4,1x),*(g14.8,1x))") dhids(j), xyz(1, j), xyz(2, j), &
               xyz(3, j), var(j), imputed(j, ngvarg + 2, i)
         end do
      end do

   end subroutine write_files

end module output_mod
