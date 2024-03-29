module output_mod

   use readpar_mod
   use geostat
   use impute_mod, only: imputed, zinit

   implicit none

contains

   subroutine write_files()

      integer :: i, j, k

      do i = 1, nreals
         do j = 1, ndata
            write (lout, "((i4,1x),*(g14.8,1x))") dhids(j), xyz(1, j), &
               xyz(2, j), xyz(3, j), var(j), imputed(j, ngvarg + 2, i), &
               zinit(j, i), (imputed(j, k, i), k=1, ngvarg + 1)
         end do
      end do

      do i = 1, nnet%nl - 1
         do j = 1, nnet%layer(i)%sb(1)
            write (lmom, "(*(g14.8,1x))") nnet%layer(i)%nnmu(j), &
               nnet%layer(i)%nnsig(j)
         end do
      end do

   end subroutine write_files

end module output_mod
