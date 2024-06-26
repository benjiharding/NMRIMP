module output_mod

   use readpar_mod
   use geostat
   use impute_mod, only: imputed, zinit, abratio

   implicit none

contains

   subroutine write_files()

      integer :: i, j, k

      do i = 1, nreals
         do j = 1, ndata
            write (lout, "((i4,1x),*(g14.8,1x))") dhids(j), xyz(1, j), &
               xyz(2, j), xyz(3, j), var(j), imputed(j, ngvarg + 2, i), &
               zinit(j, i), (imputed(j, k, i), k=1, ngvarg + 1), &
               lb(j), ub(j), abratio(j, i), rsc(j, i), seeded(j, i)
         end do
      end do

      if (nnet%norm) then
         do i = 1, nnet%nl - 1
            do j = 1, nnet%layer(i)%sb(1)
               write (lmom, "(*(g14.8,1x))") nnet%layer(i)%nnmu(j), &
                  nnet%layer(i)%nnsig(j)
            end do
         end do
      end if

   end subroutine write_files

end module output_mod
