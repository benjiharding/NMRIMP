module subs

   use mtmod
   use constants

   implicit none

contains

   subroutine nscore(nd, vr, tmin, tmax, iwt, wt, tmp, vrg, ierror)
      !-----------------------------------------------------------------------
      !
      !              Transform Univariate Data to Normal Scores
      !              ******************************************
      !
      ! This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
      ! by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
      ! as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
      ! so that the data can be returned in the same order (just in case
      ! there are associated arrays like the coordinate location).
      !
      !
      !
      ! INPUT VARIABLES:
      !
      !   nd               Number of data (no missing values)
      !   vr(nd)           Data values to be transformed
      !   tmin,tmax        data trimming limits
      !   iwt              =0, equal weighted; =1, then apply weight
      !   wt(nd)           Weight for each data (don't have to sum to 1.0)
      !   tmp(nd)          Temporary storage space for sorting
      !   lout             if > 0 then transformation table will be written
      !
      !
      !
      ! OUTPUT VARIABLES:
      !
      !   vrg(nd)          normal scores
      !   ierror           error flag (0=error free,1=problem)
      !
      !
      !
      ! EXTERNAL REFERENCES:
      !
      !   gauinv           Calculates the inverse of a Gaussian cdf
      !   sortem           sorts a number of arrays according to a key array
      !
      !
      !
      !-----------------------------------------------------------------------
      real(8), intent(in) :: vr(nd), wt(nd)
      real(8), allocatable :: vrg(:), tmp(:)
      real(8) :: twt, cumwt(nd), oldcp, cp
      real(8) :: tmin, tmax
      integer :: i, iwt, ierr, ierror, nd

      allocate (vrg(nd), tmp(nd))
      cumwt = wt

      !
      ! Sort the data in ascending order and calculate total weight:
      !
      ierror = 0
      twt = 0.d0
      do i = 1, nd
         tmp(i) = real(i)
         if (vr(i) .ge. tmin .and. vr(i) .lt. tmax) then
            if (iwt .eq. 0) then
               twt = twt + 1.
            else
               twt = twt + wt(i)
            end if
         end if
      end do
      if (nd .lt. 1 .or. twt .lt. EPSLON) then
         ierror = 1
         return
      end if
      call sortem(1, nd, vr, 2, cumwt, tmp, vr, vr, vr, vr, vr, vr)
      !
      ! Compute the cumulative probabilities:
      !
      oldcp = 0.d0
      cp = 0.d0
      do i = 1, nd
         cp = cp + cumwt(i)/twt
         cumwt(i) = (cp + oldcp)/2.0
         oldcp = cp
         call gauinv(dble(cumwt(i)), vrg(i), ierr)
      end do
      !
      ! Get the arrays back in original order:
      !
      call sortem(1, nd, tmp, 3, cumwt, vr, vrg, vr, vr, vr, vr, vr)
      !
      ! Finished:
      !
      return
   end

   subroutine locate(xx, n, is, ie, x, j)

      ! Given an array "xx" of length "n", and given a value "x", this routine
      ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
      ! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
      ! returned to indicate that x is out of range.
      !
      ! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.

      integer, intent(in) :: is, ie, n
      real(8), intent(in) :: xx(n), x

      ! local variables
      integer, intent(out) :: j
      integer :: jl, ju, jm

      ! Initialize lower and upper methods:

      !   if (is .le. 0) then
      !      is = 1
      !   end if

      jl = is - 1
      ju = ie
      if (xx(n) .le. x) then
         j = ie
         return
      end if

      ! If we are not done then compute a midpoint:

10    if (ju - jl .gt. 1) then
         jm = (ju + jl)/2

         ! Replace the lower or upper limit with the midpoint:
         if ((xx(ie) .gt. xx(is)) .eqv. (x .gt. xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
         go to 10
      end if

      ! Return with the array index:
      j = jl
   end

   subroutine sortem(ib, ie, a, iperm, b, c, d, e, f, g, h, aa)
      !-----------------------------------------------------------------------
      !
      !                      Quickersort Subroutine
      !                      **********************
      !
      ! This is a subroutine for sorting a real array in ascending order. This
      ! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
      ! in collected algorithms of the ACM.
      !
      ! The method used is that of continually splitting the array into parts
      ! such that all elements of one part are less than all elements of the
      ! other, with a third part in the middle consisting of one element.  An
      ! element with value t is chosen arbitrarily (here we choose the middle
      ! element). i and j give the lower and upper limits of the segment being
      ! split.  After the split a value q will have been found such that
      ! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
      ! performs operations on the two segments (i,q-1) and (q+1,j) as follows
      ! The smaller segment is split and the position of the larger segment is
      ! stored in the lt and ut arrays.  If the segment to be split contains
      ! two or fewer elements, it is sorted and another segment is obtained
      ! from the lt and ut arrays.  When no more segments remain, the array
      ! is completely sorted.
      !
      !
      ! INPUT PARAMETERS:
      !
      !   ib,ie        start and end index of the array to be sorteda
      !   a            array, a portion of which has to be sorted.
      !   iperm        0 no other array is permuted.
      !                1 array b is permuted according to array a
      !                2 arrays b,c are permuted.
      !                3 arrays b,c,d are permuted.
      !                4 arrays b,c,d,e are permuted.
      !                5 arrays b,c,d,e,f are permuted.
      !                6 arrays b,c,d,e,f,g are permuted.
      !                7 arrays b,c,d,e,f,g,h are permuted.
      !                8 arrays b,c,d,e,f,g,h,aa are permuted.
      !               >8 no other array is permuted.
      !
      !   b,c,d,e,f,g,h,aa  arrays to be permuted according to array a.
      !
      ! OUTPUT PARAMETERS:
      !
      !    a      = the array, a portion of which has been sorted.
      !
      !    b,c,d,e,f,g,h,aa  =arrays permuted according to array a (see iperm)
      !
      ! NO EXTERNAL ROUTINES REQUIRED:
      !
      !-----------------------------------------------------------------------
      implicit none
      real*8 a(*), b(*), c(*), d(*), e(*), f(*), g(*), h(*), aa(*)
      real*8 ta, tb, tc, td, te, tf, tg, th, taa, xa, xb, xc, xd, xe, xf, xg, xh, xaa

      ! The dimensions for lt and ut have to be at least log (base 2) n

      integer lt(64), ut(64), i, j, k, m, p, q, ie, ib, iring, iperm
      !
      ! Initialize:
      !
      j = ie
      m = 1
      i = ib
      iring = iperm + 1
      if (iperm .gt. 7) iring = 1
      !
      ! If this segment has more than two elements  we split it
      !
10    if (j - i - 1) 100, 90, 15
      !
      ! p is the position of an arbitrary element in the segment we choose the
      ! middle element. Under certain circumstances it may be advantageous
      ! to choose p at random.
      !
15    p = (j + i)/2
      ta = a(p)
      a(p) = a(i)
      go to(21, 19, 18, 17, 16, 161, 162, 163, 164), iring
164   taa = aa(p)
      aa(p) = aa(i)
163   th = h(p)
      h(p) = h(i)
162   tg = g(p)
      g(p) = g(i)
161   tf = f(p)
      f(p) = f(i)
16    te = e(p)
      e(p) = e(i)
17    td = d(p)
      d(p) = d(i)
18    tc = c(p)
      c(p) = c(i)
19    tb = b(p)
      b(p) = b(i)
21    continue
      !
      ! Start at the beginning of the segment, search for k such that a(k)>t
      !
      q = j
      k = i
20    k = k + 1
      if (k .gt. q) go to 60
      if (a(k) .le. ta) go to 20
      !
      ! Such an element has now been found now search for a q such that a(q)<t
      ! starting at the end of the segment.
      !
30    continue
      if (a(q) .lt. ta) go to 40
      q = q - 1
      if (q .gt. k) go to 30
      go to 50
      !
      ! a(q) has now been found. we interchange a(q) and a(k)
      !
40    xa = a(k)
      a(k) = a(q)
      a(q) = xa
      go to(45, 44, 43, 42, 41, 411, 412, 413, 414), iring
414   xaa = aa(k)
      aa(k) = aa(q)
      aa(q) = xaa
413   xh = h(k)
      h(k) = h(q)
      h(q) = xh
412   xg = g(k)
      g(k) = g(q)
      g(q) = xg
411   xf = f(k)
      f(k) = f(q)
      f(q) = xf
41    xe = e(k)
      e(k) = e(q)
      e(q) = xe
42    xd = d(k)
      d(k) = d(q)
      d(q) = xd
43    xc = c(k)
      c(k) = c(q)
      c(q) = xc
44    xb = b(k)
      b(k) = b(q)
      b(q) = xb
45    continue
      !
      ! Update q and search for another pair to interchange:
      !
      q = q - 1
      go to 20
50    q = k - 1
60    continue
      !
      ! The upwards search has now met the downwards search:
      !
      a(i) = a(q)
      a(q) = ta
      go to(65, 64, 63, 62, 61, 611, 612, 613, 614), iring
614   aa(i) = aa(q)
      aa(q) = taa
613   h(i) = h(q)
      h(q) = th
612   g(i) = g(q)
      g(q) = tg
611   f(i) = f(q)
      f(q) = tf
61    e(i) = e(q)
      e(q) = te
62    d(i) = d(q)
      d(q) = td
63    c(i) = c(q)
      c(q) = tc
64    b(i) = b(q)
      b(q) = tb
65    continue
      !
      ! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
      ! store the position of the largest segment in lt and ut
      !
      if (2*q .le. i + j) go to 70
      lt(m) = i
      ut(m) = q - 1
      i = q + 1
      go to 80
70    lt(m) = q + 1
      ut(m) = j
      j = q - 1
      !
      ! Update m and split the new smaller segment
      !
80    m = m + 1
      go to 10
      !
      ! We arrive here if the segment has  two elements we test to see if
      ! the segment is properly ordered if not, we perform an interchange
      !
90    continue
      if (a(i) .le. a(j)) go to 100
      xa = a(i)
      a(i) = a(j)
      a(j) = xa
      go to(95, 94, 93, 92, 91, 911, 912, 913, 914), iring
914   xaa = aa(i)
      aa(i) = aa(j)
      aa(j) = xaa
913   xh = h(i)
      h(i) = h(j)
      h(j) = xh
912   xg = g(i)
      g(i) = g(j)
      g(j) = xg
911   xf = f(i)
      f(i) = f(j)
      f(j) = xf
91    xe = e(i)
      e(i) = e(j)
      e(j) = xe
92    xd = d(i)
      d(i) = d(j)
      d(j) = xd
93    xc = c(i)
      c(i) = c(j)
      c(j) = xc
94    xb = b(i)
      b(i) = b(j)
      b(j) = xb
95    continue
      !
      ! If lt and ut contain more segments to be sorted repeat process:
      !
100   m = m - 1
      if (m .le. 0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
110   continue

   end subroutine sortem

   subroutine modsortem(a, n, iperm, b, c, d, e, f, g, h)
      !-----------------------------------------------------------------------
      !
      ! WARNING: this is a modified version of sortem - the array length
      !          is specified - this is compatible with slicing
      !
      !                      Quickersort Subroutine
      !                      **********************
      !
      ! This is a subroutine for sorting a real array in ascending order. This
      ! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
      ! in collected algorithms of the ACM.
      !
      ! The method used is that of continually splitting the array into parts
      ! such that all elements of one part are less than all elements of the
      ! other, with a third part in the middle consisting of one element.  An
      ! element with value t is chosen arbitrarily (here we choose the middle
      ! element). i and j give the lower and upper limits of the segment being
      ! split.  After the split a value q will have been found such that
      ! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
      ! performs operations on the two segments (i,q-1) and (q+1,j) as follows
      ! The smaller segment is split and the position of the larger segment is
      ! stored in the lt and ut arrays.  If the segment to be split contains
      ! two or fewer elements, it is sorted and another segment is obtained
      ! from the lt and ut arrays.  When no more segments remain, the array
      ! is completely sorted.
      !
      !
      ! INPUT PARAMETERS:
      !
      !   ib,ie        start and end index of the array to be sorteda
      !   a            array, a portion of which has to be sorted.
      !   iperm        0 no other array is permuted.
      !                1 array b is permuted according to array a
      !                2 arrays b,c are permuted.
      !                3 arrays b,c,d are permuted.
      !                4 arrays b,c,d,e are permuted.
      !                5 arrays b,c,d,e,f are permuted.
      !                6 arrays b,c,d,e,f,g are permuted.
      !                7 arrays b,c,d,e,f,g,h are permuted.
      !               >7 no other array is permuted.
      !
      !   b,c,d,e,f,g,h  arrays to be permuted according to array a.
      !
      ! OUTPUT PARAMETERS:
      !
      !    a      = the array, a portion of which has been sorted.
      !
      !    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
      !
      ! NO EXTERNAL ROUTINES REQUIRED:
      !
      !-----------------------------------------------------------------------
      !dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
      integer, intent(in) :: n, iperm
      real, dimension(n), intent(inout) :: a
      real, dimension(n), optional, intent(inout) :: b, c, d, e, f, g, h
      !
      ! The dimensions for lt and ut have to be at least log (base 2) n
      !
      integer :: lt(64), ut(64), i, j, k, m, p, q, ib, ie, iring
      real :: ta, tb, tc, td, te, tf, tg, th
      real :: xa, xb, xc, xd, xe, xf, xg, xh

      ib = 1
      ie = n
      !
      ! Initialize:
      !
      j = ie
      m = 1
      i = ib
      iring = iperm + 1
      if (iperm .gt. 7) iring = 1
      !
      ! If this segment has more than two elements  we split it
      !
10    if (j - i - 1) 100, 90, 15
      !
      ! p is the position of an arbitrary element in the segment we choose the
      ! middle element. Under certain circumstances it may be advantageous
      ! to choose p at random.
      !
15    p = (j + i)/2
      ta = a(p)
      a(p) = a(i)
      go to(21, 19, 18, 17, 16, 161, 162, 163), iring
163   th = h(p)
      h(p) = h(i)
162   tg = g(p)
      g(p) = g(i)
161   tf = f(p)
      f(p) = f(i)
16    te = e(p)
      e(p) = e(i)
17    td = d(p)
      d(p) = d(i)
18    tc = c(p)
      c(p) = c(i)
19    tb = b(p)
      b(p) = b(i)
21    continue
      !
      ! Start at the beginning of the segment, search for k such that a(k)>t
      !
      q = j
      k = i
20    k = k + 1
      if (k .gt. q) go to 60
      if (a(k) .le. ta) go to 20
      !
      ! Such an element has now been found now search for a q such that a(q)<t
      ! starting at the end of the segment.
      !
30    continue
      if (a(q) .lt. ta) go to 40
      q = q - 1
      if (q .gt. k) go to 30
      go to 50
      !
      ! a(q) has now been found. we interchange a(q) and a(k)
      !
40    xa = a(k)
      a(k) = a(q)
      a(q) = xa
      go to(45, 44, 43, 42, 41, 411, 412, 413), iring
413   xh = h(k)
      h(k) = h(q)
      h(q) = xh
412   xg = g(k)
      g(k) = g(q)
      g(q) = xg
411   xf = f(k)
      f(k) = f(q)
      f(q) = xf
41    xe = e(k)
      e(k) = e(q)
      e(q) = xe
42    xd = d(k)
      d(k) = d(q)
      d(q) = xd
43    xc = c(k)
      c(k) = c(q)
      c(q) = xc
44    xb = b(k)
      b(k) = b(q)
      b(q) = xb
45    continue
      !
      ! Update q and search for another pair to interchange:
      !
      q = q - 1
      go to 20
50    q = k - 1
60    continue
      !
      ! The upwards search has now met the downwards search:
      !
      a(i) = a(q)
      a(q) = ta
      go to(65, 64, 63, 62, 61, 611, 612, 613), iring
613   h(i) = h(q)
      h(q) = th
612   g(i) = g(q)
      g(q) = tg
611   f(i) = f(q)
      f(q) = tf
61    e(i) = e(q)
      e(q) = te
62    d(i) = d(q)
      d(q) = td
63    c(i) = c(q)
      c(q) = tc
64    b(i) = b(q)
      b(q) = tb
65    continue
      !
      ! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
      ! store the position of the largest segment in lt and ut
      !
      if (2*q .le. i + j) go to 70
      lt(m) = i
      ut(m) = q - 1
      i = q + 1
      go to 80
70    lt(m) = q + 1
      ut(m) = j
      j = q - 1
      !
      ! Update m and split the new smaller segment
      !
80    m = m + 1
      go to 10
      !
      ! We arrive here if the segment has  two elements we test to see if
      ! the segment is properly ordered if not, we perform an interchange
      !
90    continue
      if (a(i) .le. a(j)) go to 100
      xa = a(i)
      a(i) = a(j)
      a(j) = xa
      go to(95, 94, 93, 92, 91, 911, 912, 913), iring
913   xh = h(i)
      h(i) = h(j)
      h(j) = xh
912   xg = g(i)
      g(i) = g(j)
      g(j) = xg
911   xf = f(i)
      f(i) = f(j)
      f(j) = xf
91    xe = e(i)
      e(i) = e(j)
      e(j) = xe
92    xd = d(i)
      d(i) = d(j)
      d(j) = xd
93    xc = c(i)
      c(i) = c(j)
      c(j) = xc
94    xb = b(i)
      b(i) = b(j)
      b(j) = xb
95    continue
      !
      ! If lt and ut contain more segments to be sorted repeat process:
      !
100   m = m - 1
      if (m .le. 0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
110   continue
      return
   end subroutine modsortem

   subroutine dblemodsortem(a, n, iperm, b, c, d, e, f, g, h)
      !-----------------------------------------------------------------------
      !
      ! WARNING: this is a modified version of sortem - the array length
      !          is specified - this is compatible with slicing
      !
      !                      Quickersort Subroutine
      !                      **********************
      !
      ! This is a subroutine for sorting a real array in ascending order. This
      ! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
      ! in collected algorithms of the ACM.
      !
      ! The method used is that of continually splitting the array into parts
      ! such that all elements of one part are less than all elements of the
      ! other, with a third part in the middle consisting of one element.  An
      ! element with value t is chosen arbitrarily (here we choose the middle
      ! element). i and j give the lower and upper limits of the segment being
      ! split.  After the split a value q will have been found such that
      ! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
      ! performs operations on the two segments (i,q-1) and (q+1,j) as follows
      ! The smaller segment is split and the position of the larger segment is
      ! stored in the lt and ut arrays.  If the segment to be split contains
      ! two or fewer elements, it is sorted and another segment is obtained
      ! from the lt and ut arrays.  When no more segments remain, the array
      ! is completely sorted.
      !
      !
      ! INPUT PARAMETERS:
      !
      !   ib,ie        start and end index of the array to be sorteda
      !   a            array, a portion of which has to be sorted.
      !   iperm        0 no other array is permuted.
      !                1 array b is permuted according to array a
      !                2 arrays b,c are permuted.
      !                3 arrays b,c,d are permuted.
      !                4 arrays b,c,d,e are permuted.
      !                5 arrays b,c,d,e,f are permuted.
      !                6 arrays b,c,d,e,f,g are permuted.
      !                7 arrays b,c,d,e,f,g,h are permuted.
      !               >7 no other array is permuted.
      !
      !   b,c,d,e,f,g,h  arrays to be permuted according to array a.
      !
      ! OUTPUT PARAMETERS:
      !
      !    a      = the array, a portion of which has been sorted.
      !
      !    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
      !
      ! NO EXTERNAL ROUTINES REQUIRED:
      !
      !-----------------------------------------------------------------------
      !dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
      integer, intent(in) :: n, iperm
      real(kind=8), dimension(n), intent(inout) :: a
      real(kind=8), dimension(n), optional, intent(inout) :: b, c, d, e, f, g, h
      real(kind=8) th, tg, tf, te, td, tc, tb, ta, xh, xg, xf, xe, xd, xc, xb, xa
      !
      ! The dimensions for lt and ut have to be at least log (base 2) n
      !
      integer :: lt(64), ut(64), i, j, k, m, p, q, ib, ie, iring

      ib = 1
      ie = n
      !
      ! Initialize:
      !
      j = ie
      m = 1
      i = ib
      iring = iperm + 1
      if (iperm .gt. 7) iring = 1
      !
      ! If this segment has more than two elements  we split it
      !
10    if (j - i - 1) 100, 90, 15
      !
      ! p is the position of an arbitrary element in the segment we choose the
      ! middle element. Under certain circumstances it may be advantageous
      ! to choose p at random.
      !
15    p = (j + i)/2
      ta = a(p)
      a(p) = a(i)
      go to(21, 19, 18, 17, 16, 161, 162, 163), iring
163   th = h(p)
      h(p) = h(i)
162   tg = g(p)
      g(p) = g(i)
161   tf = f(p)
      f(p) = f(i)
16    te = e(p)
      e(p) = e(i)
17    td = d(p)
      d(p) = d(i)
18    tc = c(p)
      c(p) = c(i)
19    tb = b(p)
      b(p) = b(i)
21    continue
      !
      ! Start at the beginning of the segment, search for k such that a(k)>t
      !
      q = j
      k = i
20    k = k + 1
      if (k .gt. q) go to 60
      if (a(k) .le. ta) go to 20
      !
      ! Such an element has now been found now search for a q such that a(q)<t
      ! starting at the end of the segment.
      !
30    continue
      if (a(q) .lt. ta) go to 40
      q = q - 1
      if (q .gt. k) go to 30
      go to 50
      !
      ! a(q) has now been found. we interchange a(q) and a(k)
      !
40    xa = a(k)
      a(k) = a(q)
      a(q) = xa
      go to(45, 44, 43, 42, 41, 411, 412, 413), iring
413   xh = h(k)
      h(k) = h(q)
      h(q) = xh
412   xg = g(k)
      g(k) = g(q)
      g(q) = xg
411   xf = f(k)
      f(k) = f(q)
      f(q) = xf
41    xe = e(k)
      e(k) = e(q)
      e(q) = xe
42    xd = d(k)
      d(k) = d(q)
      d(q) = xd
43    xc = c(k)
      c(k) = c(q)
      c(q) = xc
44    xb = b(k)
      b(k) = b(q)
      b(q) = xb
45    continue
      !
      ! Update q and search for another pair to interchange:
      !
      q = q - 1
      go to 20
50    q = k - 1
60    continue
      !
      ! The upwards search has now met the downwards search:
      !
      a(i) = a(q)
      a(q) = ta
      go to(65, 64, 63, 62, 61, 611, 612, 613), iring
613   h(i) = h(q)
      h(q) = th
612   g(i) = g(q)
      g(q) = tg
611   f(i) = f(q)
      f(q) = tf
61    e(i) = e(q)
      e(q) = te
62    d(i) = d(q)
      d(q) = td
63    c(i) = c(q)
      c(q) = tc
64    b(i) = b(q)
      b(q) = tb
65    continue
      !
      ! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
      ! store the position of the largest segment in lt and ut
      !
      if (2*q .le. i + j) go to 70
      lt(m) = i
      ut(m) = q - 1
      i = q + 1
      go to 80
70    lt(m) = q + 1
      ut(m) = j
      j = q - 1
      !
      ! Update m and split the new smaller segment
      !
80    m = m + 1
      go to 10
      !
      ! We arrive here if the segment has  two elements we test to see if
      ! the segment is properly ordered if not, we perform an interchange
      !
90    continue
      if (a(i) .le. a(j)) go to 100
      xa = a(i)
      a(i) = a(j)
      a(j) = xa
      go to(95, 94, 93, 92, 91, 911, 912, 913), iring
913   xh = h(i)
      h(i) = h(j)
      h(j) = xh
912   xg = g(i)
      g(i) = g(j)
      g(j) = xg
911   xf = f(i)
      f(i) = f(j)
      f(j) = xf
91    xe = e(i)
      e(i) = e(j)
      e(j) = xe
92    xd = d(i)
      d(i) = d(j)
      d(j) = xd
93    xc = c(i)
      c(i) = c(j)
      c(j) = xc
94    xb = b(i)
      b(i) = b(j)
      b(j) = xb
95    continue
      !
      ! If lt and ut contain more segments to be sorted repeat process:
      !
100   m = m - 1
      if (m .le. 0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
110   continue
      return
   end subroutine dblemodsortem

   function powint(xlow, xhigh, ylow, yhigh, xval, pow) result(y)

      ! Power interpolate the value of y between (xlow,ylow) and
      ! (xhigh,yhigh) for a value of x and a power pow

      real(8) :: xlow, xhigh, ylow, yhigh, xval, pow, y

      if ((xhigh - xlow) .lt. EPSLON) then
         y = (yhigh + ylow)/2.D0
      else
         y = ylow + (yhigh - ylow)*(((xval - xlow)/(xhigh - xlow))**pow)
      end if

   end function powint

   subroutine gauinv(p, xp, ierr)
      !-----------------------------------------------------------------------
      !
      ! Computes the inverse of the standard normal cumulative distribution
      ! function with a numerical approximation from : Statistical Computing,
      ! by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
      !
      !
      !
      ! INPUT/OUTPUT:
      !
      !   p    = double precision cumulative probability value: dble(psingle)
      !   xp   = G^-1 (p) in single precision
      !   ierr = 1 - then error situation (p out of range), 0 - OK
      !
      !
      !-----------------------------------------------------------------------
      implicit none
      integer ierr
      real*8 p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, y, pp, lim, p, xp
      save p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, lim
      !
      ! Coefficients of approximation:
      !
      data lim/1.0d-10/
      data p0/-0.322232431088D0/, p1/-1.0D0/, p2/-0.342242088547D0/, &
         p3/-0.0204231210245D0/, p4/-0.0000453642210148D0/
      data q0/0.0993484626060D0/, q1/0.588581570495D0/, q2/ &
         0.531103462366D0/, &
         q3/0.10353775285D0/, q4/0.0038560700634D0/
      !
      ! Check for an error situation:
      !
      ierr = 1
      if (p .lt. lim) then
         xp = -1.0d10
         return
      end if
      if (p .gt. (1.0 - lim)) then
         xp = 1.0d10
         return
      end if
      ierr = 0
      !
      ! Get k for an error situation:
      !
      pp = p
      if (p .gt. 0.5D0) pp = 1 - pp
      xp = 0.D0
      if (p .eq. 0.5D0) return
      !
      ! Approximate the function:
      !
      y = dsqrt(dlog(1.0/(pp*pp)))
      xp = dble(y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0)/ &
                ((((y*q4 + q3)*y + q2)*y + q1)*y + q0))
      if (dble(p) .eq. dble(pp)) xp = -xp

   end subroutine gauinv

   subroutine random_stduniform(u)
      real(8), intent(out) :: u
      real(8) :: r
      call random_number(r)
      u = 1 - r
   end subroutine random_stduniform

   ! assuming a<b https://masuday.github.io/fortran_tutorial/random.html
   subroutine random_uniform(a, b, x)
      real(8), intent(in) :: a, b
      real(8), intent(out) :: x
      real(8) :: u
      call random_stduniform(u)
      x = (b - a)*u + a
   end subroutine random_uniform

   !https://stackoverflow.com/questions/44198212/a-fortran-equivalent-to-unique
   recursive Subroutine mergesort(temp, Begin, Finish, list)
      ! 1st 3 arguments are input, 4th is output sorted list
      implicit none
      integer(kind=4), intent(inout) :: Begin, list(:), temp(:)
      integer(kind=4), intent(in) :: Finish
      integer(kind=4) :: Middle
      if (Finish - Begin < 2) then    !if run size =1
         return                   !it is sorted
      else
         ! split longer runs into halves
         Middle = (Finish + Begin)/2
         ! recursively sort both halves from list into temp
         call mergesort(list, Begin, Middle, temp)
         call mergesort(list, Middle, Finish, temp)
         ! merge sorted runs from temp into list
         call Merge(temp, Begin, Middle, Finish, list)
      end if
   end subroutine mergesort

   subroutine merge(list, Begin, Middle, Finish, temp)
      implicit none
      integer(kind=4), intent(inout) :: list(:), temp(:)
      integer(kind=4), intent(in) :: Begin, Middle, Finish
      integer(kind=4) :: kx, ky, kz
      ky = Begin
      kz = Middle
    !! While there are elements in the left or right runs...
      do kx = Begin, Finish - 1
       !! If left run head exists and is <= existing right run head.
         if (ky .lt. Middle .and. (kz .ge. Finish .or. list(ky) .le. list(kz))) then
            temp(kx) = list(ky)
            ky = ky + 1
         else
            temp(kx) = list(kz)
            kz = kz + 1
         end if
      end do

   end subroutine merge

   function unique(list)
    !! usage sortedlist=Unique(list)
      implicit none
      integer(kind=4) :: strt, fin, N
      integer(kind=4), intent(inout) :: list(:)
      integer(kind=4), allocatable :: unique(:), work(:)
      logical, allocatable :: mask(:)
      ! sort
      work = list; strt = 1; N = size(list); fin = N + 1
      call mergesort(work, strt, fin, list)
      ! cull duplicate indices
      allocate (mask(N)); 
      mask = .false.
      mask(1:N - 1) = list(1:N - 1) == list(2:N)
      unique = pack(list,.not. mask)

   end function unique

   function cumsum(a) result(b)
      ! cumulative sum of integer array a
      integer, intent(in) :: a(:)
      integer :: b(size(a) + 1)
      integer :: i
      b(1) = 0
      b(2:) = a
      do i = 2, size(a) + 1
         b(i) = b(i - 1) + a(i - 1)
      end do
   end function cumsum

   subroutine shuffle(a)
      implicit none
      integer, intent(inout) :: a(:)
      integer :: i, randpos, temp
      real(8) :: r
      do i = size(a), 2, -1
         r = grnd()
         randpos = int(r*i) + 1
         temp = a(randpos)
         a(randpos) = a(i)
         a(i) = temp
      end do
   end subroutine shuffle

   subroutine get_subsample(idxs, max_sample, subidxs)

      ! return a random subsample of data indices

      integer, intent(in) :: idxs(:), max_sample
      integer, allocatable, intent(out) :: subidxs(:)
      integer :: i, num_samples, x(size(idxs))

      num_samples = size(idxs)

      if (num_samples .lt. max_sample) then
         ! just return the initial idxs
         allocate (subidxs(size(idxs)))
         subidxs = idxs
         return
      end if

      ! else return max_sample length permutation of idxs
      allocate (subidxs(max_sample))
      x = [(i, i=1, num_samples)]
      call shuffle(x)
      subidxs = idxs(x(1:max_sample))

   end subroutine get_subsample

end module subs
