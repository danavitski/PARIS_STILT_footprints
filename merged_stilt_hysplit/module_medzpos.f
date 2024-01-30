module module_medzpos

  private
  
  public :: compute_medzpos, compute_pctzpos

contains

  subroutine compute_medzpos(zpos, dyn_veg, pgrd, ptyp, medzpos, kpm, numtyp)
    ! compute median particle height (in sigma coordinates) for each
    ! pollutant/receptor (ptyp)
    implicit none
    integer, intent(in) :: kpm, numtyp
    real, intent(in) :: zpos(kpm)
    integer, intent(in) :: pgrd(kpm), ptyp(kpm)
    logical, intent(in) :: dyn_veg(numtyp)
    real, intent(out) :: medzpos(numtyp)
    ! local variables
    integer :: i, jtyp, kp, kpmtyp, nmed, ipos1, ipos2
    real :: zpostyp(kpm)
    integer :: ipostyp(kpm)
    
    medzpos(:)=0. !return 0 if there are no particles for this ptyp, or dyn_veg is off
    do jtyp=1,numtyp
       if (dyn_veg(jtyp)) then
          ! store zpos for ptyp == jtyp
          kpmtyp=0
          do kp=1,kpm
             if (pgrd(kp) .ne. 0 .and. ptyp(kp) .eq. jtyp) then
                kpmtyp=kpmtyp+1
                zpostyp(kpmtyp)=zpos(kp)
             end if
          end do
          if (kpmtyp .gt. 0) then
             ! find sort order of zpostyp
             call quicksort(kpmtyp,zpostyp,ipostyp)
             ! find median value
             medzpos(jtyp)=0.
             nmed=0
             ipos1=(kpmtyp+1)/2
             ipos2=(kpmtyp+2)/2
             do i=ipos1,ipos2
                nmed=nmed+1
                medzpos(jtyp)=medzpos(jtyp)+zpostyp(ipostyp(i))
             end do
             if (nmed .lt. 1 .or. nmed .gt. 2) &
                  stop 'internal logic error in compute_medzpos, nmed <1 or >2'
             medzpos(jtyp) = medzpos(jtyp)/nmed
          end if
       end if
    end do
  end subroutine compute_medzpos

  subroutine compute_pctzpos(zpos, dyn_veg, pgrd, ptyp, medzpos, kpm, numtyp, pct)
    ! compute percentile (pct) particle height (in sigma coordinates) for each
    ! pollutant/receptor (ptyp)
    implicit none
    integer, intent(in) :: kpm, numtyp, pct
    real, intent(in) :: zpos(kpm)
    integer, intent(in) :: pgrd(kpm), ptyp(kpm)
    logical, intent(in) :: dyn_veg(numtyp)
    real, intent(out) :: medzpos(numtyp)
    ! local variables
    integer :: i, jtyp, kp, kpmtyp, nmed, ipos1, ipos2
    real :: xpos, zpostyp(kpm)
    integer :: ipostyp(kpm)
    
    medzpos(:)=0. !return 0 if there are no particles for this ptyp, or dyn_veg is off
    do jtyp=1,numtyp
       if (dyn_veg(jtyp)) then
          ! store zpos for ptyp == jtyp
          kpmtyp=0
          do kp=1,kpm
             if (pgrd(kp) .ne. 0 .and. ptyp(kp) .eq. jtyp) then
                kpmtyp=kpmtyp+1
                zpostyp(kpmtyp)=zpos(kp)
             end if
          end do
          if (kpmtyp .gt. 0) then
             ! find sort order of zpostyp
             call quicksort(kpmtyp,zpostyp,ipostyp)
             ! find percentile (pct) value
             medzpos(jtyp)=0.
             nmed=0
             xpos = pct*(kpmtyp-1)/100.
             ipos1 = 1 + int(xpos)
             ipos2 = 1 + int(xpos+0.5)
             ipos1 = min(kpmtyp,max(1,ipos1))
             ipos2 = min(kpmtyp,max(1,ipos2))
             do i=ipos1,ipos2
                nmed=nmed+1
                medzpos(jtyp)=medzpos(jtyp)+zpostyp(ipostyp(i))
             end do
             if (nmed .lt. 1 .or. nmed .gt. 2) &
                  stop 'internal logic error in compute_pctzpos, nmed <1 or >2'
             medzpos(jtyp) = medzpos(jtyp)/nmed
          end if
       end if
    end do
  end subroutine compute_pctzpos

  SUBROUTINE quicksort(n,x,ind)

    ! source code downloaded 12/18/2015 from 
    !  http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/sorc/global_enkf.fd/quicksort.f90
    ! modified to: replace real(r_kind) by real; made ind intent(out)

    IMPLICIT NONE

    real, INTENT(IN)  :: x(n)
    INTEGER, INTENT(OUT)   :: ind(n)
    INTEGER, INTENT(IN)    :: n

    !***************************************************************************

    !                                                         ROBERT RENKA
    !                                                 OAK RIDGE NATL. LAB.

    !   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL 
    ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
    ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
    ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
    ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
    ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
    ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
    ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
    ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
    ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
    ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
    ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
    ! UNSORTED PORTION.

    ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

    !                      X - VECTOR OF LENGTH N TO BE SORTED.

    !                    IND - VECTOR OF LENGTH >= N.

    ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

    ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
    !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
    !                          X IS DEFINED BY Y(I) = X(IND(I)).

    !*********************************************************************

    ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.
    ! (OK for N up to about a billon)

    !*********************************************************************

    INTEGER   :: iu(21), il(21)
    INTEGER   :: m, i, j, k, l, ij, it, itt, indx
    real      :: r
    real      :: t

    ! LOCAL PARAMETERS -

    ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
    !            INDICES OF PORTIONS OF THE ARRAY X
    ! M =      INDEX FOR IU AND IL
    ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
    ! K,L =    INDICES IN THE RANGE I,...,J
    ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
    ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
    ! INDX =   TEMPORARY INDEX FOR X
    ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
    ! T =      CENTRAL ELEMENT OF X

    IF (n <= 0) RETURN

    ! INITIALIZE IND, M, I, J, AND R

    DO  i = 1, n
       ind(i) = i
    END DO
    m = 1
    i = 1
    j = n
    r = .375

    ! TOP OF LOOP

20  IF (i >= j) GO TO 70
    IF (r <= .5898437) THEN
       r = r + .0390625
    ELSE
       r = r - .21875
    END IF

    ! INITIALIZE K

30  k = i

    ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

    ij = i + r*(j-i)
    it = ind(ij)
    t = x(it)

    ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(i)
    IF (x(indx) > t) THEN
       ind(ij) = indx
       ind(i) = it
       it = indx
       t = x(it)
    END IF

    ! INITIALIZE L

    l = j

    ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(j)
    IF (x(indx) >= t) GO TO 50
    ind(ij) = indx
    ind(j) = it
    it = indx
    t = x(it)

    ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(i)
    IF (x(indx) <= t) GO TO 50
    ind(ij) = indx
    ind(i) = it
    it = indx
    t = x(it)
    GO TO 50

    ! INTERCHANGE ELEMENTS K AND L

40  itt = ind(l)
    ind(l) = ind(k)
    ind(k) = itt

    ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
    !   NOT LARGER THAN T

50  l = l - 1
    indx = ind(l)
    IF (x(indx) > t) GO TO 50

    ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60  k = k + 1
    indx = ind(k)
    IF (x(indx) < t) GO TO 60

    ! IF K <= L, INTERCHANGE ELEMENTS K AND L

    IF (k <= l) GO TO 40

    ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
    !   ARRAY YET TO BE SORTED

    IF (l-i > j-k) THEN
       il(m) = i
       iu(m) = l
       i = k
       m = m + 1
       GO TO 80
    END IF

    il(m) = k
    iu(m) = j
    j = l
    m = m + 1
    GO TO 80

    ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70  m = m - 1
    IF (m == 0) RETURN
    i = il(m)
    j = iu(m)

80  IF (j-i >= 11) GO TO 30
    IF (i == 1) GO TO 20
    i = i - 1

    ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90  i = i + 1
    IF (i == j) GO TO 70
    indx = ind(i+1)
    t = x(indx)
    it = indx
    indx = ind(i)
    IF (x(indx) <= t) GO TO 90
    k = i

100 ind(k+1) = ind(k)
    k = k - 1
    indx = ind(k)
    IF (t < x(indx)) GO TO 100

    ind(k+1) = it
    GO TO 90
  END SUBROUTINE quicksort

end module module_medzpos
