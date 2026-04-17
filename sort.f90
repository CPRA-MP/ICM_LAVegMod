! This is an implementation of Shell Sort and is sourced from:
!   Paul E. Black and Art S. Kagel, "Shell sort", in Dictionary of Algorithms and Data Structures [online], Paul E. Black, ed. 6 April 2023. (accessed 12 April 2024) Available from: https://www.nist.gov/dads/HTML/shellsort.html
!
!   Fortran code below originally published by the Royal Statistical Society in the journal "Applied Statistics". 
!   Shell Sort was published as algorithm 304.8 in APPL.STATIST. (1996), VOL.45, NO.3.
!   Available online at https://lib.stat.cmu.edu/apstat/304

MODULE sort

CONTAINS

SUBROUTINE SHELLSORT (X, N)
! Sorts the N values stored in array X in ascending order

      IMPLICIT NONE
      
      INTEGER N
      REAL X(N)

      INTEGER I, J, INCR
      REAL TEMP

      INCR = 1

! Loop : calculate the increment
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

! Loop : Shell sort
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20

      RETURN
      END
      
END MODULE sort