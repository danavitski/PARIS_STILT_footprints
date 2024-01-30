!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVSRT           ADVection SoRTting to minimize data I/O
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION SORTTING TO MINIMIZE DATA I/O WHEN METEOROLOGICAL
!   SUBGRID IS DEFINED.  MAY BE CALLED BEFORE EACH ADVECTION CYCLE
!   WHEN THE SUB-GRID IS SMALLER THAN THE FULL METEOROLOGICAL GRID.
!   ALL ELEMENTS THAT ARE ON THE CURRENT SUB-GRID ARE GROUPED TOGETHER
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 23 Jul 1998 - RRD
!
! USAGE:  CALL ADVSRT(KPM,LX1,LY1,XPOS,YPOS,NSORT)
!   INPUT ARGUMENT LIST:
!     KPM   - int  number of points
!     LX1,LY1     - int  lower left corner of current mete subgrid
!     XPOS,YPOS   - real puff center positions (grid units)
!   OUTPUT ARGUMENT LIST:
!     NSORT - int  sortted array by position
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: advsrt.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

      SUBROUTINE ADVSRT(KPM,LX1,LY1,XPOS,YPOS,NSORT)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)


!     particle positions
      REAL*8 XPOS(MAXPAR),YPOS(MAXPAR)
!     sort index
      INTEGER NSORT(MAXPAR)

!     need at least three elements to sort
      IF(KPM.LE.2)RETURN
!     exit if subgrid not yet set
      IF(LX1.LE.0.OR.LY1.LE.0)RETURN

!     index for points within subgrid (eventually equals KPM)
      KPT=0

!     set scan range for subgrid
      NXS2=(NXM-1)/2-2
      NYS2=(NYM-1)/2-2

!     default initial position to center of current subgrid
      XCNT=LX1+NXS2
      YCNT=LY1+NYS2

!     zero out index array to indicate unsorted particles
      DO KK=1,KPM
         NSORT(KK)=0
      END DO

!     loop through particle array until none left
      DO WHILE (KPT.LT.KPM)

!        reset the subgrid center position each pass
         XNEW=-1.0
         YNEW=-1.0

!        determine if position within range of subgrid
         DO KK=1,KPM
!           only test particles not previously identified
            IF(NSORT(KK).EQ.0)THEN
               IF(XPOS(KK).GT.XCNT-NXS2.AND.XPOS(KK).LT.XCNT+NXS2.AND.  &
     &            YPOS(KK).GT.YCNT-NYS2.AND.YPOS(KK).LT.YCNT+NYS2)THEN
!                 identify elements within current subgrid
                  KPT=KPT+1
                  NSORT(KPT)=KK
               ELSE
!                 save position of first particle not in array
!                 which will be the center of the next subgrid
                  IF(XNEW.LT.0.0.AND.YNEW.LT.0.0)THEN
                     XNEW=XPOS(KK)
                     YNEW=YPOS(KK)
                  END IF
               END IF
            END IF
         END DO

!        set new sort center point to last outside
         XCNT=XNEW
         YCNT=YNEW

!     while loop
      END DO

      RETURN
      END
