!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TMPLUS           TiMePLUS is used to add time
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIMEPLUS CONVERTS A DATE PLUS HOURS TO THE NEW DATE
!
! PROGRAM HISTORY LOG:
!   Last Revised: 14 Feb 1997 (RRD)
!                 06 Aug 1999 (RRD) - century should be 00 not 100
!
! USAGE:  CALL TMPLUS(IY,IM,ID,IH,IC)
!   INPUT ARGUMENT LIST:
!     IY,IM,ID    - int current date
!     IH    - int current hour
!     IC    - int increment in hours
!   OUTPUT ARGUMENT LIST:
!     IY,IM,ID    - int new date
!     IH    - int new hour
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE TMPLUS(IY,IM,ID,IH,IC)

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER NDM(12)
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/

!     new century year should be 00 not 100
      IY=MOD(IY,100)

!     check for leap year
      NDM(2)=28
      IYEAR=MOD(IY+99,100)+1901
      IF(IY.EQ.0)THEN
         IF(MOD(IYEAR,400).EQ.0)NDM(2)=29
      ELSE
         IF(MOD(IYEAR,4).EQ.0)NDM(2)=29
      END IF

!     current time (min) incremented by IC (+/-)
      IF(IC.EQ.0)RETURN

!     loop through by hour
      DO LOOP=1,IABS(IC)
         IH=IH+SIGN(1,IC)

!        forward clock
         IF(IC.GT.0)THEN
            IF(IH.GE.24)THEN
               IH=IH-24
               ID=ID+1
               IF(ID.GT.NDM(IM))THEN
                  ID=1
                  IM=IM+1
                  IF(IM.GT.12)THEN
                     IM=1
                     IY=IY+1
                     IY=MOD(IY,100)
                  END IF
               END IF
            END IF

!        backward clock
         ELSE
            IF(IH.LT.0)THEN
               IH=24+IH
               ID=ID-1
               IF(ID.LT.1)THEN
                  IM=IM-1
                  IF(IM.LT.1)THEN
                     IM=12
                     IY=IY-1
                  END IF
                  ID=NDM(IM)
               END IF
            END IF
         END IF

      END DO
      RETURN
      END
