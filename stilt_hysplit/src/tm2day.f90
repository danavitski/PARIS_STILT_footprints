!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TM2DAY           TiMe2DAY converts min since 1970 to time
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIME2DAY CONVERTS ACCUMULATED MINUTES SINCE BEGINNING OF
!   1970 TO TIME AS YEAR, MONTH, DAY, HOUR, MINUTE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Apr 1998 - RRD
!              07 May 1998 - RRD: y2k mod
!
! USAGE:  CALL TM2DAY(MACM,IY,IM,ID,IH,MN)
!   INPUT ARGUMENT LIST:
!     MACM        - int accumulated minutes
!   OUTPUT ARGUMENT LIST:
!     IY,IM,ID,IH,MN    - int time
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: tm2day.f90,v 1.3 2005-12-14 17:06:00 tnehrkor Exp $
!
!$$$

      SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER   IY,IM,ID,IH,MN,MACC,NADPM(12),NADPY(80)

!     accumulated days table incorporates years 1970 through 2049
      DATA NADPY/                                                       &
     &   365,  730, 1096, 1461, 1826, 2191, 2557, 2922, 3287, 3652,     &
     &  4018, 4383, 4748, 5113, 5479, 5844, 6209, 6574, 6940, 7305,     &
     &  7670, 8035, 8401, 8766, 9131, 9496, 9862,10227,10592,10957,     &
     & 11323,11688,12053,12418,12784,13149,13514,13879,14245,14610,     &
     & 14975,15340,15706,16071,16436,16801,17167,17532,17897,18262,     &
     & 18628,18993,19358,19723,20089,20454,20819,21184,21550,21915,     &
     & 22280,22645,23011,23376,23741,24106,24472,24837,25202,25567,     &
     & 25933,26298,26663,27028,27394,27759,28124,28489,28855,29220/

!     default number of accumulated days in each month (non leap year)
      DATA NADPM/0,31,59,90,120,151,181,212,243,273,304,334/

!     compute four digit year
      KK=1
      DO WHILE (NADPY(KK)*1440.LT.MACM)
         KK=KK+1
      END DO
      IYEAR=1970+(KK-1)

!     compute minutes in this year
      MACC=MACM-NADPY(KK-1)*1440

!     convert back to two digit year
      IY=MOD(IYEAR-1900,100)

!     current minute
      MHR=MACC/60
      MN=MACC-MHR*60

!     current hour
      MDA=MHR/24
      IH=MHR-MDA*24

!     current month and day
      DO K=1,12
         IDT=NADPM(K)

!        adjust accumulated days for leap year
         IF(MOD(IYEAR,4).EQ.0.AND.K.GT.2)IDT=IDT+1

         IF(IDT.LE.MDA)THEN
            MMO=IDT
            IM=K
         END IF
      END DO
      ID=MDA-MMO+1

!     check for end-of-year
      IF(ID.GT.31)THEN
         ID=ID-31
         IM=1
         IY=MOD(IY+1,100)
      END IF

      RETURN
      END
