!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TM2MIN           TiMe2MINute converts time to min since 1970
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIME2MINUTE CONVERTS TIME AS EXPRESSED IN DATE FORM:
!   YEAR, MONTH, DAY, HOUR, MINUTE  TO ACCUMULATED MINUTES
!   SINCE THE BEGINNING OF 1970
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Apr 1998 - RRD
!              07 May 1998 - RRD - y2k mod
!
! USAGE:  CALL TM2MIN(IY,IM,ID,IH,MN,MACC)
!   INPUT ARGUMENT LIST:
!     IY,IM,ID,IH,MN    - int date and time
!   OUTPUT ARGUMENT LIST:
!     MACC        - int accumulated minutes since 1 Jan 1970
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: tm2min.f90,v 1.3 2005-12-14 17:06:00 tnehrkor Exp $
!
!$$$

      SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)

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
      IYEAR=MOD(IY+99,100)+1901
      IF(IYEAR.LT.1970)IYEAR=IYEAR+100

!     number of accumulated days until this year
      IDT=NADPY(IYEAR-1970)

!     add accumulated days for this year
      IDT=IDT+NADPM(IM)+(ID-1)

!     adjust accumulated days for leap year
!     does not account for centuries where mod(iyear,400)=0
      IF(MOD(IYEAR,4).EQ.0.AND.IM.GT.2)IDT=IDT+1

!     convert to minutes
      MACC=MN+(IH+IDT*24)*60

      RETURN
      END
