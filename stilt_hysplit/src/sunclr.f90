!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNANG           SUN ANGle returns the solar angle
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SUN ANGLE - RETURNS THE SOLAR
!   ELEVATION ANGLE (90-ZENITH ANGLE) AND THE SINE OF THE ANGLE,
!   GIVEN THE LAT/LON AND TIME.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Apr 1997 - RRD
!
! USAGE:  CALL SUNCLR(JYR,JET,OLAT,OLON,EA,SEA)
!   INPUT ARGUMENT LIST:
!     JET   - int elapsed minutes since January 1st 1970
!     OLAT  - real      latitude in degrees and fraction (+ = north)
!     OLON  - real      longitude in degrees anb fraction (- = west)
!   OUTPUT ARGUMENT LIST:
!     SWF   - clear sky short wave radiative flux
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: sunclr.f90,v 1.3 2005-12-14 17:06:00 tnehrkor Exp $
!
!$$$
! CHG:(9/24/02) copied from sunang.f, to compute clear sky DSWF
      SUBROUTINE SUNCLR(JET,OLAT,OLON,SWF)

      IMPLICIT REAL*8 (A-H,O-Z)

!     radians per degree
      DATA RDPDG/0.0174532925199433d0/

! CHG(9/24/02)
!     cloud top solar const (w/m2)
      DATA SOLC/1104.0/

! CHG(9/24/02) get current year from JET
      CALL TM2DAY(JET,JYR,IM,ID,IH,MN)

!     adjust for leap year
      DTYR=365.242
      IF(MOD(JYR+1900,4).EQ.0)DTYR=366.242

!     compute elapsed time for this year
      CALL TM2MIN(JYR,1,1,0,0,JET0)
      JET1=JET-JET0

!     derive fractional part of year
      DF=360.0*JET1/(DTYR*24.0*60.0)

      SIG=279.9348+DF+1.914827*SIN(RDPDG*DF)-0.079525*COS(RDPDG*DF)     &
     &    +0.019938*SIN(2.0*RDPDG*DF)-0.00162*COS(2.0*RDPDG*DF)

!     solar declination at noon
      SD=SIN(23.4438333*RDPDG)*SIN(RDPDG*SIG)
      D=ASIN(SD)

!     time of meridian passage (solar noon time)
      SNOON=12.0+0.12357*SIN(RDPDG*DF)-0.004289*COS(RDPDG*DF)           &
     &      +0.153809*SIN(2.0*RDPDG*DF)+0.060783*COS(2.0*RDPDG*DF)

!     current hour
      HR=AMOD(JET1/60.0,24.0)

!     hour angle at current time
      HA=(HR-SNOON)*15.0+OLON

!     sine of the solar elevation angle
      SEA=COS(HA*RDPDG)*COS(OLAT*RDPDG)*COS(D)+SIN(OLAT*RDPDG)*SD

!     solar elevation angle in degrees
      EA=ASIN(SEA)/RDPDG

! CHG(9/24/02) get clear sky swf
      SWF=DMAX1(DBLE(0.),SEA*SOLC)

      RETURN
      END