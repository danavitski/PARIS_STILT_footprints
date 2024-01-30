!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNFLX           SUN FLuX incident solar radiation at sfc
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SUN FLUX RETURNS THE INCIDENT SOLAR
!   IRRADIATION AT THE SURFACE BASED UPON THE AVERAGE RH AND SOLAR
!   ELEVATION ANGLE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 02 Mar 1998 - RRD (corrected 80. limit on Rh)
!
! USAGE:  CALL SUNFLX(NLVL,SEA,QQ,SWF,TR)
!   INPUT ARGUMENT LIST:
!     NLVL  - int number of vertical levels in profile
!     SEA   - real      sine of the solar elevation angle
!     QQ    - real      RH fraction profile
!   OUTPUT ARGUMENT LIST:
!     SWF   - real      incident short wave flux (w/m2)
!     TR    - real      transmissivity
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: sunflx.f90,v 1.3 2005-12-14 17:06:00 tnehrkor Exp $
!
!$$$

      SUBROUTINE SUNFLX(NLVL,SEA,QQ,SWF,TR)

      IMPLICIT REAL*8 (A-H,O-Z)

!     relative humidity fraction profile
      REAL*8 QQ(NLVL)

!     critical transmissivity    cloud top solar const (w/m2)
      DATA TCRIT/0.5/,           SOLC/1104.0/

!=>find max RH

      RHMAX=0.0
      KMAX = 1
      DO K=1,NLVL
         IF(QQ(K).GT.RHMAX)THEN
            KMAX=K
            RHMAX=QQ(K)
         END IF
      END DO

!=>find 3 layer average about max

      RAVG=0.0
      RCNT=0.0
      DO K=MAX(1,KMAX-1),MIN(NLVL,KMAX+1)
         RAVG=RAVG+QQ(K)
         RCNT=RCNT+1.0
      END DO
      RAVG=DMAX1(DBLE(0.80),RAVG/RCNT)

!=>fractional cloud cover when Rh >0.8

      FCLD=(5.0*(RAVG-0.8))**2
      FCLD=DMAX1(DBLE(0.0),DMIN1(DBLE(1.),FCLD))

!=>transmissivity based upon critical RH (TR=0.5 when FCLD=1)

      TR=1.0-FCLD*(1.0-TCRIT)

!=>adjust solar constant for elevation angle

      SWF=DMAX1(DBLE(0.),SEA*SOLC*TR)

      RETURN
      END
