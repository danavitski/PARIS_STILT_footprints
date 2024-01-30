!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFPAR           PUFf to PARticle conversion routine
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:99-08-27
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   IS USED TO CONVERT PUFFS TO PARTICLES WHEN THE PUFF SIZE GROWS TO
!   AN ARBITRARY RADIUS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 27 Aug 1999 (RRD)
!
! USAGE:  CALL PUFPAR(KPM,SIGH,SIGV,SIGX,HDWP,PGRD,NUMPAR)
!   INPUT ARGUMENT LIST:
!     KPM   - int total number of puffs or particles
!     NUMPAR      - int maximun number of particles followed
!   OUTPUT ARGUMENT LIST:
!     SIGH,SIGV   - real      array horiz (meters) and vert sigma (sigma)
!     SIGX  - real  array alongwind velocity variance
!     HDWP  - int array Horizontal distribution within pollutant
!     PGRD  - int array meteorological grid of puff position
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: pufpar.f90,v 1.4 2007-02-16 17:53:48 tnehrkor Exp $
!
!$$$

      SUBROUTINE PUFPAR(KPM,SIGH,SIGV,SIGX,HDWP,PGRD,NUMPAR)

      use module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'

!     horizontal, vertical, and alongwind velocity variances
      REAL*8 SIGH(MAXPAR),SIGV(MAXPAR),SIGX(MAXPAR)
!     puff distribution, met grid
      INTEGER HDWP(MAXPAR),PGRD(MAXPAR)

!      COMMON /GBLGRD/ GRID, DREC, FILE

!     sigma-h in meters at which puff will be converted to particle
      DATA PDIST/25000.0/

!     go through particles again to split if needed
      DO KP=1,KPM

!        check for on-grid and puff particle
         IF(PGRD(KP).GT.0.AND.HDWP(KP).GT.0)THEN

!           determine if large enough for conversion
            IF(SIGH(KP).GE.PDIST)THEN

!              velocity variances all set to zero
               SIGH(KP)=0.0
               SIGV(KP)=0.0
               SIGX(KP)=0.0

!              distribution becomes 3d particle
               HDWP(KP)=0

!           distance test
            END IF

!        puff test
         END IF

!     particle loop
      END DO

      RETURN
      END
