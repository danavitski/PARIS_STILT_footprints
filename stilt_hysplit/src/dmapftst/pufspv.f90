!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSPV           vertical PUFf SPLitting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SPLITTING VERTICAL OCCURS WHEN PUFF EXCEEDS VERTICAL LIMITS
!   IS THEN SPLIT INTO SMALLER UNITS.  VERTICAL DISTRIBUTION IS ALWAYS
!   ASSUMED TO BE UNIFORM (TOP-HAT). NUMBER OF SPLITS IS DETERMINED
!   HOW MANY VERTICAL MODEL SIGMA LEVELS ARE ENCOMPASSED BY THE PUFF
!   THE VERTICAL PUFF STANDARD DEVIATION IS THEN EVENLY DIVIDED BY THE
!   NUMBER AND THE NEW VERTICAL POSITION IS AT THE CENTRAL POINT.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,
!              SIGV,HDWP,PAGE,PTYP,PGRD,NSORT,NUMPAR)
!   INPUT ARGUMENT LIST:
!     KPM   - int total number of puffs or particles
!     NLVL  - int number of levels in subgrid
!     ZSG   - real      array internal model sigma levels
!     NUMPAR      - int maximum number of particles in calculation
!   OUTPUT ARGUMENT LIST:
!     MASS  - real      array mass of pollutant (arbitrary units)
!     XPOS,YPOS   - real      array puff center positions (grid units)
!     ZPOS  - real      array puff center height (sigma)
!     SIGH,SIGV   - real      array horiz (meters) and vert sigma (sigma)
!     HDWP  - int array Horizontal distribution within pollutant
!     PAGE  - int array pollutant age since release (min)
!     PTYP  - int array pollutant type index number
!     PGRD  - int array meteorological grid of puff position
!     NSORT - int array sorted array index values
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

      SUBROUTINE PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,          &
     &   SIGV,HDWP,PAGE,PTYP,PGRD,NSORT,NUMPAR)

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
      INCLUDE 'DEFSIZE.INC'

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER   PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8   MASS(MAXDIM,MAXPAR)
!     sigma levels
      REAL*8   ZSG(NLVL)
!     sorted index array
      INTEGER   NSORT(MAXPAR)

!     logical flag to determine split
      LOGICAL SPLIT

!     vertical extent
      DATA SIGR /1.54/

      KPT=KPM
!     go through particles again to split if needed
      DO KP=1,KPM

!     split at less frequent intervals after 24h
      SPLIT=.FALSE.
      KHRS=PAGE(KP)/60
      IF(KHRS.LE.24.OR.MOD(KHRS,3).EQ.0)SPLIT=.TRUE.

!     check for time, valid puff on grid, and  not particle
      IF(SPLIT.AND.PGRD(KP).GT.0.AND.                                   &
     &   (HDWP(KP).EQ.1.OR.HDWP(KP).EQ.2) )THEN

!        scan extent
         RADIUS=SIGR*SIGV(KP)

!        position at bottom and top of puff
         SGT=DMAX1(ZPOS(KP)-RADIUS,DBLE(0.0))
         SGB=DMIN1(ZPOS(KP)+RADIUS,DBLE(1.0))

!        compute vertical index
         KB=1
         KT=NLVL
         DO KK=1,(NLVL-1)
            IF(ZSG(KB+1).GE.SGB)KB=KB+1
            IF(ZSG(KT-1).LE.SGT)KT=KT-1
         END DO

!        check for vertical split (when extent = grid size) and
         NSPLIT=KT-KB
         IF(NSPLIT.GE.3.AND.(KPT+NSPLIT).LE.NUMPAR)THEN

!           top-hat new vertical sigma and layers
            DEPTH=(SGB-SGT)/NSPLIT
            SIGMA=DEPTH/SIGR/2.0
            FRAC=1.0/NSPLIT

!           new layer extent
            ZBOT=SGB
            ZTOP=ZBOT-DEPTH

!           simple bifurcation split across the old sigma
            DO KK=1,NSPLIT
               XPOS(KK+KPT)=XPOS(KP)
               YPOS(KK+KPT)=YPOS(KP)

!              new vertical position
               ZNEW=0.5*(ZBOT+ZTOP)
               ZPOS(KK+KPT)=DMAX1(DMIN1(ZNEW,DBLE(1.)),DBLE(0.))

               PAGE(KK+KPT)=PAGE(KP)
               SIGH(KK+KPT)=SIGH(KP)
               SIGV(KK+KPT)=SIGMA
               HDWP(KK+KPT)=HDWP(KP)
               PTYP(KK+KPT)=PTYP(KP)
               PGRD(KK+KPT)=PGRD(KP)
               NSORT(KK+KPT)=KK+KPT

!              mass is reduced to fraction of previous total
               MASS(1,KK+KPT)=FRAC*MASS(1,KP)
               MM=MAXDIM
               DO WHILE(MM.GT.1)
                  MASS(MM,KK+KPT)=FRAC*MASS(MM,KP)
                  MM=MM-1
               END DO

!              update layer
               ZBOT=ZTOP
               ZTOP=ZBOT-DEPTH

            END DO
            KPT=KPT+NSPLIT

!           set meteo grid to (0) for old elements
            PGRD(KP)=0

!        vertical split test
         END IF

!     particle loop
      END IF
      END DO

!     reset particle counter
      KPM=KPT

      RETURN
      END
