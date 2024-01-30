!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSPH           horizontal PUFf SPLitting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SPLITTING HORIZONTAL OCCURS WHEN PUFF EXCEEDS PREDEFINED LIMIT
!   THEN SPLIT INTO SMALLER UNITS.  GAUSSIAN PUFFS ARE SPLIT INTO
!   FIVE PUFFS, WHILE TOP-HAT PUFFS ALWAYS SPLIT INTO FOUR PUFFS.
!   MASS IS EQUALLY DIVIDED FOR TOP-HAT PUFFS WHILE FOR THE GAUSSIAN
!   PUFF 60. OF THE MASS IS PUT INTO THE CENTRAL PUFF AND 10. GOES
!   INTO EACH OF THE FOUR CORNER PUFFS.  THE NEW PUFF POSITIONS ARE
!   SET BY 0.5 OF THE HORIZONTAL EXTENT (RADIUS).
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 - RRD
!
! USAGE:  CALL PUFSPH(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,
!              SIGV,HDWP,PAGE,PTYP,PGRD,NSORT,NUMPAR)
!   INPUT ARGUMENT LIST:
!     KPM   - int total number of puffs or particles
!     HGD   - real      horizontal grid distance spacing (m)
!     NUMPAR      - int maximun number of particles followed
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

      SUBROUTINE PUFSPH(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,                   &
     &   SIGV,HDWP,PAGE,PTYP,PGRD,NSORT,NUMPAR)

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
      INCLUDE 'DEFSIZE.INC'
!     meteorology grid and file
      INCLUDE 'DEFGRID.INC'

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER   PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8   MASS(MAXDIM,MAXPAR)
!     split offsets
      REAL*8   XF(5), YF(5)
!     sorted index array
      INTEGER   NSORT(MAXPAR)

      COMMON /GBLGRD/ GRID, DREC, FILE

!     horizontal split distance factors (fraction of extent)
      DATA XF/0.5,0.0,-0.5,0.0,0.0/,  YF/0.0,-0.5,0.0,0.5,0.0/

      KPT=KPM
!     go through particles again to split if needed
      DO KP=1,KPM

!     check for on-grid and puff particle
      IF(PGRD(KP).GT.0.AND.HDWP(KP).GT.0)THEN

!        distribution type determines split procedure
         IF(HDWP(KP).EQ.1.OR.HDWP(KP).EQ.3)THEN
!           gaussian
            NSPLIT=5
            SIGR=3.0
            FRAC=0.10
         ELSEIF(HDWP(KP).EQ.2.OR.HDWP(KP).EQ.4)THEN
!           top-hat
            NSPLIT=4
            SIGR=1.54
            FRAC=0.25
         END IF
!        horizontal extent
         RADIUS=SIGR*SIGH(KP)

!        check for horizontal split (when extent = 1.54 grid size)
         HGD=GRID(PGRD(KP))%SIZE*1000.0
         IF(RADIUS.GE.(1.54*HGD).AND.(KPT+NSPLIT).LT.NUMPAR)THEN

            DO KK=1,NSPLIT
!              split into new puffs offset by the sigma
               XPOS(KK+KPT)=XPOS(KP)+XF(KK)*RADIUS/HGD
               YPOS(KK+KPT)=YPOS(KP)+YF(KK)*RADIUS/HGD

!              vertical position and time remain the same
               ZPOS(KK+KPT)=ZPOS(KP)
               PAGE(KK+KPT)=PAGE(KP)

!              horizontal sigma reduced by half
               SIGH(KK+KPT)=SIGH(KP)/2.0

!              vertical sigma is used as a marker for puffs
!              using mixed mode dispersion
               IF(HDWP(KP).EQ.3.OR.HDWP(KP).EQ.4)THEN
!                 2d puffs reset vertical turbulence after split
                  SIGV(KK+KPT)=0.0
               ELSE
!                 3d puffs or particles maintain value
                  SIGV(KK+KPT)=SIGV(KP)
               END IF

!              distribution, pollutant, and meteo grid
               HDWP(KK+KPT)=HDWP(KP)
               PTYP(KK+KPT)=PTYP(KP)
               PGRD(KK+KPT)=PGRD(KP)

!              save index value until new sort
               NSORT(KK+KPT)=KK+KPM

!              central point (5th) of Gaussian gets 60. of the mass
               IF(KK.EQ.5)FRAC=0.60

!              mass is reduced depending upon distribution
               MASS(1,KK+KPT)=FRAC*MASS(1,KP)
               MM=MAXDIM
               DO WHILE(MM.GT.1)
                  MASS(MM,KK+KPT)=FRAC*MASS(MM,KP)
                  MM=MM-1
               END DO
            END DO
            KPT=KPT+NSPLIT

!           set meteo grid to (0) for old elements
            PGRD(KP)=0

         END IF

!     particle loop
      END IF
      END DO

!     reset particle counter
      KPM=KPT

      RETURN
      END
