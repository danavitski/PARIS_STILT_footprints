!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPSUS           DEPosition reSUSpension of a pollutant
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION RESUSPENSION OF A POLLUTANT - ASSUME THAT RATIO (K)
!   OF POLLUTANT IN AIR (C=MASS/M3) TO SURFACE VALUE (S=MASS/M2) =
!   10^-6 DIVIDED BY DURATION OF DEPOSITION IN DAYS.  FOR SIMPLICITY
!   WE ASSUME DAYS ALWAYS = 1.  K ALSO DEFINED AS R/S DS/DT, WHERE R
!   IS RELATED TO THE ATMOSPHERIC RESISTENCE = 1/KU*   THEREFORE
!   THE RESUSPENSION FLUX = S K / R.  NOTE THAT MULTIPLE SPECIES
!   CONCENTRATION FILES WILL RESULTS IN INDEPENDENT PARTICLES FOR EACH
!   DEFINED POLLUTANT FROM ONLY THE FIRST CONCENTRATION GRID.  MULTI
!   GRID DEFINITIONS ARE NOT SUPPORTED IN THIS APPLICATION. NOTE THAT
!   THE METEO VARIABLE FRICTION VELOCITY IS PASSED THROUGH AS THE
!   VALUE COMPUTED IN ADVPNT, THEREFORE ONLY VALID IN A SMALL DOMAIN
!   STANDARD MODEL CONFIGURATION DOES NOT SUPPORT FULL-GRID METEO.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 - RRD
!
! USAGE:  CALL DEPSUS(INITD,KGM,NUMGRD,NUMTYP,DT,KPM,CSUM,
!              MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,
!              PTYP,PGRD,NSORT)
!   INPUT ARGUMENT LIST:
!     INITD - int initial puff/particle distribution
!     KGM   - int current meteorological grid number
!     NUMGRD      - int number of concentration grids
!     NUMTYP      - int number of pollutants
!     DT    - real      time step (min)
!     CSUM  - real      concentration array
!   OUTPUT ARGUMENT LIST:
!     KPM   - int total number of puffs or particles
!     MASS  - real      mass of pollutant (arbitrary units)
!     XPOS,YPOS   - real      puff center positions (grid units)
!     ZPOS  - real      puff center height (sigma)
!     SIGH,SIGV   - real      horiz (meters) and vert sigma (sigma)
!     SIGX  - real      x component turbulence for 3D particle model
!     HDWP  - int Horizontal distribution within pollutant
!     PAGE  - int pollutant age since release (min)
!     PTYP  - int pollutant type index number
!     PGRD  - int meteorological grid of puff position
!     NSORT - int sorted array index values
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: depsus.f90,v 1.6 2008-03-26 19:16:03 tnehrkor Exp $
!
!$$$

      SUBROUTINE DEPSUS(INITD,KGM,NUMGRD,NUMTYP,DT,KPM,CSUM,            &
     &   MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,NSORT)

      use module_defmeto
      use module_defgrid
      use module_defconc

      IMPLICIT REAL*8 (A-H,O-Z)

!     default array sizes
!     meteorological grid and file definitions
!      INCLUDE 'DEFGRID.INC'
!     concentration grid and pollutant definitions
!      INCLUDE 'DEFCONC.INC'
!     meteorological information
!      INCLUDE 'DEFMETO.INC'

!     master concentration array (x,y,z,grids,species)
      REAL*8 CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!     particle positions
      REAL*8 XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8 SIGH(MAXPAR),SIGV(MAXPAR),SIGX(MAXPAR)
!     puff age, distribution, pollutant index, meteo grid, sort index
      INTEGER PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR),      &
     &   NSORT(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8 MASS(MAXDIM,MAXPAR)

!      COMMON /GBLMET/ METO
!      COMMON /GBLCON/ CONC, DIRT
!      COMMON /GBLGRD/ GRID, DREC, FILE

!         vonKarman's,  radius
      DATA VONK/0.40/,  SIGR/1.54/

      KGC=0
!     find the first grid and level which has deposition defined
      DO KG=1,NUMGRD
         IF(KGC.EQ.0)THEN
!           determine loop indicies
            NXP=CONC(KG)%NUMB_LON
            NYP=CONC(KG)%NUMB_LAT
            NZP=CONC(KG)%LEVELS
            DO KL=1,NZP
               IF(CONC(KG)%HEIGHT(KL).EQ.0)THEN
                  KGC=KG
                  KGL=KL
               END IF
            END DO
         END IF
      END DO
      IF(KGC.EQ.0)RETURN

      DO KT=1,NUMTYP
!     resuspension must be defined for this pollutant
      IF(DIRT(KT)%DOSUS)THEN

         DO JJ=1,NYP
         DO II=1,NXP

            DEPAMT=CSUM(II,JJ,KGL,KT,KGC)
            IF(DEPAMT.GT.0.0)THEN

!              increment particle/puff counter
               KPM=KPM+1
               IF(KPM.GT.MAXPAR)THEN
                  KPM=MAXPAR
                  WRITE(30,*)'WARNING depsus: exceeding puff limit'
                  RETURN
               END IF
               NSORT(KPM)=KPM

!              convert position to meteorological grid units
               PLON=(II-1)*CONC(KGC)%DELT_LON+CONC(KGC)%X1Y1_LON
               PLAT=(JJ-1)*CONC(KGC)%DELT_LAT+CONC(KGC)%X1Y1_LAT

! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
               IF(GRID(KGC)%LATLON)THEN
                 CALL GBL2XY(KGC,PLAT,PLON, XP, YP)
               ELSE
                 CALL CLL2XY(GRID(KGC)%GBASE,PLAT,PLON, XP, YP, GRID(KGC)%proj)
               END IF

!              initial position always at ground
               XPOS(KPM)=XP
               YPOS(KPM)=YP
               ZPOS(KPM)=1.0

!              alongwind and vertical variances start at zero
               SIGX(KPM)=0.0
               SIGV(KPM)=0.0

!              initial distribution (see main for definitions)
               HDWP(KPM)=INITD
!              initial age at zero
               PAGE(KPM)=0
!              pollutant type definition
               PTYP(KPM)=KT
!              initial grid is the default startup grid from main
               PGRD(KPM)=KGM

               IF(INITD.EQ.0)THEN
!                 3D particle model starts with zero horizontal
                  SIGH(KPM)=0.0
               ELSE
!                 all other combinations assume grid-cell area
!                 source (111000 m / deg - squared)
                  AREA=1.2E+10*CONC(KGC)%DELT_LAT                       &
     &                        *CONC(KGC)%DELT_LON*DCOS(PLAT/57.3)
!                 compute sigma for uniform radius
                  SIGH(KPM)=DSQRT(AREA/3.14159)/SIGR
               END IF

!              pollutant flux - mass/m2-s
               QRATE=DEPAMT*VONK*METO%USTR*DIRT(KT)%SRATE

!              determine mass lost from surface
               QTOT=DMIN1(DEPAMT, 60.0*DT*QRATE)
               CSUM(II,JJ,KGL,KT,KGC)=DEPAMT-QTOT
               MASS(KT,KPM)=QTOT*AREA

            END IF

!        horizontal grid loop
         END DO
         END DO

!     pollutant types
      END IF
      END DO

      RETURN
      END
