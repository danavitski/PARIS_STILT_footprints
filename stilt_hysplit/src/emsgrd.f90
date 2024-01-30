!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSGRD           EMiSsion GRiDded input starts puff/particle
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION GRIDDED INPUT - STARTS A NEW PUFF OR PARTICLE FROM POINT
!   DEFINED IN AN INPUT FILE. OTHER RELEASE CRITERIA ARE THE SAME.
!   INPUT FILE READ IN SUBROUTINE EMSINP.
!   EMISSIONS CAN BE STARTED AGAIN AT QCYCLE INTERVALS.
!   THE ROUTINE IS CALLED EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 26 Apr 1998 (RRD)
!                 22 Dec 1998 (RRD) - default comments for chemistry
!                 22 Jul 1999 (RRD) - variable name change
!
! USAGE:  CALL EMSGRD(NUMTYP,KPM,INITD,DT,JET,KG,
!              NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,
!              QCYCLE,NUMPAR,TBAR)
!   INPUT ARGUMENT LIST:
!     NUMTYP    - int   number of pollutant types
!     KPM       - int   number of puffs or particles
!     INITD     - int   initial distribution type
!     DT        - real  time step (min)
!     JET       - int   current elapsed time (min)
!     KG        - int   default initial grid index
!     QCYCLE    - real  optional emission cycle time in hours
!     NUMPAR    - int   maximum number of particles permitted
!     TBAR      - real  optional temperature array on conc grid
!   OUTPUT ARGUMENT LIST:
!     NSORT     - int   array index of sorted elements
!     MASS      - real  array mass of pollutant (arbitrary units)
!     XPOS,YPOS - real  array horizontal position (grid units)
!     ZPOS      - real  array puff center height (sigma)
!     SIGH,SIGV - real  array horiz (meters) and vert sigma (sigma)
!     SIGX      - real  array x component turbulence (particle model)
!     HDWP      - int   array Horizontal distribution within pollutant
!     PAGE      - int   array pollutant age since release (min)
!     PTYP      - int   array pollutant type index number
!     PGRD      - int   array meteorological grid of puff position
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: emsgrd.f90,v 1.7 2008-03-26 19:16:04 tnehrkor Exp $
!
!$$$

! JCL: add BACK as argument
      SUBROUTINE EMSGRD(NUMTYP,KPM,INITD,DT,JET,KG,                     &
     &   NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,  &
     &   QCYCLE,NUMPAR,TBAR,BACK)

!      SUBROUTINE EMSGRD(NUMTYP,KPM,INITD,DT,JET,KG,
!     :   NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,
!     :   QCYCLE,NUMPAR,TBAR)

      use module_defgrid
      use module_defconc
      use module_defspot

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology file
!      INCLUDE 'DEFGRID.INC'
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'
!     multiple source information
!      INCLUDE 'DEFSPOT.INC'

!     number of emission grid points
      PARAMETER (MQLAT=150, MQLON=150, MQVAL=4, MQHRS=24)
      REAL*8 QAREA(MQLON,MQLAT,MQVAL,MQHRS)
      CHARACTER POLID(MQVAL)*4

! JCL:
      LOGICAL BACK

!     particle positions
      REAL*8 XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8 SIGH(MAXPAR),SIGV(MAXPAR),SIGX(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8 MASS(MAXDIM,MAXPAR)
!     pollutant sort index
      INTEGER NSORT(MAXPAR)
!     optional temperature array required for chemistry applications
!     REAL*8 TBAR(MAXXP,MAXYP)
!     other flags
      LOGICAL EMIT

!      COMMON /GBLSPT/ SPOT
!      COMMON /GBLCON/ CONC, DIRT
!      COMMON /GBLGRD/ GRID, DREC, FILE
      COMMON /QARRAY/ KLAT, KLON, DLAT, DLON, QAREA, POLID
      SAVE /QARRAY/

!     top-hat radius
      DATA SIGR/1.5/

!==>compute current hour
      IH=MOD(JET/60,24)+1

!==>check to determine if emissions required for any pollutant

      EMIT=.FALSE.
      DO KK=1,NUMTYP
! JCL: differentiate between when BACK is T & F
         IF(JET.GE.DIRT(KK)%START%MACC.AND.                             &
     &      JET.LT.DIRT(KK)%START%MACC+DIRT(KK)%QHRS*60                 &
     &      .AND.(.NOT.BACK))THEN
           EMIT=.TRUE.
! JCL: when BACK is T: notice that QHRS is subtracted, not added
         ELSEIF(JET.LE.DIRT(KK)%START%MACC.AND.                         &
     &      JET.GT.DIRT(KK)%START%MACC-DIRT(KK)%QHRS*60                 &
     &      .AND.BACK)THEN
           WRITE(45,*)'EMIT! BACK IS ',BACK
           EMIT=.TRUE.
         END IF
      END DO
      IF(.NOT.EMIT)RETURN
      EMIT=.FALSE.

!==>set the number of units to emit

      IF(INITD.EQ.1.OR.INITD.EQ.2)THEN
!        gaussian or top-hat emissions
         NPAR=1
      ELSE
         QMAX=0.0
         NRPTS=KLAT*KLON
!        find maximum emission duration
         DO KK=1,NUMTYP
! JCL:      need to take ABS of DT since DT can be < 0
            QMAX=MAX(QMAX,ABS(DT)/60.0,DIRT(KK)%QHRS)
!           QMAX=MAX(QMAX,DT/60.0,DIRT(KK)%QHRS)
         END DO
!        particle emission rate: number particles per hour per source
         NPHR=NUMPAR/QMAX/NRPTS

! JCL:   need to take ABS since DT can be < 0 when BACK is T
!        particle emissions per time step
         NPAR=MAX(1, ABS(INT(NPHR*DT/60.0)))
!         NPAR=MAX(1, INT(NPHR*DT/60.0))
      END IF

!==>loop through the number of independent source locations

      DO JJ=1,KLAT
      DO II=1,KLON

!     check for non-zero emissions
      QSUM=0.0
      DO KK=1,MQVAL
         QSUM=QSUM+QAREA(II,JJ,KK,IH)
      END DO
      IF(QSUM.GT.0.0)THEN

!     compute position from lower left corner point (1,1 point)
!     assume emissions at center of grid cell (from SW corner)
      CLAT=(JJ-1)*DLAT+SPOT(1)%OLAT+DLAT/2.0
      CLON=(II-1)*DLON+SPOT(1)%OLON+DLON/2.0

! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
!     convert to grid units
      IF(GRID(KG)%LATLON)THEN
        CALL GBL2XY(KG,CLAT,CLON,XP,YP)
      ELSE
        CALL CLL2XY(GRID(KG)%GBASE,CLAT,CLON,XP,YP, GRID(KG)%proj)
      END IF

!     check for location relative to meteo grid
      IF(XP.GT.2.AND.XP.LT.(GRID(KG)%NX-1).AND.                         &
     &   YP.GT.2.AND.YP.LT.(GRID(KG)%NY-1))THEN

!     multiple emissions at each cell only for particles
      DO KP=1,NPAR

!        increment particle/puff counter
         KPM=KPM+1
         IF(KPM.GT.MAXPAR)THEN
            KPM=MAXPAR
            WRITE(30,*)'Warning: emsgrd - exceeding puff limit'
            RETURN
         END IF
         NSORT(KPM)=KPM

!        initial position
         XPOS(KPM)=XP
         YPOS(KPM)=YP

         IF(INITD.EQ.0)THEN
            SIGH(KPM)=0.0
         ELSE
!           assume area source (111000 m / deg - squared)
            AREA=1.2E+10*DLAT*DLON*COS(CLAT/57.3)
!           compute sigma for uniform radius
            SIGH(KPM)=DSQRT(AREA/3.14159)/SIGR
         END IF

!        initial depth defined from release height
!        input height assumed to define layer from ground
         DELZ=1.0-SPOT(1)%ZP

         IF(INITD.EQ.1.OR.INITD.EQ.2)THEN
!           puff variance set to layer depth
            SIGV(KPM)=ABS(DELZ/SIGR/2.0)
            ZPOS(KPM)=1.0-DELZ/2.0
         ELSE
            SIGV(KPM)=0.0
!           particles get distributed in the layer
            ZPOS(KPM)=1.0-KP*DELZ/NPAR
         END IF

!        variances start at zero
         SIGX(KPM)=0.0

!        initial distribution (see main for definitions)
         HDWP(KPM)=INITD
!        initial age at zero
         PAGE(KPM)=0
!        variable not used in this configuration (all in single element)
         PTYP(KPM)=MAXDIM
!        initial grid is the default startup grid from main
         PGRD(KPM)=KG

!==>loop through pollutants

         DO KK=1,NUMTYP

            KT=0
!           match pollutant to release with emissions table
            DO KD=1,MQVAL
               IF(DIRT(KK)%IDENT.EQ.POLID(KD))KT=KD
            END DO
!           at least one pollutant match must be found in table
            IF(.NOT.EMIT.AND.KT.NE.0)EMIT=.TRUE.

! JCL:   add BACK in condition for pollutant requiring a start
!        note reversal of sign--subtracting instead of adding
!        check if this pollutant requires a start
         IF(((.NOT.BACK).AND.(JET.GE.DIRT(KK)%START%MACC).AND.          &
     &      (JET.LT.DIRT(KK)%START%MACC+DIRT(KK)%QHRS*60)).OR.          &
     &      (BACK.AND.(JET.LE.DIRT(KK)%START%MACC).AND.                 &
     &      (JET.GT.DIRT(KK)%START%MACC-DIRT(KK)%QHRS*60)))THEN

! JCL:   following lines needed to calculate mass of each emitted particle

! JCL:          since MASS cannot be < 0, take ABS of DT
!               number of time steps in emission period
                QSTEP=MAX(DBLE(1.0),60.0*DIRT(KK)%QHRS/ABS(DT))
!                QSTEP=MAX(1.0,60.0*DIRT(KK)%QHRS/DT)
                IF(KT.NE.0)THEN
!                  total emission (qrate from input becomes a multiplier)
                   QTOT=DIRT(KK)%QHRS*QAREA(II,JJ,KT,IH)*DIRT(KK)%QRATE
                ELSE
!                  set rate as specified in the control file
                   QTOT=DIRT(KK)%QHRS*DIRT(KK)%QRATE
                END IF

!               optional temperature adjustment for isoprene
!               uncomment lines for chemistry applications
!               IF(DIRT(KK)%IDENT.EQ.'ISOP')THEN
!                  CALL TMPVAL(ATEMP,CLAT,CLON,TBAR)
!                  adjustment from jacob et al (JGR,1993,14797-14813)
!                  QTOT=QTOT*EXP(0.096*(ATEMP-298.0))
!               END IF

!               emission per time step
                QVAL=QTOT/QSTEP
!               divide amount over the number of units emitted
                MASS(KK,KPM)=QVAL/NPAR

!               optional conversion to ppm (volume factor not included)
!               assume emission in kg, molecular weight in grams
!               uncomment lines for chemistry applications
!               IF(DIRT(KK)%GPMOL.GT.0.0)
!    :             MASS(KK,KPM)=22.4E+06*MASS(KK,KPM)/DIRT(KK)%GPMOL

!           start time test
            END IF

!        pollutant type loop
         END DO

!     number of particles loop
      END DO

!     within meteo grid test
      END IF

!     non-zero emissions test
      END IF

!     number of sources loop
      END DO
      END DO

!==>check for emission cycling

!     to avoid the emission of an excessive number of particles or puffs
!     (ie the whole grid each time step), gridded emission simulations should
!     be made by cycling the emissions.  This is accomplished by emitting a
!     short burst at the cycle interval (qcycle).  The cycle interval is set
!     in SETUP.CFG.  For instance if the emission duration were set to 0.1 h,
!     then the emission rate multiplier (source term in the CONTROL file)
!     should be set to 10.0 if the cycle time is 1 h.  This would then emit the
!     correct amount of mass over the 1 hour cycle period.

      DO KK=1,NUMTYP
! JCL:   add BACK as a condition
!        test for end of emission cycle
         IF(((JET+NINT(DT)).GE.(DIRT(KK)%START%MACC                     &
     &         +NINT(DIRT(KK)%QHRS*60.0))).AND.(.NOT.BACK))THEN
!           optional restart of emissions at some later time
            KCYCL=NINT(QCYCLE)*60
            DIRT(KK)%START%MACC=DIRT(KK)%START%MACC+KCYCL
! JCL:   when BACK is T; note reversal of sign
         ELSEIF(((JET+NINT(DT)).LE.(DIRT(KK)%START%MACC                 &
     &         -NINT(DIRT(KK)%QHRS*60.0))).AND.(BACK))THEN
! JCL:      optional restart of emissions at some EARLIER time
            KCYCL=NINT(QCYCLE)*60
            DIRT(KK)%START%MACC=DIRT(KK)%START%MACC-KCYCL
         END IF
      END DO

!==>final check to determine if some emissions occured
      IF(.NOT.EMIT)THEN
         WRITE(30,*)'WARNING emsgrd: No matching pollutants found'
         WRITE(30,*)'  File content: ',POLID
         WRITE(30,*)'  Control file: ',(DIRT(KK)%IDENT,KK=1,NUMTYP)
      END IF

      RETURN
      END
