!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSPNT           EMiSsion puff/particle at a PoiNT
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INITIALIZATION STARTS A NEW PUFF OR PARTICLE IF THE
!   CURRENT MODEL TIME IS WITHIN THE TIME LIMITS SPECIFIED FOR START
!   A POLLUTANT RELEASE.  EMISSIONS CAN BE STARTED AGAIN AT QCYCLE
!   INTERVALS.  IF MULTIPLE RELEASE LOCATIONS ARE DEFINED THEN THE
!   EMISSIONS ARE UNIFORMLY DISTRIBUTED WITHIN A LAYER FROM THE LAST
!   STARTING HEIGHT TO THE CURRENT STARTING HEIGHT WHEN TWO RELEASE
!   ARE IN THE SAME LOCATION.  OTHERWISE IT IS A POINT SOURCE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 22 Apr 1999 (RRD) - added location specific emissions
!
! USAGE:  CALL EMSPNT(NLOC,NUMTYP,KPM,INITD,DT,JET,KG,
!              NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,
!              QCYCLE,NUMPAR)
!   INPUT ARGUMENT LIST:
!     NLOC      - int   total number of source locations
!     NUMTYP    - int   number of pollutant types
!     KPM       - int   number of puffs or particles
!     INITD     - int   initial distribution type
!     DT        - real  time step (min)
!     JET       - int   current elapsed time (min)
!     KG        - int   default initial grid index
!     QCYCLE    - real  optional emission cycle time in hours
!     NUMPAR    - int   maximum number of particles permitted
!   OUTPUT ARGUMENT LIST:
!     NSORT     - int   array index of sorted elements
!     MASS      - real  array mass of pollutant (arbitrary units)
!     XPOS,YPOS - real  array horizontal position (grid units)
!     ZPOS      - real  array puff center height (sigma)
!     SIGH,SIGV - real  arrau horiz (meters) and vert sigma (sigma)
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
! $Id: emspnt.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

! JCL: add argument BACK to do different tests when model running back
      SUBROUTINE EMSPNT(NLOC,NUMTYP,KPM,INITD,DT,JET,KG,                &
     &   NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,  &
     &   QCYCLE,NUMPAR,BACK)

!      SUBROUTINE EMSPNT(NLOC,NUMTYP,KPM,INITD,DT,JET,KG,
!     :   NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,
!     :   QCYCLE,NUMPAR)

      use module_defconc
      use module_defspot

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'
!     multiple source information
!      INCLUDE 'DEFSPOT.INC'

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR),SIGX(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER   PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8   MASS(MAXDIM,MAXPAR)

!     pollutant sort index
      INTEGER   NSORT(MAXPAR)

!     flag for line source emissions
      LOGICAL EMIT

! JCL:
      LOGICAL BACK

!      COMMON /GBLSPT/ SPOT
!      COMMON /GBLCON/ CONC, DIRT

!     top-hat radius
      DATA SIGR/1.5/

!==>check to determine if emissions required

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
           EMIT=.TRUE.
         END IF
      END DO
! CHG(11/04/02) don't emit particles when not wanted!
      IF(.NOT.EMIT.OR.NUMPAR.EQ.0)RETURN
!      IF(.NOT.EMIT)RETURN

!==>determine type of emission and number of locations

!   multiple starting locations can either be point sources
!   if each starting point is in a different location or vertical
!   line sources if two starting points are in the same location.
!   The line source is distributed between the given heights

      IF(NLOC.GT.1)THEN
         SPOT(1)%ZV=SPOT(1)%ZP
         DO N=2,NLOC
! JCL:(5/5/00)don't want weird line source
!            IF(INT(SPOT(N)%XP*10000.0).EQ.
!     :         INT(SPOT(N-1)%XP*10000.0).AND.
!     :         INT(SPOT(N)%YP*10000.0).EQ.
!     :         INT(SPOT(N-1)%YP*10000.0))THEN

!              when position the same move previous point release
!              height into ZV (line source defined as ZV->ZP)
!              then only emit at locations with ZV<>0
! JCL:          comment these out too
!               SPOT(N)%ZV=SPOT(N-1)%ZP
!               SPOT(N-1)%ZV=0.0
!            ELSE
! JCL:         just have this statement
!              point source bottom equals release height
               SPOT(N)%ZV=SPOT(N)%ZP
!            END IF
         END DO

!        count up the number of different x,y release locations
!        with
         NRPTS=0
         DO N=1,NLOC
!           multiple sources can be defined through control file
!           negative mass indicates point is skipped (0 uses default)
            IF(SPOT(N)%QTRM.LT.0.0)SPOT(N)%ZV=0.0
!           those set to zero are used to determine line source heights
            IF(SPOT(N)%ZV.NE.0.0)NRPTS=NRPTS+1
         END DO
      ELSE
         NRPTS=1
         SPOT(1)%ZV=SPOT(1)%ZP
      END IF

!==>set the number of units to emit

      IF(INITD.EQ.1.OR.INITD.EQ.2)THEN
!        gaussian or top-hat emissions
         NPAR=1
      ELSE
         QSUM=0.0
!        total emissions - pollutant hours is summed because each
!        pollutant type is emitted as an independent particle
         DO KK=1,NUMTYP
! JCL:      need to take ABS of DT since DT can be < 0
            QSUM=QSUM+DMAX1(ABS(DT)/60.0,DIRT(KK)%QHRS)
!            QSUM=QSUM+AMAX1(DT/60.0,DIRT(KK)%QHRS)
         END DO
!        particle emission rate: number particles per hour per source
         NPHR=NUMPAR/NRPTS/QSUM
! JCL:   need to take ABS since DT can be < 0 when BACK is T
!        particle emissions per time step
         NPAR=MAX0(1, ABS(INT(NPHR*DT/60.0)))
!         NPAR=MAX0(1, INT(NPHR*DT/60.0))
      END IF

!==>loop through the number of independent source locations

      DO N=1,NLOC

!        check for source skip (vertical line source)
         IF(SPOT(N)%ZV.NE.0.0)THEN

!     loop through pollutants (up to MAXTYP)
!     each pollutant type will start its own trajectory unless
!     this routine is modified accordingly (MAXDIM >1)

      DO KK=1,NUMTYP

!        different emissions go into index 1 unless maxdim>1
         KT=MIN0(KK,MAXDIM)

! JCL: add BACK in condition for pollutant requiring a start
!        check if this pollutant requires a start
         IF(((.NOT.BACK).AND.(JET.GE.DIRT(KK)%START%MACC).AND.          &
     &      (JET.LT.DIRT(KK)%START%MACC+DIRT(KK)%QHRS*60)).OR.          &
     &      (BACK.AND.(JET.LE.DIRT(KK)%START%MACC).AND.                 &
     &      (JET.GT.DIRT(KK)%START%MACC-DIRT(KK)%QHRS*60)))THEN

!           multiple emissions only for particles
            DO KP=1,NPAR

            KPM=KPM+1
            IF(KPM.GT.MAXPAR)THEN
               KPM=MAXPAR
               WRITE(30,*)'Warning: emsini - exceeding puff limit'
               RETURN
            END IF
            NSORT(KPM)=KPM

!           initial position from main program
            XPOS(KPM)=SPOT(N)%XP
            YPOS(KPM)=SPOT(N)%YP
            ZPOS(KPM)=SPOT(N)%ZP

!           variances all start at zero
            SIGV(KPM)=0.0
            SIGX(KPM)=0.0

            IF(SPOT(N)%AREA.LE.0.0.OR.INITD.EQ.0)THEN
!              points source defined for all species
               SIGH(KPM)=0.0
            ELSE
!              defined by source compute sigma for uniform radius
               SIGH(KPM)=SQRT(SPOT(N)%AREA/3.14159)/SIGR
            END IF


!           multiple locations get initial vertical distribution
!           when the first two starting points at the same x,y position
            DELZ=ABS(SPOT(N)%ZV-SPOT(N)%ZP)
            IF(DELZ.GT.0.0)THEN
               ZPL=DMAX1(SPOT(N)%ZV,SPOT(N)%ZP)
               IF(INITD.EQ.1.OR.INITD.EQ.2)THEN
!                 puff variance set to layer depth
                  SIGV(KPM)=ABS(DELZ/SIGR/2.0)
                  ZPOS(KPM)=ZPL-DELZ/2.0
               ELSE
!                 particles get distributed in the layer
                  ZPOS(KPM)=ZPL-KP*DELZ/NPAR
               END IF
            END IF

!           initial distribution (see main for definitions)
            HDWP(KPM)=INITD
!           initial age at zero
            PAGE(KPM)=0
!           pollutant type always defined from pollutant index
            PTYP(KPM)=KK
!           initial grid is the default startup grid from main
            PGRD(KPM)=KG

! JCL:      Following steps are necessary to calculate the mass/particle

            IF(SPOT(N)%QTRM.EQ.0.0)THEN
!              emission defined by species
               QTOT=DIRT(KK)%QHRS*DIRT(KK)%QRATE
            ELSE
!              emission defined by source location
               QTOT=DIRT(KK)%QHRS*SPOT(N)%QTRM
            END IF

! JCL:      since mass cannot be < 0, take ABS of DT
!           number of time steps in emission period
            QSTEP=DMAX1(DBLE(1.0),60.0*DIRT(KK)%QHRS/ABS(DT))
!            QSTEP=AMAX1(1.0,60.0*DIRT(KK)%QHRS/DT)
!           emission per time step
            QVAL=QTOT/QSTEP
!           divide amount over the number of units emitted
            MASS(KT,KPM)=QVAL/NPAR

!           particle loop
            END DO

!        start time test
         END IF

!     pollutant type loop
      END DO

!     emission test
      END IF

!     number of sources loop
      END DO

!==>check for emission cycling

      DO KK=1,NUMTYP
! JCL:   add BACK as a condition
!        test for end of emission cycle
         IF(((JET+NINT(DT)).GE.(DIRT(KK)%START%MACC                     &
     &         +NINT(DIRT(KK)%QHRS*60.0))).AND.(.NOT.BACK))THEN
!           optional restart of emissions at some LATER time
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

      RETURN
      END
