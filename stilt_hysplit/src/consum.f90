!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONSUM           CONcentration SUMmation puffs or particles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION SUMMATION FOR PUFFS OR PARTICLES IS CALLED EACH TIME
!   PARTICLES ARE SUMMED TO A CELL WHILE PUFF INTERSECTIONS TO EACH
!   ARE COMPUTED AND SUMMED.  NO INTERPOLATION IS PERFORMED AND THE
!   TIME STEP SHOULD BE SUFFICIENTLY SMALL WITH RESPECT TO THE
!   MINIMUM CONCENTRATION GRID CELL SIZE. ROUTINE CALLED FOR EACH
!   PARTICLE OR PUFF AFTER ITS NEW ADVECTION POSITION IS COMPUTED.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 24 Nov 1997 (RRD)
!                 20 Apr 1999 (RRD) - added terrain compression adjustment
!                 20 Oct 1999 (RRD) - correction for dateline
!                 08 Nov 1999 (RRD) - grid scan limits test simplified
!
! USAGE:  CALL CONSUM(NUMGRD,PLAT,PLON,DT,JET,ZMDL,ZSFC,
!              CGSIZE,MASS,DEPT,ZPOS,SIGH,SIGV,HDWP,PTYP,CSUM)
!   INPUT ARGUMENT LIST:
!     NUMGRD    - int   number of concentration grids
!     PLAT      - real  particle position latitude
!     PLON      - real  particle position longitude
!     DT        - real  time step (min)
!     JET       - int   elapsed time (min)
!     ZMDL      - real  model domain top
!     ZSFC      - real  height of ground surface (m)
!     CGSIZE    - real  minimum sampling grid size (km)
!     MASS      - real  mass of pollutant (arbitrary units)
!     DEPT      - real  deposition amount of mass
!     ZPOS      - real  puff center height (sigma)
!     SIGH,SIGV - real  horiz (meters) and vert sigma (sigma)
!     HDWP      - int   Horizontal distribution within pollutant
!     PTYP      - int   pollutant type index number
!     CSUM      - real  concentration summation matrix
!   OUTPUT ARGUMENT LIST:
!     CSUM      - real  concentration summation matrix
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: consum.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

! JCL: add BACK as argument and have different tests when BACK is T
      SUBROUTINE CONSUM(NUMGRD,PLAT,PLON,DT,JET,ZMDL,ZSFC,              &
     &   CGSIZE,MASS,DEPT,ZPOS,SIGH,SIGV,HDWP,PTYP,CSUM,BACK)

!      SUBROUTINE CONSUM(NUMGRD,PLAT,PLON,DT,JET,ZMDL,ZSFC,
!     :   CGSIZE,MASS,DEPT,ZPOS,SIGH,SIGV,HDWP,PTYP,CSUM)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

!     pollutant distribution and type
      INTEGER HDWP,PTYP
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8 MASS(MAXDIM)

!     master concentration array (x,y,z,grids,species)
      REAL*8 CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!     deposition array
      REAL*8 DEPT(MAXTYP)

! JCL:
      LOGICAL BACK

!      COMMON /GBLCON/ CONC, DIRT

!     distribution     meters/deg        deg/radian
      DATA SIGR/1.54/, YMPD/111198.323/, DGPR/57.295828/

!     Gaussian constants (2 PI)
      DATA TWOPI / 6.283185 /

!     meters per degree longitude
      XMPD=YMPD*DCOS(PLAT/DGPR)

!     set default pollutant type (over-ride if MAXDIM>1)
      KT=PTYP

!     convert particle position to meters
! CHG(09/10/03) adapt change from version 46
!      ZPAR=ZMDL*(1.0-ZPOS)
!     terrain compression adjustment
!      ZPAR=ZMDL*ZPAR/(ZMDL-ZSFC)
      ZPAR=(ZMDL-ZSFC)*(1.0-ZPOS)

!     loop through the number of concentration grids
      DO KG=1,NUMGRD

! JCL:  change conditions if BACK is TRUE
!     keep track of active minimum for delt-t calculation
      IF (((JET.LT.CONC(KG)%STOP%MACC).AND.(.NOT.BACK)).OR. &
          ((JET.GT.CONC(KG)%STOP%MACC).AND.BACK))           &
         CGSIZE=DMIN1(CONC(KG)%SIZE,CGSIZE)
! JCL: when BACK is F
      IF (.NOT.BACK)THEN
        IF((JET.LT.CONC(KG)%START%MACC).OR.                             &
     &     (JET.GE.CONC(KG)%STOP%MACC))THEN
! JCL:
           PRINT *, 'CONSUM PROBLEM!', JET.LT.CONC(KG)%START%MACC
           GO TO 100
        END IF
! JCL: when BACK is T
      ELSEIF((JET.GT.CONC(KG)%START%MACC).OR.                           &
     &     (JET.LE.CONC(KG)%STOP%MACC))THEN
! JCL:
           PRINT *, 'CONSUM PROBLEM!2', JET.LT.CONC(KG)%STOP%MACC,      &
     &                                   JET.GT.CONC(KG)%START%MACC
           GO TO 100
! JCL: end tests with BACK
      END IF

!     compute radius (m) of horizontal plume extent
      IF(HDWP.EQ.0)THEN
!           particle model - no distribution
            RADIUS=1.0E+25
!           normalized dispersion factor = cell area
            DFXY0=1.0/(XMPD*CONC(KG)%DELT_LON*YMPD*CONC(KG)%DELT_LAT)
!           vertical component depends upon cell size (see below)

      ELSEIF(HDWP.EQ.1.OR.HDWP.EQ.3)THEN
!           gaussian model - scan to 3 sigma
            RADIUS=3.0*SIGH
!           plume vertical extent
! CHG(09/10/03) adapt change from version 46
!            ZBOT=ZMDL*(1.0-DMAX1(ZPOS+SIGR*SIGV,DBLE(0.0)))
!            ZTOP=ZMDL*(1.0-DMIN1(ZPOS-SIGR*SIGV,DBLE(1.0)))
!           terrain compression adjustment
!            ZBOT=ZMDL*ZBOT/(ZMDL-ZSFC)
!            ZTOP=ZMDL*ZTOP/(ZMDL-ZSFC)
            ZBOT=(ZMDL-ZSFC)*(1.0-DMAX1(ZPOS+SIGR*SIGV,DBLE(0.0)))
            ZTOP=(ZMDL-ZSFC)*(1.0-DMIN1(ZPOS-SIGR*SIGV,DBLE(1.0)))

!           vertical distribution uniform over layer
            VDIST=DMAX1(ZTOP-ZBOT,DBLE(1.0))
!           horizontal dispersion factor for area deposition
            DFXY0=1.0/(TWOPI*SIGH*SIGH)
!           uniform concentration within volume
            DFXYZ=DFXY0/VDIST

      ELSEIF(HDWP.EQ.2.OR.HDWP.EQ.4)THEN
!           top hat distribution - scan 1.54 sigma
            RADIUS=SIGR*SIGH
!           plume vertical extent
            ZBOT=ZMDL*(1.0-DMAX1(ZPOS+SIGR*SIGV,DBLE(0.0)))
            ZTOP=ZMDL*(1.0-DMIN1(ZPOS-SIGR*SIGV,DBLE(1.0)))
!           terrain compression adjustment
            ZBOT=ZMDL*ZBOT/(ZMDL-ZSFC)
            ZTOP=ZMDL*ZTOP/(ZMDL-ZSFC)
!           vertical distribution uniform over layer
            VDIST=DMAX1(ZTOP-ZBOT,DBLE(1.0))
!           horizontal dispersion factor for area deposition
            DFXY0=1.0/(3.14159*RADIUS*RADIUS)
!           uniform concentration within volume
            DFXYZ=DFXY0/VDIST
      END IF

      IF(CONC(KG)%SNAP.EQ.0)THEN
!        factors for integrations per sampling interval
! JCL:   take ABS of DT
         SFACT=ABS(DT)/CONC(KG)%DELTA%MACC
!         SFACT=DT/CONC(KG)%DELTA%MACC
      ELSE
!        no factor required for snapshot maps
         SFACT=1.0
      END IF

!     longitude correction to avoid dateline problems (RRD - 10/20/99)
      IF(CONC(KG)%X1Y1_LON.GE.0.0.AND.PLON.LT.0.0)THEN
         GLON=PLON+360.0
      ELSE
         GLON=PLON
      END IF

!     find the grid position on current concentration grid
      XP=1.0+(GLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON
      YP=1.0+(PLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT

!     for particles (no distribution) insure no scan
      IF(HDWP.EQ.0)THEN
!        round position so that cell centered over point
         XI1=NINT(XP)
         YJ1=NINT(YP)
         XI2=XI1
         YJ2=YJ1

      ELSE
!        find lower left corner based upon scan radius
         CLAT=PLAT-RADIUS/YMPD
         CLON=GLON-RADIUS/XMPD

!        convert corner to grid index number
         XI1=1.0+(CLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON
         YJ1=1.0+(CLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT

!        opposite corner from puff position for scan
         XI2=2.0*XP-XI1
         YJ2=2.0*YP-YJ1
      END IF

!     convert real grid to nearest integer value
      I1=INT(XI1)
      I2=NINT(XI2)
      J1=INT(YJ1)
      J2=NINT(YJ2)

!     loop through grid (sampling) points within plume
      DO II=I1,I2
      DO JJ=J1,J2

!        test limits
         IF(II.LT.1.OR.II.GT.CONC(KG)%NUMB_LON)GOTO 50
         IF(JJ.LT.1.OR.JJ.GT.CONC(KG)%NUMB_LAT)GOTO 50

!        distance (meters) of sampling point to plume center
         CLON=(II-1.0)*CONC(KG)%DELT_LON+CONC(KG)%X1Y1_LON
         CLAT=(JJ-1.0)*CONC(KG)%DELT_LAT+CONC(KG)%X1Y1_LAT
         DELY=(PLAT-CLAT)*YMPD
         DELX=(GLON-CLON)*XMPD
         DIST=DSQRT(DELY*DELY+DELX*DELX)

!        grid point within plume radius
         IF(DIST.LT.RADIUS)THEN

!           sum through all levels within plume
            NLVL=CONC(KG)%LEVELS

!           vertical particles with various horizontal distributions
            IF(HDWP.EQ.0.OR.HDWP.EQ.3.OR.HDWP.EQ.4)THEN

               GFACT=1.0
!              horizontal gaussian - vertical particle
               IF(HDWP.EQ.3)                                            &
     &            GFACT=DEXP(-0.5*DIST*DIST/SIGH/SIGH)

!              all variations use particles in vertical
               ZBOT=0
               DO KL=1,NLVL
!                 determine vertical cell sizes (input defines top)
                  ZTOP=CONC(KG)%HEIGHT(KL)
                  VDIST=DMAX1(ZTOP-ZBOT,DBLE(1.0))

!                 check for simultaneous species at position
                  DO KK=1,MAXDIM
!                    use default if multiple species not defined
                     IF(MAXDIM.GT.1)KT=KK

!                    deposition output requested
                     IF(ZBOT.EQ.0.AND.ZTOP.EQ.0.0)THEN
                        CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)+      &
     &                     DEPT(KK)*DFXY0*GFACT

!                    particle falls within vertical extent of cell
                     ELSEIF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)THEN
!                       air concentration option
                        CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)+      &
     &                     MASS(KK)*SFACT*DFXY0*GFACT/VDIST

                     END IF
                  END DO

                  ZBOT=ZTOP
!              level loop
               END DO

!           gaussian or top-hat distributions
            ELSE

               GFACT=1.0
!              compute horizontal gaussian factor if required
               IF(HDWP.EQ.1)                                            &
     &            GFACT=DEXP(-0.5*DIST*DIST/SIGH/SIGH)

               DO KL=1,NLVL
!                 vertical concentration factors for Gaussian
                  ZLVL=CONC(KG)%HEIGHT(KL)

                  DO KK=1,MAXDIM
                     IF(MAXDIM.GT.1)KT=KK

!                    deposition output requested (level=0 defined)
                     IF(ZLVL.EQ.0.0)THEN
                        CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)+      &
     &                     DEPT(KK)*DFXY0*GFACT

!                    output grid node within vertical plume
                     ELSEIF(ZLVL.GE.ZBOT.AND.ZLVL.LT.ZTOP)THEN
                        CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)+      &
     &                     MASS(KK)*SFACT*DFXYZ*GFACT
                     END IF

                  END DO
               END DO

!           distribution test
            END IF

!        within radius test
         END IF

!        plume within grid test
   50    CONTINUE

!     i,j loop
      END DO
      END DO

  100 CONTINUE

!     grid loop
      END DO

      RETURN
      END
