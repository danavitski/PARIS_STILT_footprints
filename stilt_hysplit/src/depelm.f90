!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPELM           DEPosition of a pollutant ELeMent
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION OF A POLLUTANT ELEMENT COMPUTES GRAVITATIONAL SETTLING,
!   DRY DEPOSITION EITHER EXPLICIT OR VIA THE RESISTANCE METHOD,
!   WET REMOVAL, AND RADIOACTIVE DECAY AS APPLIED TO ONE POLLUTANT
!   PARTICLE OR PUFF EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 27 Oct 1998 (MDC) - Below cloud scavenging correction
!                                   - Added explicit gravitational settling
!                 20 Apr 1999 (RRD) - terrain compression factor
!                 18 Jun 1999 (RRD) - consistent gas definition
!
! USAGE:  CALL DEPELM(NLVL,ZX,DT,ZMDL,ZSFC,ZSG,MASS,DEPT,ZPOS,SIGV,KTYP,
!              LAND,ROUG,SFCL,USTR,PSI,SFLX,HDWP,RAIN,DD,TT,QQ)
!   INPUT ARGUMENT LIST:
!     NLVL  - int number of vertical levels
!     ZX        - real  vertical array position (index)
!     DT        - real  time step (min)
!     ZMDL      - real  model domain top (meters)
!     ZSFC  - real      height of terrain surface (m)
!     ZSG   - real      model sigma levels
!     MASS      - real  mass of pollutant (arbitrary units)
!     ZPOS      - real  puff center height (sigma)
!     SIGV  - real      vert sigma (sigma)
!     KTYP      - int   pollutant type index number
!     LAND  - int land use category (1-11)
!     ROUG  - real      aerodynamic roughness length (m)
!     SFCL  - real      height of the surface layer (m)
!     USTR  - real      friction velocity (m/s)
!     PSI   - real      integrated stability function for heat
!     SFLX  - real      incident short wave flux (w/m2)
!     HDWP      - int   pollutant distribution type (index)
!     RAIN  - real      precipitation value (m/min)
!     DD        - real  air density profile (kg/m3)
!     TT    - real      temperature profile
!     QQ    - real      humidity profile (fraction 0-1)
!   OUTPUT ARGUMENT LIST:
!     DEPT      - real  deposition total (mass units)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: depelm.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

      SUBROUTINE DEPELM(NLVL,ZX,DT,ZMDL,ZSFC,ZSG,MASS,DEPT,ZPOS,SIGV,   &
     &   KTYP,LAND,ROUG,SFCL,USTR,PSI,SFLX,HDWP,RAIN,DD,TT,QQ)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     concentration and pollutant structure
!      INCLUDE 'DEFCONC.INC'

      INTEGER   HDWP
!     pollutant mass array, profile variables and sigma levels
      REAL*8   MASS(MAXDIM),  DD(NLVL),TT(NLVL),QQ(NLVL),ZSG(NLVL)
!     return deposition array
      REAL*8   DEPT(MAXTYP)

!      COMMON /GBLCON/ CONC, DIRT

!          GRAVITY      DYNAMIC VISCOSITY  MEAN FREE PATH
!          (m/s2)       (g m-1 s-1)        (m at stp)
      data grav/9.801/, dmvc/1.789E-02/,   frep/6.53E-08/

!          STP DENSITY    gas constant
!          (g/m3)         (atm-liter / deg-mole)
      data dstp/1.2E+03/, rgas/0.082057/

!     vertical puff distribution scan factor
      data SIGR/1.54/

!     rounded vertical index position for meteorology profiles
      KLVL=NINT(ZX)

      IF(HDWP.EQ.1.OR.HDWP.EQ.2)THEN
!        height value at bottom and top of puff in meters
         PBOT=(ZMDL-ZSFC)*(1.0-DMIN1(ZPOS+SIGR*SIGV,DBLE(1.0)) )
         PTOP=(ZMDL-ZSFC)*(1.0-DMAX1(ZPOS-SIGR*SIGV,DBLE(0.0)) )
         PDEPTH=DMAX1(PTOP-PBOT,DBLE(1.0))
      ELSE
!        particle distribution for vertical
         PBOT=(ZMDL-ZSFC)*(1.0-ZPOS)
         PTOP=PBOT
         KK=MIN0(INT(ZX)+1,NLVL)
!        assume depth equal to meteorological cell size
         PDEPTH=(ZMDL-ZSFC)*(ZSG(KK-1)-ZSG(KK))
         PDEPTH=MAX(SFCL,PDEPTH)
      END IF

!     set default pollutant type (over-ride if MAXDIM>1)
      KT=KTYP

!     check for simultaneous species at position
      DO KK=1,MAXDIM

!        use default if multiple species not defined
         IF(MAXDIM.GT.1)KT=KK

         IF(DIRT(KT)%DODRY)THEN
!           explicit definition of the dry deposition velocity
            VD=DIRT(KT)%DRYVL
!           set gravitational settling if defined as particle
            IF(DIRT(KT)%DOGAS)THEN
               VG=0.0
            ELSE
               VG=VD
            END IF

         ELSEIF(DIRT(KT)%DOGRV.OR.DIRT(KT)%DORES)THEN
!           local air density
            AIRD=DD(KLVL)
!           convert (kg/m3) --> (g/m3)
            AIRD=AIRD*1000.0

!           compute gravitational settling for particles
            IF(DIRT(KT)%DOGRV)THEN
!              particle density from g/cc g/m3
               DENS=DIRT(KT)%PDENS*1.0E+06
!              particle diameter (um) to (m)
               PDIA=DIRT(KT)%PDIAM*1.0E-06
!              base settling velocity
               VB=PDIA*PDIA*GRAV*(DENS-AIRD)/(18.0*DMVC)
!              slip correction (mean free path = particle diameter)
               FREA=FREP*(DSTP/AIRD)
               SC=1.0+(2.0*FREA/PDIA)*(1.26+0.4*DEXP(-0.55*PDIA/FREA))
!              final value apply shape correction factor
               VG=VB*SC/DIRT(KT)%SHAPE
            ELSE
               VG=0.0
            END IF

            IF(DIRT(KT)%DORES)THEN
!              compute resistance based VD
               CALL DEPDRY(KT,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,TT(1),  &
     &            PDIA,VG,VD)
            ELSE
!              without resistance computation assume Vd = settling
               VD=VG
            END IF

         ELSE
            VD=0.0
            VG=0.0
         END IF

!        change in vertical position due to settling
         IF(VG.GT.0.0)THEN
            DROP=60.0*VG*DT
            PTOP=DMAX1(DBLE(0.0),PTOP-DROP)
            PBOT=DMAX1(DBLE(0.0),PBOT-DROP)
            ZPOS=1.0-0.5*(PTOP+PBOT)/(ZMDL-ZSFC)
!           adjust values for puffs - don't overwrite for particles
            IF(HDWP.EQ.1.OR.HDWP.EQ.2)THEN
               SIGV=(PTOP-PBOT)/SIGR/2.0/(ZMDL-ZSFC)
               PDEPTH=DMAX1(PTOP-PBOT,DBLE(1.0))
            END IF
         END IF

!        zero out removal rate constant (1/min)
         BETA=0.0

!        if puff within first layer compute dry removal time constant
         IF(VD.GT.0.0.AND.PBOT.LT.SFCL)BETA=60.0*VD/PDEPTH

!        test for wet removal processes
         IF(DIRT(KT)%DOWET.AND.RAIN.GT.0.0)THEN

!           determine bottom and top of the precip layer (80. to 60.)
            KBOT=0
            KTOP=NLVL
            DO K=1,NLVL
               KRH=QQ(K)*100.0+0.5
               IF(KBOT.EQ.0.AND.KRH.GE.80)KBOT=K
               IF(KBOT.NE.0.AND.KTOP.EQ.NLVL.AND.KRH.LE.60)KTOP=K
            END DO

!           default layer (to be consistent with rain even if no cloud)
            IF(KBOT.EQ.0)THEN
               DO K=1,NLVL
                  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(K))
                  IF(ZLVL.LT. 500.0)KBOT=K
                  IF(ZLVL.LT.3000.0)KTOP=K
               END DO
            END IF

!           check if any part of pollutant within cloud layer
            CTOP=(ZMDL-ZSFC)*(1.0-ZSG(KTOP))
            CBOT=(ZMDL-ZSFC)*(1.0-ZSG(KBOT))

            IF(PBOT.LT.CTOP)THEN
!              fraction of pollutant below cloud top
               FRBCT=1.0
               IF(PTOP.GT.CTOP)FRBCT=1.0-(PTOP-CTOP)/PDEPTH

!              fraction of pollutant above cloud bottom
               FRACB=1.0
               IF(PBOT.LT.CBOT)FRACB=1.0-(CBOT-PBOT)/PDEPTH
               IF(PTOP.LT.CBOT)FRACB=0.0

               IF(DIRT(KT)%DOGAS)THEN
!                 equilibrium concentration (mass units) for gases
!                 applies as long as material is below cloud top
                  IF(PBOT.LT.CTOP)THEN
!                    deposition velocity
                     DEPV=DIRT(KT)%WETGAS*RGAS*TT(KLVL)*RAIN
!                    rate constant
                     BETA=BETA+DEPV*FRBCT/PDEPTH
                  END IF

               ELSE
!                 for particles below cloud use scavenging coefficient
!                 only fraction of mass below cloud removed
!                 correction 10/27/98
                  IF(PBOT.LT.CBOT)                                      &
     &               BETA=BETA+DIRT(KT)%WETLO*(1.0-FRACB)*60.0

!                 for particles within cloud use scavenging ratio
                  IF(PTOP.GT.CBOT)THEN
!                    deposition velocity
                     DEPV=DIRT(KT)%WETIN*RAIN
!                    rate constant
                     BETA=BETA+DEPV*FRBCT*FRACB/PDEPTH
                  END IF

!              gas or particle
               END IF

!           layer test
            END IF

!        rain test
         END IF

!        test for radioactive decay
         IF(DIRT(KT)%DORAD)THEN
!           convert half-life (days) to time constant (1/min)
            RTC=ALOG(0.5)/(DIRT(KT)%RHALF*1440.0)
!           apply immediately to mass since it doesn't deposit
            MASS(KK)=MASS(KK)*DEXP(DT*RTC)
         END IF

         IF(BETA.EQ.0.0)THEN
!           no deposition
            DEPT(KK)=0.0
         ELSEIF(DT*BETA.LT.0.01)THEN
!           small removal values assume linear approximation
            DEPT(KK)=MASS(KK)*DT*BETA
            MASS(KK)=MASS(KK)-DEPT(KK)
         ELSE
!           apply exponential for removal > 1.
            DEPT(KK)=MASS(KK)*(1.0-DEXP(-DT*BETA))
!           can't remove more than exists
            DEPT(KK)=DMIN1(MASS(KK),DEPT(KK))
            MASS(KK)=MASS(KK)-DEPT(KK)
         END IF

!     species loop
      END DO

      RETURN
      END
