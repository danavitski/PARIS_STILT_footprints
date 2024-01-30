!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFECM           PRoFile ECMwf processes ECMWF data
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE ECMWF CONVERTS A SOUNDING ON HYBRID ECMWF PRESSURE-SIGMA
!   COORDINATES (TYPE 4) TO MODEL TERRAIN FOLLOWING SIGMA LEVEL.
!   INTERPOLATES ALL METEO VARIABLES TO MODEL VERTICAL GRID
!   HYBRID COORDINATES CONSIST OF A PRESSURE OFFSET AND SIGMA VALUE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  16 Apr 1999 (RRD) - added terrain variable
!                  20 Apr 1999 (RRD) - option to use terrain height field
!
! USAGE:  CALL PRFECM(VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,
!              P0,U0,V0,T0,Z0,NZ,PSG,ESG,
!              U,V,W,T,Q,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,ZZ,RH,DEN)
!   INPUT ARGUMENT LIST:
!     VMIX  - log       vertical mixing flag
!     QFLG  - log       specific humidity indicator
!     UFLG,TFLG - log   low level wind and temp indicator
!     PFLG,SFLG - log   surface pressure and terrain indicator
!     ZSFC  - real      terrain height (m)
!     Z0    - real      roughness length (m)
!     Z0    - real      roughness length (m)
!     NZ    - int       number of input levels
!     PSG   - real      data sigma-p profile
!     U,V   - real      horizontal wind components
!     W     - real      vertical motion component (dp/dt)
!     T     - real      data temperature profile (deg K)
!     Q     - real      data specific humidity (kg/kg)
!     NL    - int       number of output sigma levels
!     ZMDL  - real      internal model top (meters)
!     ZSG   - real      internal model output sigma levels
!   OUTPUT ARGUMENT LIST:
!     P0    - real      surface pressure at data terrain (mb)
!     U0,V0 - real      low level horizontal wind components
!     T0    - real      low level temperaure (deg K)
!     PP    - real      pressure at sigma level (mb)
!     UU,VV - real      horizontal wind components
!     WW    - real      vertical motion term (sigma/time)
!     TT    - real      virtual potential temperature (pot K)
!     ZZ    - real      internal model sigma height (meters)
!     RH    - real      relative humidity fraction (0-1)
!     DEN   - real      air density (kg/m3)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: prfecm.f90,v 1.4 2006-09-29 16:05:49 skoerner Exp $
!
!$$$

      SUBROUTINE PRFECM (VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,             &
                         P0,U0,V0,T0,Z0,NZ,PSG,ESG,                      &
                         U,V,W,T,Q,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,ZZ,RH,DEN)

      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL VMIX,QFLG,UFLG,TFLG,PFLG,SFLG

!     temporary profile arrays
      REAL*8 U(NZ),V(NZ),W(NZ),T(NZ),Q(NZ),PSG(NZ),ESG(NZ)
      REAL*8 ZSG(NL),PP(NL),TT(NL),ZZ(NL),RH(NL),DEN(NL),             &
             UU(NL),VV(NL),WW(NL)

      SAVE KFLAG

!     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
      DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/

!     diagnostic
      DATA KFLAG/0/

!---------------------------------------------------------------------------------------------------
!=>remap sigma coordinates based upon surface pressure
!  zsg composed of integer offset (10's of mb) plus fraction sigma
!  use esg variable to hold corrected sigma from zsg

      DO K=1,NZ
         KPRES = INT(PSG(K))
         SIGMA = PSG(K)-KPRES
         ESG(K) = (KPRES+P0*SIGMA)/P0
      END DO

!=>compute height of ground surface if not otherwise available

      IF (.NOT.SFLG) THEN
!         computation assumes surface pressure field
!         estimate height of 1013 hPa using lowest upper level temp
          ZSFC = DLOG(1013.0/P0)*RDRY*T(1)/GRAV
      END IF

!=>compute the log of the surface pressure

      IF (PFLG) THEN
!        surface pressure from input data
         PBOT = DLOG(P0)
      ELSE
!        use lowest upper level temp to estimate surface pressure
         PBOT = DLOG(DBLE(1013.))-ZSFC*GRAV/RDRY/T(1)
         P0 = EXP(PBOT)
      END IF

!=>use adiabatic profile to estimate surface temperature

      IF (TFLG) THEN
!        available use low level value
         TBOT = T0
      ELSE
!        drop adiabatically to surface
         TBOT = T(1)*ESG(1)**(-0.286)
      END IF
!     convert all temperatures to virtual
      IF (QFLG) TBOT = TBOT*(1.0+0.61*Q(1))

!=>convert specific humidity to fraction of RH. to fraction

      IF (QFLG) THEN
         ESAT = EXP(21.4-(5351.0/T(1)))
         RBOT = Q(1)*P0/(0.622*ESAT)
      ELSE
         RBOT = Q(1)/100.0
      END IF

!=>estimate remaining surface values

!     initial vertical index
      KL = 1
!     initial output sigma level
! CHG(09/10/03) correct transformation between sigma and agl
!      ZLVL = ZMDL*(1.0-ZSG(KL))
!     adjust for terrain compression
!      ZLVL = ZMDL*ZLVL/(ZMDL-ZSFC)
      ZLVL = (1.0-ZSG(KL))*(ZMDL-ZSFC)

!     estimate surface values
      ZBOT = 0.0
      WBOT = 0.0

!=>go through each input level and interpolate to model level

      DO KZ=1,NZ

!        log of pressure at level
         PTOP = DLOG(P0*ESG(KZ))

!        virtual temperature
         IF (QFLG) THEN
            TTOP = (1.0+0.61*Q(KZ))*T(KZ)
         ELSE
            TTOP = T(KZ)
         END IF

!        use layer average for hypsometric equation
         TBAR = 0.5*(TTOP+TBOT)
         DELZ = (PBOT-PTOP)*RDRY*TBAR/GRAV
         ZTOP = ZBOT+DELZ

!        only for the first input data level
         IF (KZ == 1) THEN
!           set low level (ground) wind data if available
            IF (UFLG) THEN
               UBOT = U0
               VBOT = V0
            ELSE
!              use neutral log-law when output below data level
               IF (ZLVL < ZTOP) THEN
                  ATOP = DLOG(ZTOP/Z0)
                  ABOT = DLOG(ZLVL/Z0)
                  UBOT = U(1)*ABOT/ATOP
                  VBOT = V(1)*ABOT/ATOP
               ELSE
                  UBOT = U(1)
                  VBOT = V(1)
               END IF
            END IF
         END IF

!        convert to rh fraction
         IF (QFLG) THEN
            ESAT = EXP(21.4-(5351.0/T(KZ)))
            RTOP = Q(KZ)*ESG(KZ)*P0/(0.622*ESAT)
         ELSE
            RTOP = Q(KZ)/100.0
         END IF

         UTOP = U(KZ)
         VTOP = V(KZ)
         WTOP = W(KZ)

         DO WHILE (ZLVL <= ZTOP)
!           height <0 flags extrapolated levels not to be
!           used in stability calculations
            IF (KZ == 1 .AND. (.NOT.TFLG) .AND. VMIX) THEN
               ZZ(KL)=-ZLVL
            ELSE
               ZZ(KL) = ZLVL
            END IF

!           basic linear interpolation
            FRAC = (ZLVL-ZBOT)/DELZ
            TT(KL) = FRAC*(TTOP-TBOT)+TBOT
            RH(KL) = FRAC*(RTOP-RBOT)+RBOT
            UU(KL) = FRAC*(UTOP-UBOT)+UBOT
            VV(KL) = FRAC*(VTOP-VBOT)+VBOT

!           linear interpolation of log of pressure
            PP(KL) = EXP(FRAC*(PTOP-PBOT)+PBOT)

!           density and potential temperature from local value
            DEN(KL) = P2JM*PP(KL)/(TT(KL)*RDRY)
            TT(KL) = TT(KL)*(1000.0/PP(KL))**0.286

!           vertical velocity term converted to sigma/time
            OMEGA = FRAC*(WTOP-WBOT)+WBOT
            WW(KL) = P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))

            KL = KL+1
            IF (KL > NL) RETURN
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL = ZMDL*(1.0-ZSG(KL))
!           adjust for terrain compression
!            ZLVL = ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL = (1.0-ZSG(KL))*(ZMDL-ZSFC)
         END DO

         PBOT = PTOP
         ZBOT = ZTOP
         RBOT = RTOP
         TBOT = TTOP
         UBOT = UTOP
         VBOT = VTOP
         WBOT = WTOP

      END DO

!     sounding ends but levels remain to be filled
      IF (KL <= NL) THEN
         IF (KFLAG == 0) THEN
            WRITE(30,*)'WARNING prfecm: data extrapolation'
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL = ZMDL*(1.0-ZSG(KL))
!           adjust for terrain compression
!            ZLVL = ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL = (1.0-ZSG(KL))*(ZMDL-ZSFC)
            WRITE(30,*)'from level (m): ',ZLVL
            KFLAG = 1
         END IF

         DO KK=KL,NL
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL = ZMDL*(1.0-ZSG(KK))
!           adjust for terrain compression
!            ZLVL = ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL = (1.0-ZSG(KK))*(ZMDL-ZSFC)
            DELZ = ZLVL-ZZ(KK-1)
            ZZ(KK) = ZLVL
            TT(KK) = TT(KK-1)
            RH(KK) = RTOP

!           use previous pressure and density to find pressure
            PP(KK) = PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM

!           convert potential back to local temperature
            TEMP = TT(KK)*(PP(KK)/1000.0)**0.286
            DEN(KK) = P2JM*PP(KK)/(TEMP*RDRY)
            UU(KK) = UU(KK-1)
            VV(KK) = VV(KK-1)

!           diminish magnitude when no data
            WW(KK) = WW(KK-1)/2.0
         END DO
      END IF
      RETURN
      END
