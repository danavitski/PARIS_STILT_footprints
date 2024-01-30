!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFPRS           processes meteo PRoFile on PReSsure surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE PRESSURE CONVERTS INPUT DATA ON AN ABSOLUTE PRESSURE
!   COORDINATE SYSTEM TO INTERNAL MODEL TERRAIN FOLLOWING COORDINATE
!   INPUT DATA HEIGHTS ASSUMED TO BE RELATIVE TO MSL, INPUT TEMPERATURE
!   ARE ASSUME TO BE VIRTUAL.  ADDITIONAL DIAGNOSTIC LOCAL VARIABLES
!   COMPUTED ARE DENSITY AND POTENTIAL TEMPERATURE.  VERTICAL DATA
!   LINEARLY INTERPOLATED.  REQUIRED VALUES OUTSIDE OF THE INPUT DATA
!   HEIGHT RANGE ARE EXTRAPOLATED AS APPROPRIATE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  20 Apr 1999 (RRD) - option to use terrain height field
!
! USAGE:  CALL PRFPRS(VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,
!              P0,U0,V0,T0,Z0,NZ,PSG,
!              Z,U,V,W,T,R,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,ZZ,RH,DEN)
!   INPUT ARGUMENT LIST:
!     VMIX  - log vertical mixing flag
!     QFLG  - log specific humidity indicator
!     UFLG,TFLG - log   low level data indicator
!     PFLG,SFLG - log   surface pressure and terrain indicator
!     ZSFC  - real      terrain height (m)
!     Z0    - real      roughness length (m)
!     NZ    - int number of input levels
!     PSG   - real      data pressure levels (mb)
!     Z         - real  absolute heights of surfaces (m)
!     U,V   - real      horizontal wind components
!     W         - real  vertical motion component (dp/dt)
!     T           - real      data temperature profile (deg K)
!     R           - real      relative humidity
!     NL        - int   number of output sigma levels
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
! $Id: prfprs.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

      SUBROUTINE PRFPRS(VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,             &
     &   P0,U0,V0,T0,Z0,NZ,PSG,                                         &
     &   Z,U,V,W,T,R,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,ZZ,RH,DEN)

      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL VMIX,QFLG,UFLG,TFLG,PFLG,SFLG

!     temporary arrays for processing sounding
      REAL*8   U(NZ),V(NZ),W(NZ),T(NZ),Z(NZ),R(NZ),PSG(NZ)
      REAL*8   ZSG(NL),PP(NL),TT(NL),ZZ(NL),RH(NL),DEN(NL),             &
     &       UU(NL),VV(NL),WW(NL)

      SAVE KFLAG

!     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
      DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/

!     extrapolation error message flag to avoid multiple messages
      DATA KFLAG/0/


!=>compute height of ground surface if not otherwise available

      IF(.NOT.SFLG)THEN
         IF(P0.GT.PSG(1))THEN
!           build down to surface using lowest upper level temp
            ZSFC=Z(1) - DLOG( P0/PSG(1) ) *RDRY*T(1)/GRAV
         ELSE
!           just interpolate between data levels
            K=1
            DO WHILE (PSG(K).GE.P0)
               K=K+1
            END DO
            FRAC=(PSG(K-1)-P0)/(PSG(K-1)-PSG(K))
            ZSFC=FRAC*(Z(K)-Z(K-1))+Z(K-1)
         END IF
      END IF

!=>compute the log of the surface pressure

      IF(PFLG)THEN
!        surface pressure from input data
         PBOT=DLOG(P0)
      ELSE
!        use lowest upper level temp to estimate surface pressure
         PBOT=(Z(1)-ZSFC)*GRAV/RDRY/T(1)+DLOG(PSG(1))
         P0=EXP(PBOT)
      END IF

!=>use adiabatic profile to estimate surface temperature

      IF(TFLG)THEN
         TBOT=T0
      ELSE
         TBOT=T(1)*(PSG(1)/P0)**(-0.286)
      END IF
!     convert temperature to virtual
      IF(QFLG)TBOT=TBOT*(1.0+0.61*R(1))

!=>convert specific humidity to fraction of RH. to fraction

      IF(QFLG)THEN
         ESAT=EXP(21.4-(5351.0/T(1)))
         RBOT=R(1)*P0/(0.622*ESAT)
      ELSE
         RBOT=R(1)/100.0
      END IF

!=>estimate remaining surface values

      ZBOT=0.0
      WBOT=0.0

!     initial vertical index of internal grid
      KL=1
!     first output height level on internal grid system
! CHG(09/10/03) correct transformation between sigma and agl
!      ZLVL=ZMDL*(1.0-ZSG(KL))
!     terrain compression adjustment
!      ZLVL=ZLVL*ZMDL/(ZMDL-ZSFC)
      ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)

!=>go through each input level and interpolate to model level

      DO KZ=1,NZ

!        set low level wind data if available
         IF(KZ.EQ.1)THEN
! jcl:      EDAS has 10-m winds, so simply use these values
            IF(UFLG)THEN
               UBOT=U0
               VBOT=V0
            ELSE
!              use neutral log-law when output below data level
               IF(ZLVL.LT.(Z(1)-ZSFC))THEN
                  ATOP=DLOG((Z(1)-ZSFC)/Z0)
                  ABOT=DLOG(ZLVL/Z0)
                  UBOT=U(1)*ABOT/ATOP
                  VBOT=V(1)*ABOT/ATOP
               ELSE
                  UBOT=U(1)
                  VBOT=V(1)
               END IF
            END IF
         END IF

!        log of pressure at level
         PTOP=DLOG(PSG(KZ))
         ZTOP=Z(KZ)-ZSFC
         UTOP=U(KZ)
         VTOP=V(KZ)
         WTOP=W(KZ)

!        virtual temperature
         IF(QFLG)THEN
            TTOP=(1.0+0.61*R(KZ))*T(KZ)
         ELSE
            TTOP=T(KZ)
         END IF

!        convert to rh fraction
         IF(QFLG)THEN
            ESAT=EXP(21.4-(5351.0/T(KZ)))
            RTOP=R(KZ)*PSG(KZ)/(0.622*ESAT)
         ELSE
            RTOP=R(KZ)/100.0
         END IF

         DO WHILE (ZLVL.LE.ZTOP)
!           height <0 flags extrapolated levels not to be
!           used in stability calculations
            IF(KZ.EQ.1.AND.(.NOT.TFLG).AND.VMIX)THEN
               ZZ(KL)=-ZLVL
            ELSE
               ZZ(KL)=ZLVL
            END IF

!           basic linear interpolation
            FRAC=(ZLVL-ZBOT)/(ZTOP-ZBOT)
            TT(KL)=FRAC*(TTOP-TBOT)+TBOT
            RH(KL)=FRAC*(RTOP-RBOT)+RBOT
            UU(KL)=FRAC*(UTOP-UBOT)+UBOT
            VV(KL)=FRAC*(VTOP-VBOT)+VBOT

!           linear interpolation of log of pressure
            PP(KL)=EXP(FRAC*(PTOP-PBOT)+PBOT)

!           density and potential temperature from local value
            DEN(KL)=P2JM*PP(KL)/(TT(KL)*RDRY)
            TT(KL)=TT(KL)*(1000.0/PP(KL))**0.286

!           vertical velocity term converted to sigma/time
            OMEGA=FRAC*(WTOP-WBOT)+WBOT
            WW(KL)=P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))

            KL=KL+1
            IF(KL.GT.NL)RETURN
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KL))
!           adjust for terrain compression
!            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)

         END DO

!        update bottom definition when above surface
         IF(KL.GT.1)THEN
            ZBOT=ZTOP
            PBOT=PTOP
            TBOT=TTOP
            UBOT=UTOP
            VBOT=VTOP
            WBOT=WTOP
            RBOT=RTOP
         END IF

!     input data level loop
      END DO

!     sounding ends but levels remain to be filled
      IF(KL.LE.NL)THEN
         IF(KFLAG.EQ.0)THEN
            WRITE(30,*)'WARNING prfprs: data extrapolation'
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KL))
!           terrain compression
!            ZLVL=ZLVL*ZMDL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)
            WRITE(30,*)'from level (m): ',ZLVL
            KFLAG=1
         END IF

         DO KK=KL,NL
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KK))
!           terrain compression
!            ZLVL=ZLVL*ZMDL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KK))*(ZMDL-ZSFC)
            DELZ=ZLVL-ZZ(KK-1)
            ZZ(KK)=ZLVL
            TT(KK)=TT(KK-1)
            RH(KK)=RH(KK-1)

!           use previous pressure and density to find pressure
            PP(KK)=PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM

!           convert potential back to local temperature
            TEMP=TT(KK)*(PP(KK)/1000.0)**0.286
            DEN(KK)=P2JM*PP(KK)/(TEMP*RDRY)
            UU(KK)=UU(KK-1)
            VV(KK)=VV(KK-1)

!           diminish magnitude when no data
            WW(KK)=WW(KK-1)/2.0
         END DO
      END IF

      RETURN
      END
