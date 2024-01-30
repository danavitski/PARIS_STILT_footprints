!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFTER           processes PRoFile on TERrain surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE TERRAIN CONVERTS METEO INPUT DATA ON TERRAIN FOLLOWING
!   (TYPE 3) TO MODEL TERRAIN FOLLOWING SIGMA LEVELS
!   SEE DOCBLOCK OF PRFPRS FOR A MORE COMPLETE DESCRIPTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 20 Jan 1999 (RRD) - coamps vertical velocity in m/s
!                  20 Apr 1999 (RRD) - added terrain variable
!
! USAGE:  CALL PRFTER(VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,
!              P0,U0,V0,T0,Z0,NZ,PSG,
!              P,U,V,W,T,Q,NL,ZMDL,ZSG,
!              PP,UU,VV,WW,TT,ZZ,RH,DEN)
!   INPUT ARGUMENT LIST:
!     VMIX  - log vertical mixing flag
!     QFLG  - log specific humidity indicator
!     UFLG,TFLG - log   low level wind and temp indicator
!     PFLG,SFLG - log   surface pressure and terrain indicator
!     ZSFC  - real      terrain height (m)
!     Z0    - real      roughness length (m)
!     NZ    - int number of input levels
!     PSG   - real      data height profile
!     P         - real  pressure at data level (mb)
!     U,V   - real      horizontal wind components
!     W         - real  vertical motion component (dp/dt)
!     T           - real      data temperature profile (deg K)
!     Q           - real      data specific humidity (kg/kg)
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
!     ZZ    - real      internal model sigma height agl (meters)
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
! $Id: prfter.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

! JCL(03/27/03): pass on flag:  whether data from RAMS or not
      SUBROUTINE PRFTER(VMIX,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,             &
     &   P0,U0,V0,T0,Z0,NZ,PSG,                                         &
     &   P,U,V,W,T,Q,NL,ZMDL,ZSG,                                       &
     &   PP,UU,VV,WW,TT,ZZ,RH,DEN,RAMSFLG)

      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL VMIX,QFLG,UFLG,TFLG,PFLG,SFLG

! JCL(03/27/03): set flag: whether data from RAMS or not
      LOGICAL RAMSFLG

!     temporary profile arrays
      REAL*8   U(NZ),V(NZ),W(NZ),T(NZ),Q(NZ),P(NZ),PSG(NZ)
      REAL*8   ZSG(NL),PP(NL),TT(NL),ZZ(NL),RH(NL),DEN(NL),             &
     &       UU(NL),VV(NL),WW(NL)

      SAVE KFLAG

!     model top of input data used for scaling terrain surfaces, the value
!     is obtained from meteorological model -  RAMS: 20000  COAMPS: 34800
!     Mesoscale models define a Z* coordinate system where the height
!     above ground of the models levels varies according to terrain height
!     such that  Z* = Ztop ( Zmsl - Zsfc ) / ( Ztop - Zsfc ).  Here we need
!     to convert back from Z* given with the input data to the actual
!     Zagl at each grid point before interpolation of data to internal
!     hysplit grid.  Above equation is then: Zagl = Z* ( 1 - Zsfc / Ztop )

!     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
      DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/

!     diagnostic
      DATA KFLAG/0/

!      WRITE(45,*)'in PRFTER:',RAMSFLG

! CHG&JCL(10/07/03): synchronize constants with values used in RAMS (see rconstants.h)
      IF(RAMSFLG)THEN
         RDRY=287.0
         GRAV=9.80
      END IF

!=>set model top of input data according to maximum data level
!  assume <20km = RAMS   >20km = COAMPS
      IF(PSG(NZ).LE.20000.0)THEN
         ZMDLT=20000.0
      ELSE
         ZMDLT=34800.0
      END IF

! CHG(09/10/03) use RAMS USP model top, now in ZMDL
      IF(RAMSFLG)ZMDLT=ZMDL

!=>compute height of ground surface if not otherwise available

      IF(.NOT.SFLG)THEN
!         computation assumes surface pressure field
!         estimate height of 1013 hPa using lowest upper level temp
          ZSFC=DLOG(1013.0/P0)*RDRY*T(1)/GRAV
      END IF

!=>compute the surface pressure if needed

      IF(.NOT.PFLG)THEN
!        use lowest upper level temp to estimate surface pressure
         P0=EXP(DLOG(DBLE(1013.))-ZSFC*GRAV/RDRY/T(1))
      END IF

!=>use adiabatic profile to estimate surface temperature

      IF(TFLG)THEN
!        available use low level value
         TBOT=T0
      ELSE
!        first compute lowest data level in agl units
         ZDAT=PSG(1)*(1.0-ZSFC/ZMDLT)
!        convert data level to sigma units
         SIGZ=1.0-ZDAT/(ZMDLT-ZSFC)
!        drop adiabatically to surface (ratio in sigma units)
         TBOT=T(1)*SIGZ**(-0.286)
      END IF
!     convert all temperatures to virtual
! CHG&JCL(10/07/03) rams: have T in K and moisture as sphu (want: T in K as non-virtual?)
      IF(QFLG.AND..NOT.RAMSFLG)TBOT=TBOT*(1.0+0.61*Q(1))
      IF(QFLG.AND.RAMSFLG)TBOT=TBOT*(1.0+0.61*Q(1))
!      IF(QFLG.AND.RAMSFLG)TBOT=TBOT

!=>convert specific humidity to fraction of RH. to fraction

      IF(QFLG)THEN
         ESAT=EXP(21.4-(5351.0/T(1)))
         RBOT=Q(1)*P0/(0.622*ESAT)
      ELSE
         RBOT=Q(1)/100.0
      END IF

!=>interpolate the input data to the internal model grid

!     initial vertical index for output
      KL=1
!     initial internal sigma level for output
! CHG(09/10/03) correct transformation between sigma and agl
!      ZLVL=ZMDL*(1.0-ZSG(KL))
!     adjust for terrain compression
!      ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
      ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)

!     set surface values
      ZBOT=0.0
      PBOT=P0
      WBOT=0.0

      DO KZ=1,NZ

!        set low level wind data if available
         IF(KZ.EQ.1)THEN
            IF(UFLG)THEN
               UBOT=U0
               VBOT=V0
            ELSE
!              use neutral log-law when output below data level
!              convert normalized units to agl for level
               ZTOP=PSG(1)*(1.0-ZSFC/ZMDLT)
               IF(ZLVL.LT.ZTOP)THEN
                  ATOP=DLOG(ZTOP/Z0)
                  ABOT=DLOG(ZLVL/Z0)
                  UBOT=U(1)*ABOT/ATOP
                  VBOT=V(1)*ABOT/ATOP
               ELSE
                  UBOT=U(1)
                  VBOT=V(1)
               END IF
            END IF
         END IF

!        virtual temperature
! CHG&JCL(10/07/03) rams: have T in K and moisture as sphu (want: T in K as non-virtual?)
         IF(QFLG.AND..NOT.RAMSFLG)THEN
            TTOP=(1.0+0.61*Q(KZ))*T(KZ)
         ELSEIF(QFLG.AND.RAMSFLG)THEN
            TTOP=(1.0+0.61*Q(KZ))*T(KZ)
!            TTOP=T(KZ)
         ELSE
            TTOP=T(KZ)
         END IF

         PTOP=P(KZ)
         UTOP=U(KZ)
         VTOP=V(KZ)
         WTOP=W(KZ)

!        convert to rh fraction
         IF(QFLG)THEN
            ESAT=EXP(21.4-(5351.0/T(KZ)))
            RTOP=Q(KZ)*PTOP/(0.622*ESAT)
         ELSE
            RTOP=Q(KZ)/100.0
         END IF

!        convert the normalized height to actual agl of level
         ZTOP=PSG(KZ)*(1.0-ZSFC/ZMDLT)

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
            PP(KL)=FRAC*(PTOP-PBOT)+PBOT

!           density and potential temperature from local value
            DEN(KL)=P2JM*PP(KL)/(TT(KL)*RDRY)

! CHG(10/07/03) use different constants for RAMS (consistent with rconstants.h)
            IF(.NOT.RAMSFLG)TT(KL)=TT(KL)*(1000.0/PP(KL))**0.286
            IF(RAMSFLG)TT(KL)=TT(KL)*(1000.0/PP(KL))**(RDRY/1004.0)

!           vertical velocity term converted to sigma/time
            OMEGA=FRAC*(WTOP-WBOT)+WBOT

            IF(ZMDLT.LE.20000.0)THEN
!              assume data from RAMS where omega in mb/sec
!              and then ds/dt = omega / (row g Ztop)
              WW(KL)=P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))
            ELSE
!              assume data from coamps where omega in m/sec
!              and then ds/dt = - dz/dt / Ztop
              WW(KL)=-OMEGA/(ZMDL-ZSFC)
            END IF
! JCL(03/27/03) simple linear interpolation; convert to dsigma/dt
! CHG but don't apply density normalization yet, do this when using winds locally
! WWND corresponds to "w*-flux", not dZ/dt!
! w* = dz*/dt, z*=Zt * (z - Zsurf)/(Zt-Zsurf)= zagl*Zt/(Zt-Zsurf)
! i.e. z* = Zt * (1.0 - sigma_z)
!            IF(RAMSFLG)WW(KL)=-1.0*(FRAC*(WTOP-WBOT)+WBOT)/(ZMDL-ZSFC)
            IF(RAMSFLG)WW(KL)=-1.0*(FRAC*(WTOP-WBOT)+WBOT)/ZMDL
! TEST
!            IF(RAMSFLG)WW(KL)=-1.0*(FRAC*(WTOP-WBOT)+WBOT)/(ZMDL)
!           write(45,*)'prfter:',RAMSFLG,KL,WW(KL)
            KL=KL+1
            IF(KL.GT.NL)RETURN
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KL))
!           adjust for terrain compression
!            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)
         END DO

         PBOT=PTOP
         ZBOT=ZTOP
         RBOT=RTOP
         TBOT=TTOP
         UBOT=UTOP
         VBOT=VTOP
         WBOT=WTOP

!         IF(KZ.EQ.2)WRITE(45,*)'zb,zsfc:',ZTOP,ZSFC
!         IF(KZ.EQ.1)WRITE(45,*)'zlvl,frac:',ZLVL,frac

      END DO

!     sounding ends but levels remain to be filled
      IF(KL.LE.NL)THEN
         IF(KFLAG.EQ.0)THEN
            WRITE(30,*)'WARNING prfter: data extrapolation'
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KL))
!           adjust for terrain compression
!            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)
            WRITE(30,*)'from level (m): ',ZLVL
            KFLAG=1
         END IF

         DO KK=KL,NL
! CHG(09/10/03) correct transformation between sigma and agl
!            ZLVL=ZMDL*(1.0-ZSG(KK))
!           adjust for terrain compression
!            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
            ZLVL=(1.0-ZSG(KK))*(ZMDL-ZSFC)
            DELZ=ZLVL-ZZ(KK-1)
            ZZ(KK)=ZLVL
            TT(KK)=TT(KK-1)
            RH(KK)=RTOP

!           use previous pressure and density to find pressure
            PP(KK)=PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM

!           convert potential back to local temperature
! CHG(10/07/03) use different constants for RAMS (consistent with rconstants.h)
            IF(.NOT.RAMSFLG)TEMP=TT(KK)*(PP(KK)/1000.0)**0.286
            IF(RAMSFLG)TEMP=TT(KK)*(PP(KK)/1000.0)**(RDRY/1004.0)
            DEN(KK)=P2JM*PP(KK)/(TEMP*RDRY)
            UU(KK)=UU(KK-1)
            VV(KK)=VV(KK-1)

!           diminish magnitude when no data
            WW(KK)=WW(KK-1)/2.0

         END DO
      END IF

      RETURN
      END
