!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBANL           STaBility ANaLysis from meteo profiles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY ANALYSIS FROM MODEL SURFACE FLUX PARAMETERS
!   OR FROM PROFILE IF THOSE OUTPUTS ARE NOT AVAILABLE.
!   FOR COMPUTATIONAL PURPOSES LEVEL 2 IS THE TOP OF SURFACE
!   SURFACE LAYER.  ROUTINE RETURNS THE HEIGHT NORMALIZED
!   MONIN-OBUKHOV LENGTH, FRICTION TERMS, AND MIXED LAYER DEPTH.
!
! PROGRAM HISTORY LOG:
!   Last Revision: 01 Jul 1997 - RRD
!
! USAGE:  CALL STBANL(KS,SFLG,EFLG,STAR,Z0,FMU,FMV,FHS,NL,
!              UU,VV,TT,ZZ,DEN,USTR,TSTR,WSTR,ZMIX,SLEN,PSI)
!   INPUT ARGUMENT LIST:
!     KS    - int index that indicates top of surface layer
!     SFLG  - log surface fluxes available (heat & momentum)
!     EFLG  - log momentum flux as scalar exchange
! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!     EFLX      - log   momentum flux as scalar exchange
!     HFLX,UFLX - log   surface fluxes available (heat & momentum)
!     USTAR,TSTAR       - log   friction velocity and temperature available
!     Z0    - real      roughness length (m)
!     FMU   - real      u-component momentum flux (N/m2) or ...
!                 scalar exchange coefficient (kg/m2-s) EFLG=TRUE
!                 friction velocity (m/s) when STAR=TRUE
!     FMV   - real      v-component momentum flux
!     FHS   - real      sensible heat flux (W/m2)
!                 friction temperature  when STAR=TRUE
!     NL        - int   number of output sigma levels
!     UU,VV - real      array horizontal wind components (m/s)
!     TT    - real      array virtual potential temperature (deg K)
!     ZZ    - real      array height at levels (m)
!     DEN   - real      array air density (kg/m3)
!   OUTPUT ARGUMENT LIST:
!     USTR  - real      friction velocity (m/s)
!     TSTR  - real      friction temperature (deg K)
!     WSTR  - real      convective velocity scale (m/s)
!     ZMIX  - real      mixed layer depth (m)
!     SLEN  - real      Obukhov stability length (m)
!     PSI   - real      integrated stability function for heat
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: stbanl.f90,v 1.9 2007-02-16 17:53:48 tnehrkor Exp $
!
!$$$

! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!       SUBROUTINE STBANL(KS,SFLG,EFLG,STAR,Z0,FMU,FMV,FHS,NL,
!     :   UU,VV,TT,ZZ,DEN,USTR,TSTR,WSTR,ZMIX,SLEN,PSI)
! CHG:(11/19/01) add humidity as input argument, to allow calculation of LCL, LFC and LOC
! CHG:(12/04/01) add ZLOC (limit of convection) as output argument
! CHG:(9/17/02) add 'ICONVECT' as convection flag
      SUBROUTINE STBANL(KS,EFLX,HFLX,UFLX,USTAR,TSTAR,Z0,FMU,FMV,FHS,NL,&
     &   UU,VV,TT,ZZ,DEN,RH,USTR,TSTR,WSTR,ZMIX,SLEN,PSI,ZLOC,ICONVECT)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)

! CHG:(9/13/02) out of bound error for local arrays, use max size
!     array sizes

      REAL*8   UU(NL),VV(NL),TT(NL),ZZ(NL),DEN(NL),RH(NL)

! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!      LOGICAL SFLG,EFLG,ZFLG,STAR
      LOGICAL EFLX,HFLX,UFLX,USTAR,TSTAR,ZFLG

!     gravity (m/s2)       Karman's    heat (J/Kg-K)
      DATA GRAV/9.80616/,  VONK/0.4/,  CP/1005.0/

!     coefficients for Beljaars-Holtslag data (1991)
      DATA A/1.0/, B/0.66667/, C/5.0/, D/0.35/, PR/0.923/

!     coefficients for phi Ri => z/L valid at 75m over land
      DATA A1/0.005/, B1/41.2/, C1/1.185/, D1/1.5/, E1/1.37/, F1/0.50/

!     coefficients for psi integral
      DATA P1/0.1164E-04/, P2/-2.7188/, P3/-2.1551/, P4/-0.9859/,       &
     &     P5/-0.1765/

! CHG:(11/19/01)
!     dry air (J/Kg-K)   mb to j/m3
      DATA RDRY/287.04/,      P2JM/100.0/

! CHG:(12/01/01) needed for convective parameterization
!      REAL*8   TK(NL),PRS(NL),E(NL),ES(NL),R(NL),RS(NL),TP(NL)
!      REAL*8   TDRY(NL),TPV(NL),TEV(NL)
! CHG:(9/13/02) out of bound error for local arrays, use max size
      REAL*8   TK(NZM),PRS(NZM),E(NZM),ES(NZM),R(NZM),RS(NZM),TP(NZM)
      REAL*8   TDRY(NZM),TPV(NZM),TEV(NZM)

      REAL(KIND(1d0)), EXTERNAL :: RTSAFE, RSAT_T

!---------------------------------------------------------------------------------------------------
!=>analyze sounding to exclude extrapolated data

!     go up the sounding to find the first level of real data
!     prf??? routines set height to -height if internal model levels
!     are below the first data level.  Prevents using interpolated
!     levels for stability analysis

      KDAT=1
      DO WHILE (ZZ(KDAT).LT.0.0)
!        convert back to positive value and save index
         ZZ(KDAT)=-ZZ(KDAT)
         KDAT=KDAT+1
      END DO
      IF(KDAT.GE.NL)THEN
         WRITE(*,*)'Error: stbanl - no observed data'
         STOP
      END IF

!     level used as top point for sfc layer stability calculations
      KDAT=MAX0(2,KS,KDAT)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG:(11/19/01) BUG, TT is NOT virtual pot. Temp, but only pot Temp
!     at least for EDAS, FNL, (not NGM)

!=>find mixed layer depth as the height T(Zi) > 2 + Tsfc
! CHG(8/11/01): replace by more detailed algorithm, at end of routine
!      ZFLG=.FALSE.
!      ZMIX=ZZ(NL)
!      DO K=NL,1,-1
!         IF(.NOT.ZFLG.AND.TT(K).GT.TT(1)+2.0)THEN
!            ZMIX=ZZ(K)
!         ELSE
!            ZFLG=.TRUE.
!         END IF
!      END DO
!     minimum mixed layer height (day or night)
!      ZMIX=DMAX1(ZMIX,DBLE(250.))


!=>determine z/L from surface flux fields
! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!      IF(SFLG)THEN
      IF(UFLX.OR.USTAR)THEN

         IF(USTAR)THEN
!           friction velocity available as input field
            FVEL=DMAX1(FMU,DBLE(0.0001))

! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!        only scalar exchange coefficient available
!         IF(EFLG)THEN
         ELSEIF(EFLX)THEN

! JCL:
!           take ABS of FMU b/c VERY occasionally (~10 out of 100000)
!           FMU is NEGATIVE. . .then model crashes b/c taking SQRT of NEG #
!            USTR=DMAX1(DBLE(0.01)
!     &           ,DSQRT(DABS(FMU)*DABS(UU(KDAT))/DEN(KDAT)))
!            VSTR=DMAX1(DBLE(0.01)
!     &           ,DSQRT(DABS(FMU)*DABS(VV(KDAT))/DEN(KDAT)))
!           friction velocity (m/s)
!            USTR=AMAX1(0.01,DSQRT(FMU*DABS(UU(KDAT))/DEN(KDAT)))
!            VSTR=AMAX1(0.01,DSQRT(FMU*DABS(VV(KDAT))/DEN(KDAT)))
! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!           only scalar exchange coefficient available
            USTR=DMAX1(DBLE(0.01)                                       &
     &           ,DSQRT(DABS(FMU)*DABS(UU(KDAT))/DEN(KDAT)))
            VSTR=DMAX1(DBLE(0.01)                                       &
     &           ,DSQRT(DABS(FMU)*DABS(VV(KDAT))/DEN(KDAT)))
            USTR=DSIGN(USTR,UU(KDAT))
            VSTR=DSIGN(VSTR,VV(KDAT))
!           scalar friction velocity
            FVEL=DSQRT(USTR*USTR+VSTR*VSTR)

!        both u and v momentum fluxes are available
!         ELSE
!            USTR=DSIGN(DMAX1(DBLE(0.01)
!     &        ,DSQRT(DABS(FMU)/DEN(KDAT))),FMU)
!            VSTR=DSIGN(DMAX1(DBLE(0.01)
!     &        ,DSQRT(DABS(FMV)/DEN(KDAT))),FMV)
!         END IF

! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
         ELSEIF(UFLX)THEN
!           both u and v momentum fluxes are available
            USTR=DSIGN(DMAX1(DBLE(0.01)                                 &
     &        ,DSQRT(DABS(FMU)/DEN(KDAT))),FMU)
            VSTR=DSIGN(DMAX1(DBLE(0.01)                                 &
     &        ,DSQRT(DABS(FMV)/DEN(KDAT))),FMV)
!           scalar friction velocity
            FVEL=DSQRT(USTR*USTR+VSTR*VSTR)
         END IF


! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!        scalar friction velocity
!         FVEL=DSQRT(USTR*USTR+VSTR*VSTR)
!
!        friction temperature from sensible heat flux
!         FHS=DMAX1(DBLE(-25.),FHS)
!         TSTR=DSIGN(DMAX1(DBLE(0.0001)
!     &     ,DABS(FHS/(DEN(KS)*CP*FVEL))),-FHS)


! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
         IF(TSTAR)THEN
!           friction temperature available
            TSTR=DSIGN(DMAX1(DBLE(0.0001)                               &
     &        ,DABS(FHS)),FHS)
         ELSE
!           friction temperature from sensible heat flux
            FHS=DMAX1(DBLE(-25.0),FHS)
            TSTR=DSIGN(DMAX1(DBLE(0.0001)                               &
     &        ,DABS(FHS/(DEN(KS)*CP*FVEL))),-FHS)
         END IF

!        normalized Monin-Obukhov length (for top of sfc layer)
         ZL=ZZ(KS)*VONK*GRAV*TSTR/(FVEL*FVEL*TT(KS))

!        replace component friction velocity with scalar
         USTR=FVEL


! JCL:(5/7/01) implement Draxler's 'generalized flux computation' changes made in July 28, 2000
!=>normalized profiles friction temperature and velocity available
!
!      ELSEIF(STAR)THEN
!
!        copy to normal names that are returned
!         USTR=FMU
!         TSTR=FHS
!
!        normalized Monin-Obukhov length (for top of sfc layer)
!         ZL=ZZ(KS)*VONK*GRAV*TSTR/(USTR*USTR*TT(KS))


!=>determine z/L from wind and temperature sounding

      ELSE
!        Bulk Richardson number uses only non-extrapolated levels
         TBAR=(TT(KDAT)+TT(1))/2.0
         DELZ=ZZ(KDAT)-ZZ(1)
         DELT=TT(KDAT)-TT(1)
! JCL:   (5/19/00)remove **2--may run faster
         DELU=(UU(KDAT)-UU(1))*(UU(KDAT)-UU(1))+                        &
     &        (VV(KDAT)-VV(1))*(VV(KDAT)-VV(1))
!         DELU=(UU(KDAT)-UU(1))**2+(VV(KDAT)-VV(1))**2
         DELU=DMAX1(DBLE(0.0001),DELU)
         RI=GRAV*DELT*DELZ/TBAR/DELU

!        correction for excessive height (see Golder Ri~Z^2)
!        when data level above assumed surface layer top
         DZ=ZZ(KDAT)*ZZ(KDAT)/ZZ(KS)/ZZ(KS)
         RI=DMAX1(DBLE(-1.),DMIN1(RI/DZ,DBLE(1.)))

!==>convert bulk Ri to z/L at level KS using Hess formula

         S=DLOG(ZZ(KS)/Z0+1.0)

!        roughness length for heat
         Z0H=0.1*Z0
         T=DLOG(ZZ(KS)/Z0H +1.0)
         V=DLOG(Z0/Z0H)

         IF(RI.GT.0.AND.RI.LT.0.08)THEN
            ZL=(-T+10.0*S*RI+DSQRT(T*T-20.0*S*T*RI+20.0*S*S*RI))/       &
     &         (10.0*(1.0-5.0*RI))
         ELSEIF(RI.GE.0.08)THEN
            ZL=(A1*S+B1)*RI*RI+(C1*S-D1*V-E1)*RI
         ELSE
            ZL=RI*(S*S/T-F1)
         END IF
         ZL=DSIGN(ZL,RI)

         IF(ZL.GE.0.0)THEN
!           stable surface layer Beljaars-Holtslag
            ETRM=B*DEXP(-D*ZL)*(1.0-D*ZL+C)
            PHIM=1.0+(A+ETRM)*ZL
            PHIH=PR*(1.0+(A*SQRT(1.0+A*B*ZL)+ETRM)*ZL)
         ELSE
!           unstable surface layer Betchov-Yaglom / Kadar-Perepelkin
            PHIM=((1.0+0.625*ZL*ZL)/(1.0-7.5*ZL))**0.3333
            PHIH=0.64*((3.0-2.5*ZL)/(1.0-10.0*ZL+50.0*ZL*ZL))**0.3333
         END IF

!        compute friction terms
         USTR=VONK*ZZ(KS)*SQRT(DELU)/PHIM/DELZ
         TSTR=VONK*ZZ(KS)*DELT/PHIH/DELZ

!        recompute Z/L from friction terms to be consistent
         ZL=ZZ(KS)*VONK*GRAV*TSTR/(USTR*USTR*TT(KS))
      END IF


!=>check limits on Z/L (-2 <= z/L <= 15 )

      ZL=DMAX1(DBLE(-2.0),DMIN1(DBLE(15.),ZL))
      IF(ZL.EQ.0.0)ZL=0.001
      SLEN=ZZ(KS)/ZL

!=>compute integral (psi) required for deposition calculations

      IF(ZL.LT.-0.001)THEN
!        unstable
         PSI=P1+ZL*(P2+ZL*(P3+ZL*(P4+P5*ZL)))
      ELSEIF(ZL.GE.-0.001.AND.ZL.LT.0.0)THEN
!        neutral
         PSI=-2.7283*ZL
      ELSE
!        stable
         PSI=-(1.0+A*B*ZL)**1.5-B*(ZL-C/D)*EXP(-D*ZL)-B*C/D+1.0
      END IF

! CHG(8/11/01): replace by more detailed algorithm
!      use excess temperature as function of virt. pot. tmp
!      flux, Ustr and Wstr after Holtslag&Boville

      ZFLG=.FALSE.
      ZMIX=ZZ(1)
      ZMIXNEW=ZZ(1)
      RI=0.0
      DO K=2,NL
! CHG(02/18/04) back to Ri=0.25 as crit. Rich. number
         IF (.NOT.ZFLG.AND.RI.LT.0.25) THEN
! CHG(11/26/02) use Ri=0.5 as crit. Rich. number
!         IF(.NOT.ZFLG.AND.RI.LT.0.5)THEN
            ZMIX=ZZ(K)
            DELZ=ZZ(K)-ZZ(1)
!           calculate excess temp for convective cases
            IF(ZL.LT.0.0)THEN
!             surface layer temperature scale
              TSTK=TSTR*DEN(KS)/DEN(K)
              WSTR=ABS(GRAV*USTR*TSTK*ZMIX/TT(1))**0.3333
              TEX=DABS(8.5*USTR*TSTK/                                   &
     &          (USTR*USTR*USTR+0.6*WSTR*WSTR*WSTR)**0.3333)
            ELSE
              TEX=0.0
            END IF
            TT1=TT(1)+TEX
            DELT=TT(K)-TT1
! CHG(02/18/04) back to version with surface U,V = UU(1) etc.
! CHG(11/1902) Try using surface with U,V==0
            DELU=(UU(K)-UU(1))*(UU(K)-UU(1))+                           &
     &           (VV(K)-VV(1))*(VV(K)-VV(1))+                           &
     &           100*USTR*USTR
!            DELU=UU(K)*UU(K)+
!     &           VV(K)*VV(K)+
!     &           100*USTR*USTR
            DELU=DMAX1(DBLE(0.0001),DELU)
            RIOLD=RI
            RI=GRAV*DELT*DELZ/TT1/DELU
!           linear interpolation to altitude where Ri=0.25
! CHG(02/18/04) back to Ri=0.25 as crit. Rich. number
            ZMIX=(ZZ(K)-ZZ(K-1))/(RI-RIOLD)*(0.25-RIOLD)+ZZ(K-1)
! CHG(11/26/02) CHANGE TO Ri=0.5
!            ZMIX=(ZZ(K)-ZZ(K-1))/(RI-RIOLD)*(0.5-RIOLD)+ZZ(K-1)
! CHG(10/30/02): no linear interpolation to altitude where Ri=0.25,
! but use NEAREST level
            ZMIXNEW=ZZ(K)
! CHG(10/30/02): always use next higher level
! CHG(02/25/04): use nearest level to remove bias betw. mod and observations of zi
            IF(DABS(ZZ(K)-ZMIX).GE.DABS(ZZ(K-1)-ZMIX))THEN
              ZMIXNEW=ZZ(K-1)
            END IF
         ELSE
            ZFLG=.TRUE.
         END IF
      END DO

      !  check wether a 'good' mixing layer height has been found
      IF (ZMIXNEW >= ZZ(NL-1)) THEN
         !(050725) temporary for WRF test 
         !(050725) STOP 'Subroutine stbanl: Couldn''t find ZMIX. Wrong met data? STOP.'
         ZMIXNEW=ZZ(NL-1)
         WRITE(*,*)'ZMIXNEW=',ZMIXNEW
      END IF

!     minimum mixed layer height (day or night)
! CHG:use ZZ(2) as minimum, more realistic during night
! ZZ(2) is around 75 meters above ground
!      ZMIX=DMAX1(ZMIXNEW,ZZ(2))
! use level 3 as minimum, ZZ(3) is around 202 meters above ground
      ZMIX=DMAX1(ZMIXNEW,ZZ(3))

!=>compute convective velocity scale
      IF(ZL.GE.0.0)THEN
         WSTR=0.0
      ELSE
         WSTR=ABS(GRAV*USTR*TSTR*ZMIX/TT(KS))**0.3333
      END IF


! CHG:(9/17/02) add 'ICONVECT' as convection flag
      IF(ICONVECT.EQ.0)THEN
         ZLOC=-999.0
      ELSE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCC    FOR CONVECTIVE REDISTRIBUTION   CCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG:(11/19/01) compute lifting condensation level
!     First: get TLCL, then get PLCL, then first level with P<PLCL

!     get TLCL from surface temperature and humidity, after STULL
!     from TT (potential temperature) to TK (temperature in K)
!     (There is a bug in prfsig.f etc, no virtual temperatures
!     calculated for met data with humidity in relative humidity,
!     only for spec. humidity)
!     use parcel starting at lowest level (1), around 10m, has
!     temperature closest to T02M (EDAS and FNL)

!     initialize
!     lifting condensation level
      LCL=NL
!     level of free convection
      LFC=NL
!     limit of convection
      LOC=1
!     index for level just above zmix
      KMX=NL

!     loop over all levels, from bottom to top
      DO K=1,NL
!CCCC get p, T(K), e, esat, and rsat
!     get pressure (mbar), from TT (K) and DEN (kg/m3)
!     p=DEN*TK*RDRY/P2JM; TK=TT*(p/1000)**0.286
!     i.e. p=DEN*TT*(p/1000)**0.286*RDRY/P2JM
         PRS(K)=(DEN(K)*TT(K)*(1.0/1000.0)**0.286*RDRY/P2JM)            &
     &        **(1/(1.0-0.286))
!     get T (K) from pot. temp. (K)
         TK(K)=TT(K)*(PRS(K)/1000)**0.286

!     get saturation vap. pres. (mbar) from T (K)
         IF(TK(K).GE.273.15)THEN
!     liquid
            ES(K)=6.1078*exp(17.29694*(TK(K)-273.16)/(TK(K)-35.86))
         ELSE
!     ice
            ES(K)=exp(23.33086-6111.72784/TK(K)+0.15215*log(TK(K)))
         END IF

!     from RH (rel.humidity fraction) to vapor pressure (mbar)
         E(K)=RH(K)*ES(K)

!     get spec. hum. r (g/g) from e (mbar) and p (mbar)
         R(K)=0.622*E(K)/(PRS(K)-E(K))

!     get sat. spec. hum. rsat (g/g) from e (mbar) and p (mbar)
         RS(K)=0.622*ES(K)/(PRS(K)-ES(K))

!     get dry adiabate ( T in K of dry adiab. lifted parcel...)
         TDRY(K)=TK(1)* (PRS(K)/PRS(1))**0.286

!CCCC end get p, TK(K), e, esat, and rsat

!CCCCCCC get lifting condensation level
         IF(K.EQ.1)THEN
!     get lifting condensation Temperature (K) from level (1) vap. pres. (mbar) and T (K)
            TLCL=2840.0/(3.5*DLOG(TK(1))-DLOG(E(1)/10.0)-7.108)+55.0
!     get lifting condensation pressure (mbar) from p (mbar) and TLCL (K)
            PLCL=PRS(1)*(TLCL/TK(1))**3.5
         END IF
!     find level LCL just above condensation
         IF(TDRY(K).GE.TLCL)THEN
            LCL=K+1
         END IF
!CCCCCCC end get lifting condensation level

!CCCCCC get parcel temperature TP (cloud temperature) (in K)
!     starting with parcel from lowest model layer (1);

!     below tlcl: dry adiabat;
         IF(TDRY(K).GT.TLCL)THEN
            TP(K)=TDRY(K)
!     Use specific humidity of parcel to get virt. temp
            TPV(K)=TP(K)*(1+0.61*R(1))
         END IF

!     above TLCL: wet adiabat; use pressure levels from model
         IF(TDRY(K).LE.TLCL)THEN

!     get saturation equivalent potential temperature
!     1. of parcel at true lifting condensation level
            TTES1=TT(1)+R(1)*2500.0*TT(1)/TLCL
!     2. of parcel at level K above LCL, as function of TP
!     (which will be the temperature needed to match the
!     saturation equivalent potential temperature in cloud)
!     rtsafe() is used to find TP above LCL
!     Max. value, needed for iteration (min. = TDRY)
            TMAX=TDRY(K)+30.0
            T_ACC=0.01
            TP(K)=rtsafe(TDRY(K)-1,TMAX,PRS(K),TTES1,T_ACC)
!     Use sat. specific humidity of parcel to get virt. temp
            TPV(K)=TP(K)*(1+0.61*RSAT_T(PRS(K),TP(K)))
         END IF
!CCCCCC end get parcel temperature

!     Get virtual temperature of environment
         TEV(K)=TK(K)*(1+0.61*R(K))

!     Check whether parcel is buoyant:
!     first time: LFC (level of free convection)
!     last time (of buoyant layer): LOC (limit of convection)
!     level of free convection
         IF(TPV(K).GT.TEV(K).AND.LFC.EQ.NL)LFC=K
!     level of free convection
         IF(TPV(K).LT.TEV(K).AND.LFC.LT.NL &
     &      .AND.K.GT.LFC.AND.LOC.EQ.1) LOC=K

!     index for level just above zmix
         IF(ZZ(K).GT.ZMIX.AND.KMX.EQ.NL)KMX=K+1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCC for testing
!         IF(K.EQ.1)THEN
!            WRITE(45,*) 'LVL P T TP R RS TEV TPV, DTV'
!         END IF
!         WRITE(45,*)K,PRS(K),TK(K),TP(K),R(K),RS(K),TEV(K),TPV(K),
!     &   TPV(K)-TEV(K)
!         IF(K.EQ.NL)THEN
!            WRITE(45,*)'LCL LFC LOC KMX ZLCL ZLFC ZLOC ZKMX'
!            WRITE(45,*)LCL,LFC,LOC,KMX,ZZ(LCL), ZZ(LFC),
!     &      ZZ(LOC), ZZ(KMX)
!         END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCC end for testing

      END DO
!     end loop over all levels, from bottom to top

!     assign altitude for limit of convection
      ZLOC=ZZ(LOC)

!     check if positively buoyant area is not too far above zi,
!     otherwise change ZLOC to neg. value
      IF(ZZ(LFC).GT.ZMIX+2000.0)LOC=1

!     check if no positively buoyant area was found
      IF(LOC.EQ.1)ZLOC=-999.0

!     check if ZLOC is just above ZMIX (when levels are equal)
      IF(LOC.LE.KMX)ZLOC=-999.0

!     TEST ONLY
!     calculate enhancement factor for conv. redistribution area
!     multiplied with cloud fraction, this gives the fraction to be vertically redistributed
!     Get cloud base mass flux in kg/m2
!     Assume 2 m/s updraft averaged over cloudbase
!     Assume 3 h duration (repetition cycle for met fields)
!      FMCB=2.0*3.0*3600.0*DEN(LCL)
!     get integrated density up to cloud base, in kg/m2
!      DENT=0.0
!      DO K=1,LCL-1,1
!         DENT=DENT+DEN(K)*(ZZ(K+1)-ZZ(K))
!      END DO
!     get enhancement factor as ratio FMCB/DENT
!      FENH=FMCB/DENT
!            WRITE(45,*)'FMCB DENT FENH'
!            WRITE(45,*)FMCB,DENT,FENH
! CHG:(12/05/01) calculated FENH are in the range of 8-60, so redistribute whole gridcell independent of cloud cover

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCC END FOR CONVECTIVE REDISTRIBUTION   CCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG: (9/17/02) end of IF(ICONVECT.EQ.1)
      END IF

      RETURN
      END
