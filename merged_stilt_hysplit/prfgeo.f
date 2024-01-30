!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! $Id: prfgeo.f,v 1.5 2013/04/09 20:08:20 trn Exp $
! SUBPROGRAM:  PRFGEO           processes GEOS profiles on pressure sigma
!   PRGMMR:    Thomas Nehrkorn   ORG: AER, Inc      DATE: 2013/04
!   based on code PRFSIG written by: ROLAND DRAXLER   ORG: R/E/AR
!
! ABSTRACT:  THIS CODE WRITTEN AT AER
!   PROFILE SIGMA CONVERTS A SOUNDING ON PRESSURE-SIGMA COORDINATES
!   (TYPE 1) TO MODEL TERRAIN FOLLOWING SIGMA-Z LEVEL.
!   INTERPOLATES ALL METEO VARIABLES TO MODEL VERTICAL GRID.
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

SUBROUTINE PRFGEO(VMIX,TKEN,VELV,ZFLG,QFLG,UFLG,TFLG,PFLG,DZDT,  &
                  ZSFC,P0,U0,V0,T0,Z0,ap0,NZ,PSG,P,U,V,W,T,Q,zoraq,ap,azort,NL,ZMDL, &
                  ZSG, PP,UU,VV,WW,TT,ZZ,RH,DEN,AA)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables 
!-------------------------------------------------------------------------------

  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: zflg     ! pressure available (nonhydrostatic)
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: dzdt     ! vertical velocity in height units
  REAL,       INTENT(IN)    :: zsfc     ! terrain height (m)
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  REAL,       INTENT(INOUT) :: ap0      ! surface pressure (time-averaged) at data terrain (mb)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: p  (:)   ! pressure data (non-hydrostatic)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: q  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: zoraq(:) ! height or time-averaged spef hum (not: turbulent kinetic energy (m2/s2))
  REAL,       INTENT(IN)    :: ap (:)   ! pressure (time-averaged) (not: u-component velocity var (m2/s2))
  REAL,       INTENT(IN)    :: azort(:) ! time-averaged height or temperature (not: w-component velocity var (m2/s2))
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL,    PARAMETER :: grav  = 9.80616    ! gravity (m/s2)
  REAL,    PARAMETER :: rdry  = 287.04     ! dry air (J/Kg-K)
  REAL,    PARAMETER :: p2jm  = 100.0      ! mb to j/m3
  INTEGER            :: kflag  = 0         ! diagnostic  

  REAL               :: zlvl
  REAL               :: pbot,zbot,rbot,tbot,ubot,vbot,wbot,ebot,hbot,xbot
  REAL               :: ptop,ztop,rtop,ttop,utop,vtop,wtop,etop,htop,xtop
  REAL               :: frac,abot,omega,atop,esat,delz,tbar
  real               :: apbot, aptop, azbot, aztop, adelz, afrac !for time-averaged pressures/heights
  INTEGER            :: kk,kz,kl
  logical            :: have_z
  logical :: lprint  !for debugging type output
  integer,save :: iprint=0, icall=0, iskip=100  !for debugging type output
  integer,save :: maxprint=100  !for debugging type output

  SAVE KFLAG

!-------------------------------------------------------------------------------
! compute the log of the surface pressure
!-------------------------------------------------------------------------------

  if (tken .or. velv) then
     write (*,*) 'PRFGEO: unexpected .TRUE. for tken and/or velv, conflicts with usage of e/h/x arrays'
     STOP 'PRFGEO: unexpected .TRUE. for tken and/or velv'
  end if
  
  if (.not. zflg) then
     write (*,*) 'PRFGEO: unexpected .FALSE. for ZFLG: pressure values are required'
     STOP 'PRFGEO: unexpected .FALSE. for ZFLG'
  end if
  
  have_z = maxval(zoraq) .gt. 200.  !this will evaluate to .FALSE. if average SPHU or RH is stored in zoraq

  icall = icall + 1
  lprint = icall .eq. 1 .or. mod(icall,iskip) .eq. 0
  if (lprint) iprint = 1+iprint
  if (lprint .and. iprint .le. maxprint) &
       & WRITE(KFJCLMSG,'(a,i8,a/(a,l10))') 'Entered prfgeo with icall=',icall, &
       & ' and flags for','pressures (zflg)=',zflg, &
       & 'p0 (pflg)=',pflg,'u,v at sfc (uflg)=',uflg,'t0 (tflg)=',tflg, &
       & 'q (vs rh, QFLG)=',qflg,'have_z=',have_z

  IF(PFLG)THEN
!    surface pressure from input data
     PBOT=LOG(P0)
     aPBOT=LOG(aP0)
  ELSE
!    use lowest upper level temp to estimate surface pressure
!    assume that mslp always 1013 hPa
     PBOT=LOG(1013.0)-ZSFC*GRAV/RDRY/T(1)
     P0=PBOT
     ap0=pbot
  END IF

!-------------------------------------------------------------------------------
! use adiabatic profile to estimate surface temperature
!-------------------------------------------------------------------------------

  IF(TFLG)THEN
!    available use low level value
     TBOT=T0
  ELSE
!    drop adiabatically to surface
     TBOT=T(1)*PSG(1)**(-0.286)
  END IF
! convert all temperatures to virtual (rrd - 12/16/2002)
  IF(QFLG)TBOT=TBOT*(1.0+0.61*Q(1))

!-------------------------------------------------------------------------------
! convert specific humidity to fraction of RH% to fraction
!-------------------------------------------------------------------------------

  IF(QFLG)THEN
     ESAT=EXP(21.4-(5351.0/T(1)))
     RBOT=Q(1)*P0/(0.622*ESAT)
  ELSE
     RBOT=Q(1)/100.0
  END IF

! Vertical interpolation to the instantaneous pressures/heights
! initial output vertical index
  KL=1
! initial output sigma level
  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
! starting height (agl) and vertical velocity
  ZBOT=0.0

  if (lprint .and. iprint .le. maxprint) then
     WRITE(KFJCLMSG,'(a/(a,g15.6))') 'Input and derived sfc variables:','zsfc=',zsfc, &
          & 'p0 =',p0,'t0=',t0,'tv (tbot)=',tbot, &
          & 'rh (rbot)=',rbot,'u0=',u0,'v0=',v0,'zbot=',zbot,'zmdl=',zmdl
     WRITE(KFJCLMSG,'(a/a5,10a15)') 'Input profile:', &
          & 'k','p','zoraq','z(agl)','temp','tv','q','rh'
  endif

  INSTK_LOOP: DO KZ=1,NZ
     
!    log of pressure at level
     PTOP=LOG(P(kz))

     !    reactivate virtual temperature effects (disabled in prfsig: rrd - 12/16/2002)
     IF(QFLG)THEN
        TTOP=(1.0+0.61*Q(KZ))*T(KZ)
     ELSE
        TTOP=T(KZ)
     END IF

!    use layer average for hypsometric equation
     if (have_z) then
        ztop=zoraq(kz)-zsfc ! Convert from height (MSL) to height (AGL)
        delz=ztop-zbot
     else
        TBAR=0.5*(TTOP+TBOT)
        DELZ=(PBOT-PTOP)*RDRY*TBAR/GRAV
        ZTOP=ZBOT+DELZ
     end if

!    convert to rh fraction
     IF(QFLG)THEN
        ESAT=EXP(21.4-(5351.0/T(KZ)))
        RTOP=Q(KZ)*P(kz)/(0.622*ESAT)
     ELSE
        RTOP=Q(KZ)/100.0
     END IF

     if (lprint .and. iprint .le. maxprint) &
          & WRITE(KFJCLMSG,'(i5,10g15.6)') &
          & kz,p(kz),zoraq(kz),ztop,t(kz),ttop,q(kz),rtop

     INSTKL_LOOP: DO WHILE (ZLVL.LE.ZTOP)
!       height <0 flags extrapolated levels not to be
!       used in stability calculations
        IF(KZ.EQ.1.AND.(.NOT.TFLG).AND.VMIX)THEN
           ZZ(KL)=-ZLVL
        ELSE
           ZZ(KL)=ZLVL
        END IF

!       basic linear interpolation
        FRAC=(ZLVL-ZBOT)/DELZ
        TT(KL)=FRAC*(TTOP-TBOT)+TBOT
        RH(KL)=FRAC*(RTOP-RBOT)+RBOT

!       linear interpolation of log of pressure
        PP(KL)=EXP(FRAC*(PTOP-PBOT)+PBOT)

!       density and potential temperature from local value
        AA(KL)=TT(KL)
        DEN(KL)=P2JM*PP(KL)/(AA(KL)*RDRY)
        TT(KL)=AA(KL)*(1000.0/PP(KL))**0.286

        KL=KL+1
        IF(KL.GT.NL) exit INSTK_LOOP
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
     end DO INSTKL_LOOP
     
     PBOT=PTOP
     ZBOT=ZTOP
     RBOT=RTOP
     TBOT=TTOP
  end DO INSTK_LOOP

!-------------------------------------------------------------------------------
! sounding ends but levels remain to be filled
!-------------------------------------------------------------------------------

  IF(KL.LE.NL)THEN
     IF(KFLAG .lt. 2)THEN
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
        WRITE(KF21,*)'--------------------------------------------------------'
        WRITE(KF21,*)' NOTICE prfgeo: extrapolation from inst level (k,m): ',KL,ZLVL
        WRITE(KF21,*)'Input data levels: ',NZ,'    Internal Sigma levels: ',NL
        WRITE(KF21,*)'--------------------------------------------------------'
        KFLAG=kflag+1
     END IF

     DO KK=KL,NL
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KK))
        DELZ=ZLVL-ZZ(KK-1)
        ZZ(KK)=ZLVL
        TT(KK)=TT(KK-1)
        AA(KK)=AA(KK-1)
        RH(KK)=RTOP

!       use previous pressure and density to find pressure
        PP(KK)=PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM
        DEN(KK)=P2JM*PP(KK)/(AA(KK)*RDRY)
     END DO
  END IF

! Vertical interpolation to time-averaged pressures/heights:

! initial output vertical index
  KL=1
! initial output sigma level
  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
! starting height (agl) and vertical velocity
  AZBOT=0.0
  WBOT=0.0

  if (.not. have_z) then
! Need time-averaged temperatures for hypsometric eqn
! (but use instantaneous T02M for the bottom value if available)
     IF(TFLG)THEN
!    available use low level value
        TBOT=T0
     ELSE
!    drop adiabatically to surface
        TBOT=azort(1)*PSG(1)**(-0.286)
     end if
! convert all temperatures to virtual (rrd - 12/16/2002)
     TBOT=TBOT*(1.0+0.61*zoraq(1))
  END IF

  if (lprint .and. iprint .le. maxprint) then
     WRITE(KFJCLMSG,'(a/a5,10a15)') 'Input profile:', &
          & 'k','ap','azort(msl,K)','az(agl)','u','v','omega'
  endif
  TAVGK_LOOP: DO KZ=1,NZ

!    use layer average for hypsometric equation
     if (have_z) then
        aztop=azort(kz) - zsfc ! Convert from height (MSL) to height (AGL)
        adelz=aztop-azbot
     else
!    log of pressure at level
        aPTOP=LOG(aP(kz))
!    reactivate virtual temperature effects (disabled in prfsig: rrd - 12/16/2002)
        IF(QFLG)THEN
           TTOP=(1.0+0.61*zoraq(KZ))*azort(KZ)
        ELSE
           TTOP=azort(KZ)
        END IF
        TBAR=0.5*(TTOP+TBOT)
        aDELZ=(aPBOT-aPTOP)*RDRY*TBAR/GRAV
        aZTOP=aZBOT+aDELZ
     end if

!    only for the first input data level
     IF(KZ.EQ.1)THEN
!       set low level (ground) wind data if available
        IF(UFLG)THEN
           UBOT=U0
           VBOT=V0
        ELSE
!          use neutral log-law when output below data level
           IF(ZLVL.LT.ZTOP)THEN
              ATOP=LOG(aZTOP/Z0)
              ABOT=LOG(ZLVL/Z0)
              UBOT=U(1)*ABOT/ATOP
              VBOT=V(1)*ABOT/ATOP
           ELSE
              UBOT=U(1)
              VBOT=V(1)
           END IF
        END IF

     END IF

     UTOP=U(KZ)
     VTOP=V(KZ)
     WTOP=W(KZ)

     if (lprint .and. iprint .le. maxprint) &
          & WRITE(KFJCLMSG,'(i5,10g15.6)') &
          & kz,ap(kz),azort(kz),aztop,u(kz),v(kz),w(kz)

     TAVGKL_LOOP: DO WHILE (ZLVL .LE. aZTOP)
!       height <0 flags extrapolated levels not to be
!       used in stability calculations
        IF(KZ.EQ.1.AND.(.NOT.TFLG).AND.VMIX)THEN
           ZZ(KL)=-ZLVL
        ELSE
           ZZ(KL)=ZLVL
        END IF

!       basic linear interpolation
        aFRAC=(ZLVL-aZBOT)/aDELZ
        UU(KL)=aFRAC*(UTOP-UBOT)+UBOT
        VV(KL)=aFRAC*(VTOP-VBOT)+VBOT

!       vertical velocity term converted to sigma/time
!       uses previously computed (inst) density
        OMEGA=aFRAC*(WTOP-WBOT)+WBOT
        IF(.NOT.DZDT)THEN
!          input velocity in pressure units
           WW(KL)=P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))
        ELSE
!          input velocity in height units (some wrf data)
           WW(KL)=-OMEGA/(ZMDL-ZSFC)
        END IF

        KL=KL+1
        IF(KL.GT.NL) exit TAVGK_LOOP
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
     end DO TAVGKL_LOOP

     if (.not. have_z) aPBOT=aPTOP
     aZBOT=aZTOP
     UBOT=UTOP
     VBOT=VTOP
     WBOT=WTOP

  end DO TAVGK_LOOP
  
!-------------------------------------------------------------------------------
! sounding ends but levels remain to be filled
!-------------------------------------------------------------------------------

  IF(KL.LE.NL)THEN
     IF(KFLAG .LT. 2)THEN
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
        WRITE(KF21,*)'--------------------------------------------------------'
        WRITE(KF21,*)' NOTICE prfgeo: extrapolation from tavg level (k,m): ',KL,ZLVL
        WRITE(KF21,*)'Input data levels: ',NZ,'    Internal Sigma levels: ',NL
        WRITE(KF21,*)'--------------------------------------------------------'
        KFLAG=KFLAG+1
     END IF

     DO KK=KL,NL
        UU(KK)=UU(KK-1)
        VV(KK)=VV(KK-1)
!       diminish magnitude when no data
        WW(KK)=WW(KK-1)/2.0
     END DO
  END IF

  if (lprint .and. iprint .le. maxprint) then
     WRITE(KFJCLMSG,'(a/a5,10a15)') 'Output profile:', &
          & 'k','pp','zsg','zz','temp','thetv','rh','u','v','w','den'
     do kk=1,nl
        WRITE(KFJCLMSG,'(i5,10g15.6)') &
             & kk,pp(kk),zsg(kk),zz(kk),aa(kk),tt(kk),rh(kk),uu(kk),vv(kk),ww(kk),den(kk)
     enddo
  endif

END SUBROUTINE prfgeo
