!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! $Id: prfwrf.f,v 1.5 2012/09/25 14:44:16 trn Exp $
! SUBPROGRAM:  PRFWRF           processes PRoFile on WRF mass-based eta
!                              (essentially pressure sigma) levels
!   PRGMMR:    Thomas Nehrkorn   ORG: AER, Inc      DATE: 2005/08/12
!   based on codes PRFSIG and PRFTER written by: ROLAND DRAXLER   ORG: R/E/AR
!
! ABSTRACT:  THIS CODE WRITTEN AT AER, Inc.
!   CONVERTS A SOUNDING ON WRF half- and full eta COORDINATES
!   TO MODEL TERRAIN FOLLOWING SIGMA LEVEL.
!   INTERPOLATES ALL METEO VARIABLES TO MODEL VERTICAL GRID.
!   SEE DOCBLOCK OF PRFPRS FOR A MORE COMPLETE DESCRIPTION
!
! PROGRAM HISTORY LOG:
! $Log: prfwrf.f,v $
! Revision 1.5  2012/09/25 14:44:16  trn
! Fixed bug in extrapolation above model top
!
! Revision 1.4  2012/09/05 16:38:27  trn
! Merge changes from stilt_sigmap branch, at tag SLC_hum_tested
!
! Revision 1.3.2.5  2012/09/05 16:23:36  trn
! Moved compute_zsgfull calls outside module_wrfeta for speedup
!
! Revision 1.3.2.4  2012/08/31 14:43:32  trn
! Bug fixes, replace allocatable by automatic arrays
!
! Revision 1.3.2.3  2012/08/27 20:53:05  trn
! Bug fixes for agreement between old run and vcoord=1 new run
!
! Revision 1.3.2.2  2012/08/09 20:40:39  trn
! Intermediate checkin: successfull compilation and linking
!
! Revision 1.3.2.1  2012/08/09 17:06:59  trn
! Initial set of mods for wrf-sigmap support - before compilation and testing
!
! Revision 1.3  2009/12/04 20:23:55  trn
! Activate diagnostic message
!
! Revision 1.2  2009/11/25 15:48:31  trn
! Bug fix
!
! Revision 1.1  2009/10/26 15:36:55  jel
! first working version on beehive - single meteorology
!
! Revision 1.7  2008/06/04 17:15:31  tnehrkor
! Added comments explaining initialization of DFXUP1_2 etc
!
! Revision 1.6  2006/09/26 14:02:46  tnehrkor
! Decreased diagnostic output
!
! Revision 1.5  2006/07/19 16:41:03  tnehrkor
! Fixed bug for uninitialized values of TKEN_2
!
! Revision 1.4  2006/03/07 16:47:59  tnehrkor
! Fixed array overindexing print statement
!
! Revision 1.3  2006/02/28 15:43:07  tnehrkor
! Removed extraneous if-statement
!
! Revision 1.2  2006/02/17 15:06:23  skoerner
! fixed several incorrect I/O statements, fixed wrong INTENT for variable ZSG
!
! Revision 1.1  2005/12/14 17:05:59  tnehrkor
! Added support for WRF model output; includes a bug fix for dmass computation in prfcom
!
! Revision 1.8  2005/12/12 15:32:46  trn
! Added checks for array overindexing
!
! Revision 1.7  2005/11/04 17:39:17  trn
! Changes for WRF-Grell convective fluxes
!
! Revision 1.6  2005/10/11 17:17:03  trn
! Bug fix
!
! Revision 1.5  2005/10/11 16:36:17  trn
! Changes for vertical staggering of AWRF w
!
! Revision 1.3  2005/10/03 15:16:19  trn
! Changes for mom flux output from WRF, corrected use of awrf winds
!
! Revision 1.2  2005/09/30 18:20:19  trn
! Intermediate changes for wrf winds/fluxes.  Problems: dmass, stblanl for fluxflg=T
!
! Revision 1.1  2005/09/27 13:45:30  trn
! Initial version of AWRF met files: instantaneous winds
!
!
! USAGE:  see below
!   INPUT ARGUMENT LIST: see below
!   OUTPUT ARGUMENT LIST: see below
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

SUBROUTINE PRFwrf(vmix,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,        &
     P0,U0,V0,T0,Z0,w0,alt1,mu,muu,muv,msfu,msfv,msft,fluxflg,NZ,eta, &
     phalf,U,V,W,T,Q,NL,ZMDL,ZSG,zsg_full,dzsg_full,PP,UU,VV,WW,TT,ZZ,RH,DEN, &
     deepflg, shallflg, &
     CFXUP1_1,       &
     CFXUP2_1,CFXDN1_1,DFXUP1_1,         &
     DFXUP2_1,EFXUP1_1,EFXUP2_1,         &
     DFXDN1_1,EFXDN1_1,TKEN_1,           &
     CFXUP1_2,       &
     CFXUP2_2,CFXDN1_2,DFXUP1_2,         &
     DFXUP2_2,EFXUP1_2,EFXUP2_2,         &
     DFXDN1_2,EFXDN1_2,TKEN_2, alt2, vcoord, &
     wrfvcoords, wrfptop, npar_wrfhyb)

  USE funits
  USE module_wrfeta, ONLY: eta2agl
  implicit none 

  !   INPUT ARGUMENT LIST:
  logical, intent(in) :: vmix !ignored
  LOGICAL, intent(in) :: QFLG !flag for mix ratio vs RH data
  LOGICAL, intent(in) :: UFLG,TFLG,PFLG,SFLG !flags for presence of sfc data
  real, intent(inout) :: ZSFC,P0  !these may be set in this routine
  real, intent(in) :: T0,U0,V0,w0,Z0,mu,muu,muv,msfu,msfv,msft
  integer, intent(in) :: nz,nl !number of levels in input and output profiles
  !     ZSFC  - real      terrain height (m)
  !     P0    - real      surface pressure at data terrain (mb)
  !     U0,V0 - real      low level horizontal wind components
  !     T0    - real      low level temperature (deg K)
  !     Z0    - real      aerodynamic roughness length (m)
  !     mu    - real      dry hydrostatic sfc pressure - ptop (WRF "dry mass" variable) (Pa)

  REAL, intent(in) ::  U(NZ),V(NZ),W(NZ),T(NZ),Q(NZ),eta(nz),PHALF(NZ), ALT1(NZ), ZSG(NL)
  REAL,      INTENT(IN)    :: zsg_full (:)
  REAL,      INTENT(IN)    :: dzsg_full (:)
  real, intent(in) :: zmdl
  !     eta   - real      eta-values at full-levels 2-(nz+1) (0-1)
  !     Phalf - real      pressure at half-levels (mb)
  !     U,V   - real      horizontal wind components (half-levels, 
  !                       already at unstaggered (xy) grid point
  !     W     - real      vertical motion term (dzdt, m/s), at full-level
  !     T     - real      potential temperature (pot K)
  !     Q     - real      relative humidity fraction (0-1) or mix ratio
  !     alt1  - real      dry inverse density for WRF momentum flux decoupling and height computation
  logical, intent(in) :: fluxflg !flag for WRF momentum flux input
  logical, intent(in) :: deepflg, shallflg
  REAL, intent(in) ::  CFXUP1_1(nz),       &
       &      CFXUP2_1(nz),CFXDN1_1(nz),DFXUP1_1(nz),         &
       &      DFXUP2_1(nz),EFXUP1_1(nz),EFXUP2_1(nz),         &
       &      DFXDN1_1(nz),EFXDN1_1(nz),TKEN_1(nz)

  REAL,intent(out) ::   PP(NL),TT(NL),ZZ(NL),RH(NL),DEN(NL),             &
       &       UU(NL),VV(NL),WW(NL)
  !   OUTPUT ARGUMENT LIST:
  !     P0    - real      surface pressure at data terrain (mb)
  !     T0    - real      low level temperaure (deg K)
  !     PP    - real      pressure at sigma level (mb)
  !     UU,VV - real      horizontal wind components
  !     WW    - real      vertical motion term (sigma/time)
  !     TT    - real      virtual potential temperature (pot K)
  !     ZZ    - real      internal model sigma height (meters)
  !     RH    - real      relative humidity fraction (0-1)
  !     DEN   - real      air density (kg/m3)
  REAL,intent(out) ::   CFXUP1_2(nl),       &
       &      CFXUP2_2(nl),CFXDN1_2(nl),DFXUP1_2(nl),         &
       &      DFXUP2_2(nl),EFXUP1_2(nl),EFXUP2_2(nl),         &
       &      DFXDN1_2(nl),EFXDN1_2(nl),TKEN_2(nl), alt2(nl)
  integer, intent(in) :: vcoord !Vertical coordinate code (1:sigma-z)
  real, intent(in) :: wrfvcoords(:,:), wrfptop
  integer, intent(in) :: npar_wrfhyb
  
  real :: dzdt, maxdeep, maxshall
  real, save :: epsmax=1.e-5 !threshold: convective mass fluxes > 0
  real, save :: epstotm=0.05 !threshold: fractional mass conservation violation of convective mass fluxes
  real, dimension(NZ) :: deta_full, deta_half, full_agl, half_agl
  real :: vclvl, vcbot, vctop, sign_eta, pdh_top
  real :: z0l

  integer kflag
  real grav, rdry, p2jm
  SAVE KFLAG

  integer :: k, kl, kz, kk, klw, kzalt, kzw
  real :: pbot, tbot, esat, rbot, zlvl, zbot, wbot
  real :: ptop, ttop, tbar, delz, ztop, ubot, vbot, atop, abot, rtop
  real :: utop, vtop, wtop, frac, omega, temp, alttop, altbot, alt_int
  real :: p1013, pref, frac_alt, frac_cnv, totmas1, totmas2
  real :: CFXUP1_bot,       &
       &      CFXUP2_bot,CFXDN1_bot,TKEN_bot,           &
       &      CFXUP1_top,       &
       &      CFXUP2_top,CFXDN1_top,TKEN_top

  real :: hm(nl), hw(nl+1)
  real :: pfu, pfd, phm, pfm

  logical :: lprint  !for debugging type output
  integer,save :: iprint=0, icall=0, iskip=1, iprcnv=0  !for debugging type output
  integer,save :: maxprint=10  !for debugging type output
!!$      integer,save :: maxprint=0  !for suppressing debugging type output

  !     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
  DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/, p1013/1013.0/, pref/1000./

  !     diagnostic
  DATA KFLAG/0/

  !=>compute height of ground surface if not otherwise available

  icall = icall + 1
  lprint = icall .eq. 1 .or. mod(icall,iskip) .eq. 0
  if (lprint) iprint = 1+iprint
  if (lprint .and. iprint .le. maxprint) &
       & WRITE(KFJCLMSG,'(a,i8,a/(a,l10))') 'Entered PRFwrf with icall=',icall, &
       & ' and flags for','zsfc (sflg)=',sflg, &
       & 'p0 (pflg)=',pflg,'u,v at sfc (uflg)=',uflg,'t0 (tflg)=',tflg, &
       & 'q (vs rh, QFLG)=',qflg,'fluxflg=',fluxflg

  !     Compute eta layer thicknesses, and distance between half-levels:
  !     Note that deta will be negative
  sign_eta=-1.  !use this to convert eta to a coordinate that increases with height
  deta_full(1) = eta(1) - 1.0
  deta_half(1) = 0.5*deta_full(1)
  do kk=2,nz
     deta_full(kk) = eta(kk)-eta(kk-1)
  enddo
  do kk=2,nz
     deta_half(kk) = 0.5*(deta_full(kk-1)+deta_full(kk))
  enddo
  if (vcoord .eq. 1) then
  !     Compute z-sigma mass and w-level heights AGL
     do kl=1,nl
        hm(kl) = (1.0-zsg(kl))*(zmdl-zsfc)
     end do
     hw(1)=0
     do kl=2,nl+1
        hw(kl)=hm(kl-1) + hm(kl-1) - hw(kl-1)
     end do
  endif
  pfd=0.
  DO KZ=1,NZ
     if (npar_wrfhyb .gt. 0) then
        !  now use the alternate hydrostatic eqn used in dyn_em_module_initialize_real.F
        !  and calc_p_rho_phi (module_big_step_utilities_em.F):
        !   dPHI = -alt*p*dLOG(p), where p is dry pressure
        ! c3f*MU + c4f + ptop (note offset for storage of full_eta vs full_wrfvcoords
        !  (For TF coordinates, c3f=eta c4f=0)
        !  can only do this when npar_wrfhyb .gt. 0, because wrfptop is read in from same file
        pfu=wrfptop + wrfvcoords(kz+1,7)*mu+wrfvcoords(kz+1,8)
     endif
     if (kz .eq. 1) then
        !        lowest layer: use lowest layer inverse density alt1:
        if (npar_wrfhyb .gt. 0) then
           pfd=wrfptop + wrfvcoords(1,7)*mu+wrfvcoords(1,8)
           ! c3h*MU + c4h + ptop
           !  (For TF coordinates, c3h=eta_half c4h=0)
           phm=wrfptop + wrfvcoords(1,3)*mu+wrfvcoords(1,4)
           pfm=0.5*(pfd+phm)
           full_agl(1) = phm*log(pfd/pfu)*alt1(1)/grav
           half_agl(1) = pfm*log(pfd/phm)*alt1(1)/grav
        else
        !        lowest layer: use lowest layer inverse density alt1:
           delz = -mu*deta_half(1)*alt1(1)/grav
           half_agl(1) = delz
        !        use WRF hydrostatic equation (eq. 2.9 of SkaKD+05) to compute full-level AGL heights:
           delz = -mu*deta_full(1)*alt1(1)/grav
           full_agl(1) = delz
        endif
     else
        !        upper layers: use averaged inverse density alt1:
        !        average(alt) = 0.5*(deta_full(kz-1)*alt1(kz-1)+deta_full(kz)*alt1(kz))/deta_half(kz)
        if (npar_wrfhyb .gt. 0) then
        ! p from eta, using hybrid WRF vertical coordinates:
           pfd = wrfptop + wrfvcoords(kz,7)*mu+wrfvcoords(kz,8)
           phm=wrfptop + wrfvcoords(kz,3)*mu+wrfvcoords(kz,4)
           pfm=0.5*(pfd+phm)
        ! delz from alt and p*dLog(p):
           full_agl(kz) = full_agl(kz-1) + (1/grav)*phm*log(pfd/pfu)* alt1(kz)
           half_agl(kz) = full_agl(kz-1) + (1/grav)*pfm*log(pfd/phm)* &
             0.5*(deta_full(kz-1)*alt1(kz-1)+deta_full(kz)*alt1(kz))/deta_half(kz)
        else
        !        Note that deta_half(kz) in denominator cancels with factor out front
           delz = -mu*0.5*(deta_full(kz-1)*alt1(kz-1)+deta_full(kz)*alt1(kz))/grav
           half_agl(kz) = half_agl(kz-1) + delz
        !        use WRF hydrostatic equation (eq. 2.9 of SkaKD+05) to compute full-level AGL heights:
           delz = -mu*deta_full(kz)*alt1(kz)/grav
           full_agl(kz) = full_agl(kz-1) + delz
        end if
     endif
  enddo

  IF(.NOT.SFLG)THEN
     if (.not. pflg) then
        WRITE(KFJCLMSG,*) 'prfwrf: SFLG and PFLG are FALSE, cannot recover'
        stop 'prfwrf: SFLG and PFLG are FALSE, cannot recover'
     endif
     !         computation assumes surface pressure field
     !         estimate height of 1013 hPa using lowest upper level temp
     ZSFC=LOG(p1013/P0)*RDRY*T(1)/GRAV
     WRITE(KFJCLMSG,*) "prfwrf: SFLG=FALSE, computed zsfc from P0,T(1):",zsfc,P0,T(1)
  END IF


  !=>compute the log of the surface pressure

  IF(PFLG)THEN
     !        surface pressure from input data
     PBOT=LOG(P0)
  ELSE
     !        use lowest upper level temp to estimate surface pressure
     PBOT=LOG(p1013)-ZSFC*GRAV/RDRY/T(1)
     P0=EXP(PBOT)
  END IF


  !=>use adiabatic profile to estimate surface temperature

  IF(TFLG)THEN
     !        available use low level value
     TBOT=T0
  ELSE
     !        drop adiabatically to surface
     TBOT=T(1)*(PHALF(1)/pref)**(-0.286)
  END IF
  !     convert all temperatures to virtual
  IF(QFLG)TBOT=TBOT*(1.0+0.61*Q(1))

  !=>convert specific humidity to fraction of RH

  IF(QFLG)THEN
     ESAT=EXP(21.4-(5351.0/T(1)))
     RBOT=Q(1)*P0/(0.622*ESAT)
  ELSE
     RBOT=Q(1)/100.0
  END IF

  !=>set bottom dry inverse density
  if (fluxflg) altbot=alt1(1)

  !     starting height (agl) and vertical velocity
  ZBOT=0.0
  if (deepflg .or. shallflg) TKEN_bot=0.

  !=>integrate hypsometric equation from near-ground
  !     initial output vertical index
  KL=1
  !     initial output sigma level
  ! CHG(09/10/03) correct transformation between sigma and agl
  !      ZLVL=ZMDL*(1.0-ZSG(KL))
  !     adjust for terrain compression
  !      ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
  if (vcoord .eq. 1) then
     ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC) !use layer height AGL as the vertical coordinate for looping/interpolation
     vclvl=zlvl
     vcbot=ZBOT
  elseif (vcoord .eq. 2) then
     if (npar_wrfhyb .gt. 0) then
      write(KFJCLMSG,*) 'prfwrf: unsupported vcoord=2 for hybrid vertical coordinate WRF fields'
      write(*,*) 'prfwrf: unsupported vcoord=2 for hybrid vertical coordinate WRF fields'
      stop 'prfwrf: unsupported vcoord=2 for hybrid vertical coordinate WRF fields'
     end if
     vclvl=sign_eta*zsg(kl) ! use layer wrf eta (sigmap) as the vertical coordinate for looping/interpolation
     ! compute the needed height agl here:
     ZLVL=mu*(1.-abs(vclvl))*alt1(1)/grav
     vcbot=sign_eta*1.
  else
     write(KFJCLMSG,*) 'prfwrf: invalid vcoord=',vcoord
     write(*,*) 'prfwrf: invalid vcoord=',vcoord
     stop 'prfwrf: invalid vcoord'
  end if

  if (lprint .and. iprint .le. maxprint) then
     WRITE(KFJCLMSG,'(a/(a,g15.6))') 'Input and derived sfc variables:','zsfc=',zsfc, &
          & 'p0 =',p0,'t0=',t0,'tv (tbot)=',tbot, &
          & 'rh (rbot)=',rbot,'u0=',u0,'v0=',v0,'zbot=',zbot,'w0=',w0,'zmdl=',zmdl, &
          & 'mu (Pa)=',mu,'msfu=',msfu,'msfv=',msfv,'msft=',msft
     if (fluxflg) WRITE(KFJCLMSG,'(a,g15.6)') 'muu=',muu, &
          & 'muv=',muv
     WRITE(KFJCLMSG,'(a/a5,10a15)') 'Input profile:', &
          & 'k','phalf','zhalf','thet','temp','tv','q','rh','u','v','w_full'
  endif

  KZ_LOOP: DO KZ=1,NZ

     !        log of pressure at level
     PTOP=LOG(PHALF(KZ))
     !        convert input theta to temperature:
     temp=t(kz) * (phalf(kz)/pref) ** 0.286
     !        virtual temperature
     IF(QFLG)THEN
        TTOP=(1.0+0.61*Q(KZ))*Temp
     ELSE
        TTOP=Temp
     END IF

     if (vcoord .eq. 1) then
        vctop=half_agl(kz) !use layer height AGL as the vertical coordinate for looping/interpolation
     elseif (vcoord .eq. 2) then
        vctop=vcbot + sign_eta*deta_half(kz) ! use layer wrf eta (sigmap) as the vertical coordinate for looping/interpolation
        if (kz .eq. 1) pdh_top=phalf(kz)-mu*sign_eta*vctop/p2jm ! approximate WRF PTOP (in mb)
     end if

     !        only for the first input data level
     IF(KZ.EQ.1)THEN
        !           set low level (ground) wind data if available
        IF(UFLG)THEN
           UBOT=U0
           VBOT=V0
        ELSE
           !              use neutral log-law when output below data level
           IF(vcLVL.LT.vcTOP)THEN

              ! JCL:            occasionally, Z0=0, and model crashes b/c dividing by 0
              !                 so reset Z0 to a very small value (0.01)

              z0l=max(Z0,0.01)

              ATOP=LOG(half_agl(kz)/Z0l)
              ABOT=LOG(ZLVL/Z0l)

              UBOT=U(1)*ABOT/ATOP
              VBOT=V(1)*ABOT/ATOP

           ELSE
              UBOT=U(1)
              VBOT=V(1)
           END IF
           if (fluxflg) then !decouple u,v-winds
              if (npar_wrfhyb .gt. 0) then
                 ! replace mu(half-level) by (c1h*mu+c2h):
                 ubot=ubot*msfu/(wrfvcoords(1,1)*muu+wrfvcoords(1,2))
                 vbot=vbot*msfv/(wrfvcoords(1,1)*muv+wrfvcoords(1,2))
              else
                 ubot=ubot*msfu/muu
                 vbot=vbot*msfv/muv
              end if
           end if
        END IF
     END IF  !kz.eq.1

     !        convert to rh fraction
     IF(QFLG)THEN
        ESAT=EXP(21.4-(5351.0/temp))
        RTOP=Q(KZ)*PHALF(KZ)/(0.622*ESAT)
     ELSE
        RTOP=Q(KZ)/100.0
     END IF

     UTOP=U(KZ)
     VTOP=V(KZ)

     if (fluxflg) then
        !decouple u,v-winds
        if (npar_wrfhyb .gt. 0) then
           ! replace mu(half-level) by (c1h*mu+c2h):
           utop=utop*msfu/(wrfvcoords(kz,1)*muu+wrfvcoords(kz,2))
           vtop=vtop*msfv/(wrfvcoords(kz,1)*muv+wrfvcoords(kz,2))
        else
           utop=utop*msfu/muu
           vtop=vtop*msfv/muv
        end if
        !=>set top dry inverse density
        alttop=alt1(KZ)
     end if

     if (deepflg .or. shallflg) TKEN_top=TKEN_1(kz)

     if (lprint .and. iprint .le. maxprint) &
          & WRITE(KFJCLMSG,'(i5,10g15.6)') &
          & kz,phalf(kz),half_agl(kz),t(kz),temp,ttop,q(kz),rtop,utop,vtop,w(kz)

     KL_LOOP: DO WHILE (vcLVL.LE.vcTOP)
        ZZ(KL)=ZLVL

        !           basic linear interpolation
        FRAC=(vcLVL-vcBOT)/(vcTOP-vcBOT)
        TT(KL)=FRAC*(TTOP-TBOT)+TBOT
        RH(KL)=FRAC*(RTOP-RBOT)+RBOT
        UU(KL)=FRAC*(UTOP-UBOT)+UBOT
        VV(KL)=FRAC*(VTOP-VBOT)+VBOT
        alt2(kl)=FRAC*(altTOP-altBOT)+altBOT
        if (deepflg .or. shallflg) TKEN_2(kl)=FRAC*(TKEN_top-TKEN_bot) + TKEN_bot

        !           linear interpolation of log of pressure
        PP(KL)=EXP(FRAC*(PTOP-PBOT)+PBOT)

        !           density and virtual potential temperature from local value
        DEN(KL)=P2JM*PP(KL)/(TT(KL)*RDRY)
        TT(KL)=TT(KL)*(pref/PP(KL))**0.286

        KL=KL+1
        IF(KL.GT.NL) exit KZ_LOOP
        ! CHG(09/10/03) correct transformation between sigma and agl
        !            ZLVL=ZMDL*(1.0-ZSG(KL))
        !           adjust for terrain compression
        !            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
        if (vcoord .eq. 1) then
           ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC) !use layer height AGL as the vertical coordinate for looping/interpolation
           vclvl=ZLVL
        else
           vclvl=sign_eta*zsg(kl) ! use layer wrf eta (sigmap) as the vertical coordinate for looping/interpolation
           zlvl=zz(kl-1) - mu*(zsg(kl)-zsg(kl-1))*alt2(kl-1)/grav
        end if
     end DO KL_LOOP

     PBOT=PTOP
     vcBOT=vcTOP
     RBOT=RTOP
     TBOT=TTOP
     UBOT=UTOP
     VBOT=VTOP
     altbot = alttop

  end DO KZ_LOOP

  ! Separate loop for interpolation to w-levels
  maxdeep = 0.
  maxshall = 0.
  zbot=0.0
  wbot=w0
  if (vcoord .eq. 2) then
     wbot=0.
     ! Compute output full-level heights AGL
     hw(1)=0.
     call eta2agl(1.,zbot,mu,alt2,zsg,nl,zsg_full,dzsg_full,hm,hw(2:nl+1),.TRUE.)
     zbot=0.
  end if
  
  kzalt=1
  klw=1 !Initialize output full level klw=1: at surface
  zlvl=hw(klw) !=0
  if (deepflg) then
     CFXUP1_bot=0.
     CFXDN1_bot=0.
  endif
  if (shallflg) CFXUP2_bot=0.
  ! Loop over input full levels: kzw
  KZW_LOOP: do kzw=1,nz
     !      use previously computed full_agl heights:
     ZTOP=full_agl(kzw)
     if (vcoord .eq. 1) then
        vcbot=zbot !use full level height AGL as the vertical coordinate for looping/interpolation
        vctop=ztop
        vclvl=zlvl
     elseif (vcoord .eq. 2) then
        vcbot=sign_eta*1. ! use full level wrf eta (sigmap) as the vertical coordinate for looping/interpolation
        if (kzw .gt. 1) vcbot=sign_eta*eta(kzw-1)
        vctop=sign_eta*eta(kzw)
        if (klw .eq. 1) then
           vclvl=sign_eta*1.
        elseif (klw .gt. 1 .and. klw .le. nl) then
           vclvl=sign_eta*0.5*(zsg(klw-1)+zsg(klw))
        elseif (klw .eq. nl+1) then
           vclvl=sign_eta*(zsg(nl)+0.5*(zsg(nl)-zsg(nl-1)))
        else
           stop 'prfwrf - klw internal logic error 1'
        end if
        
     end if
     !        interpolate w from full levels to the w-levels
     !        (full-level variables are stored at sfc (for k=1) and at k-1 for k=2,n_full)
     wtop = w(kzw)
     if (deepflg) then
        if (kzw .lt. nz) then
           CFXUP1_top=CFXUP1_1(kzw)
           CFXDN1_top=CFXDN1_1(kzw)
        else
           CFXUP1_top=0.
           CFXDN1_top=0.
        endif
        maxdeep=max(maxdeep,abs(CFXUP1_top))
        maxdeep=max(maxdeep,abs(CFXDN1_top))
     endif
     if (shallflg) then
        if (kzw .lt. nz) then
           CFXUP2_top=CFXUP2_1(kzw)
        else
           CFXUP2_top=0.
        endif
        maxshall=max(maxshall,abs(CFXUP2_top))
     endif
     KLW_LOOP: DO WHILE (vcLVL.LE.vcTOP)
        FRAC=(vcLVL-vcBOT)/(vcTOP-vcBOT)
        !           w-interpolation:
        !           interpolate alt1 from half levels to w-level:
        do while (zlvl .ge. half_agl(kzalt) .and. kzalt .lt. nz)
           kzalt=kzalt+1
        end do
        if (kzalt .le. 1) then
           frac_alt=0.
        else
           !dwen(20090825)               frac_alt=min(dble(1),(zlvl-half_agl(kzalt-1))/(half_agl(kzalt)-half_agl(kzalt-1)))
           frac_alt=min(1.0,(zlvl-half_agl(kzalt-1))/(half_agl(kzalt)-half_agl(kzalt-1)))
        endif
        alt_int=alt1(max(1,kzalt-1)) + frac_alt*(alt1(kzalt)-alt1(max(1,kzalt-1)))

        !           interpolate w from full levels:            
        dzdt=FRAC*(WTOP-WBOT)+WBOT
        if (klw .le. nl) then
           if (.not. fluxflg) then
              !              vertical velocity (dz/dt) term converted to sigma/time:
              !              ww=d(zsg)/dt; zsg = (Zmdl-Z)/(Zmdl-Zsfc) [note sign change!]
              WW(KLw) = -dzdt/(ZMDL-ZSFC)
           else
              if (vcoord .eq. 1) then
                 !              WRF vertical flux term is [coupled omega + (eta/m)d(mu)/dt]
                 !              multiply by (-alt1/grav) to obtain dzdt (still divided by map scale factor m)
                 !              divide by -(ZMDL-ZSFC) and multiply by map scale factor m to convert to ww
                 WW(KLw) = dzdt * msft * (alt_int/grav) / (ZMDL-ZSFC)
              else
                 !              compute d(eta)/dt=(m/mu)*omega
                 ww(klw) = (msft/mu)*dzdt
              end if
           endif
        endif
        !
        ! NOTE on DFXUP1_2 etc ([CFD]FX(UP|DN)1_2) initialization:
        !     These are intent(out), and thus initialized in this routine.  
        !
        !     The assignment statements below (BLOCK A) make use of previously
        !     initialized values, even though they appear before the
        !     initialization assignment statements (BLOCK C).  This is because
        !     the array index klw is incremented (BLOCK B) right before BLOCK C.
        !
        !     Illustration of program flow for the first two KLW_LOOP iterations:
        !      first iteration of KLW_LOOP:
        !         klw=1: BLOCK A skipped; 
        !                BLOCK B: klw incremented (=2)
        !         klw=2: BLOCK C executed, DFXUP1_2(1) initialized
        !      second iteration of KLW_LOOP:
        !         klw=2: BLOCK A executed (DFXUP1_2(1) updated); 
        !                BLOCK B: klw incremented (=3)
        !         klw=3: BLOCK C executed, DFXUP1_2(2) initialized
        !
        ! ------Beginning of BLOCK A----------------------------------------------------
        if (klw .gt. 1) then
           frac_cnv = (zlvl-max(zbot,hw(klw-1)))/(ztop-zbot)
           if (deepflg) then
              ! straightforward interpolation of up/downdraft mass fluxes, but note that staggering follows
              ! cgrell convention: hw(klw=2) is first flux (sigw) level above ground, fluxes stored in CFXUP1_2(1)
              CFXUP1_2(klw-1)=FRAC*(CFXUP1_top-CFXUP1_bot) + CFXUP1_bot
              CFXDN1_2(klw-1)=FRAC*(CFXDN1_top-CFXDN1_bot) + CFXDN1_bot
              ! en/de-trainment fluxes are summed over layer:
              !   input(wrf levels):    dfxup1_1(kzw)   is detrainment between zbot and ztop; 
              !   output(STILT levels): dfxup1_2(klw-1) is detrainment between hw(klw-1)     and zlvl=hw(klw)
              !                  klw=2: dfxup1_2(1)     is detrainment between sfc [hw(1)=0] and zlvl=hw(2)
              DFXUP1_2(klw-1)=DFXUP1_2(klw-1)+frac_cnv*DFXUP1_1(kzw)
              DFXDN1_2(klw-1)=DFXDN1_2(klw-1)+frac_cnv*DFXDN1_1(kzw)
              EFXUP1_2(klw-1)=EFXUP1_2(klw-1)+frac_cnv*EFXUP1_1(kzw)
              EFXDN1_2(klw-1)=EFXDN1_2(klw-1)+frac_cnv*EFXDN1_1(kzw)
           endif
           if (shallflg) then
              CFXUP2_2(klw-1)=FRAC*(CFXUP2_top-CFXUP2_bot) + CFXUP2_bot
              DFXUP2_2(klw-1)=DFXUP2_2(klw-1)+frac_cnv*DFXUP2_1(kzw)
              EFXUP2_2(klw-1)=EFXUP2_2(klw-1)+frac_cnv*EFXUP2_1(kzw)
           endif
        endif !klw > 1
        ! ------End of BLOCK A----------------------------------------------------------
        ! ------Beginning of BLOCK B----------------------------------------------------
        KLw=KLw+1
        ! ------End of BLOCK B----------------------------------------------------------
        ! ------Beginning of BLOCK C----------------------------------------------------
        ! Initialize de/en-trainment fluxes
        if (klw .le. nl+1) then
           if (deepflg) then
              DFXUP1_2(klw-1)=0.
              DFXDN1_2(klw-1)=0.
              EFXUP1_2(klw-1)=0.
              EFXDN1_2(klw-1)=0.
           endif
           if (shallflg) then
              DFXUP2_2(klw-1)=0.
              EFXUP2_2(klw-1)=0.
           endif
           ! ------End of BLOCK C----------------------------------------------------------
        else
           exit KZw_LOOP
        endif
        zlvl = hw(klw) ! set zlvl=output full-level AGL height
        if (vcoord .eq. 1) then
           vclvl=zlvl !use full level height AGL as the vertical coordinate for looping/interpolation
        elseif (vcoord .eq. 2) then
           ! use full level wrf eta (sigmap) as the vertical coordinate for looping/interpolation
           if (klw .le. nl) then
              vclvl=sign_eta*0.5*(zsg(klw-1)+zsg(klw))
           elseif (klw .eq. nl+1) then
              vclvl=sign_eta*(zsg(nl)+0.5*(zsg(nl)-zsg(nl-1)))
           else
              stop 'prfwrf - klw internal logic error 2'
           end if
        end if
        
        if (zlvl .gt. ztop) then
           ! klw layer extends over old and new kzw layer: initialize with contribution from previous kzw layer
           frac_cnv = (ztop-hw(klw-1))/(ztop-zbot)
           if (deepflg) then
              DFXUP1_2(klw-1)=frac_cnv*DFXUP1_1(kzw)
              DFXDN1_2(klw-1)=frac_cnv*DFXDN1_1(kzw)
              EFXUP1_2(klw-1)=frac_cnv*EFXUP1_1(kzw)
              EFXDN1_2(klw-1)=frac_cnv*EFXDN1_1(kzw)
           endif
           if (shallflg) then
              DFXUP2_2(klw-1)=frac_cnv*DFXUP2_1(kzw)
              EFXUP2_2(klw-1)=frac_cnv*EFXUP2_1(kzw)
           endif
        endif !zlvl > ztop
     end DO KLW_LOOP
     wbot=wtop
     zbot=ztop
     vcBOT=vcTOP
     if (deepflg) then
        CFXUP1_bot=CFXUP1_top
        CFXDN1_bot=CFXDN1_top
     endif
     if (shallflg) CFXUP2_bot=CFXUP2_top
  end do KZW_LOOP

  !     sounding ends but levels remain to be filled (separately check klw)
  IF(KL .LE. NL .or. klw .le. nl+1)THEN
     IF(KFLAG.EQ.0)THEN
        WRITE(*,*)'WARNING prfwrf: data extrapolation because kl or klw < nl(+1): ', &
             &  kl, klw, nl
        ! CHG(09/10/03) correct transformation between sigma and agl
        !            ZLVL=ZMDL*(1.0-ZSG(KL))
        !           adjust for terrain compression
        !            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
        !            ZLVL=(1.0-ZSG(KL))*(ZMDL-ZSFC)
        if (kl .le. nl .and. vcoord .eq. 1) WRITE(*,*)'from level (m): ', &
             & (1.0-ZSG(KL))*(ZMDL-ZSFC)
        if (kl .le. nl .and. vcoord .eq. 2) WRITE(*,*)'from level (eta): ', &
             & ZSG(KL)
        KFLAG=1
     END IF

     DO KK=KL,NL
        ! CHG(09/10/03) correct transformation between sigma and agl
        !            ZLVL=ZMDL*(1.0-ZSG(KK))
        !           adjust for terrain compression
        !            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
        if (vcoord .eq. 1) then
           ZLVL=(1.0-ZSG(KK))*(ZMDL-ZSFC)
           DELZ=ZLVL-ZZ(KK-1)
           !           use previous pressure and density to find pressure
           PP(KK)=PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM
        elseif (vcoord .eq. 2) then
           PP(KK)=mu*zsg(kk)/p2jm + pdh_top !convert eta to pressure (approximate)
           delz=(pp(kk-1)-pp(kk))*p2jm/(den(kk-1)*GRAV) !invert formula above
           zlvl=zz(kk-1)+delz
        endif
        ZZ(KK)=ZLVL
        TT(KK)=TT(KK-1)
        RH(KK)=RTOP

        !           convert potential back to local temperature
        TEMP=TT(KK)*(PP(KK)/1000.0)**0.286

        DEN(KK)=P2JM*PP(KK)/(TEMP*RDRY)
        UU(KK)=UU(KK-1)
        VV(KK)=VV(KK-1)
        if (deepflg .or. shallflg) TKEN_2(kk)=TKEN_2(kk-1)

     END DO
     DO KK=KLw,NL+1
        ! CHG(09/10/03) correct transformation between sigma and agl
        !            ZLVL=ZMDL*(1.0-ZSG(KK))
        !           adjust for terrain compression
        !            ZLVL=ZMDL*ZLVL/(ZMDL-ZSFC)
        !           diminish magnitude when no data
        if (kk .le. NL) WW(KK)=WW(KK-1)/2.0
        if (kk .gt. 2) then
           if (deepflg) then
              CFXUP1_2(kk-1)=CFXUP1_2(kk-2)/2.0
              CFXDN1_2(kk-1)=CFXDN1_2(kk-2)/2.0
              DFXUP1_2(kk-1)=DFXUP1_2(kk-2)/2.0
              DFXDN1_2(kk-1)=DFXDN1_2(kk-2)/2.0
              EFXUP1_2(kk-1)=EFXUP1_2(kk-2)/2.0
              EFXDN1_2(kk-1)=EFXDN1_2(kk-2)/2.0
           endif
           if (shallflg) then
              CFXUP2_2(kk-1)=CFXUP2_2(kk-2)/2.0
              DFXUP2_2(kk-1)=DFXUP2_2(kk-2)/2.0
              EFXUP2_2(kk-1)=EFXUP2_2(kk-2)/2.0
           endif
        endif !kk > 2
     END DO
  END IF

  if (lprint .and. iprint .le. maxprint) then
     WRITE(KFJCLMSG,'(a/a5,9a15)') 'Output profile:', &
          & 'k','pp','zsg','zz','thetv','rh','u','v','w','den'
     do kk=1,nl
        WRITE(KFJCLMSG,'(i5,9g15.6)') &
             & kk,pp(kk),zsg(kk),zz(kk),tt(kk),rh(kk),uu(kk),vv(kk),ww(kk),den(kk)
     enddo
  endif

  if (maxdeep .lt. epsmax) then
     do kk=1,nl
        CFXUP1_2(kk)=0.0
        CFXDN1_2(kk)=0.0
        DFXUP1_2(kk)=0.0
        DFXDN1_2(kk)=0.0
        EFXUP1_2(kk)=0.0
        EFXDN1_2(kk)=0.0
     enddo
  endif
  if (maxshall .lt. epsmax) then
     do kk=1,nl
        CFXUP2_2(kk)=0.0
        DFXUP2_2(kk)=0.0
        EFXUP2_2(kk)=0.0
     enddo
  endif
  if (maxshall .ge. epsmax .or. maxdeep .ge. epsmax) then
     iprcnv = iprcnv + 1
     ! Check mass conservation of convective fluxes
     ! Input:
     if (iprcnv .le. maxprint) then
        WRITE(KFJCLMSG,'(a/a5,11a15)') 'Input convective fluxes:', &
             & 'kz','cfu1','cfd1','dfu1','dfd1','efu1','efd1','totmas1', &
             &      'cfu2','dfu2','efu2','totmas2'
        do kz=2,nz
           if (maxdeep .ge. epsmax) then
              totmas1=cfxup1_1(kz)-abs(cfxdn1_1(kz)) & !subin
                   & - (cfxup1_1(kz-1)-abs(cfxdn1_1(kz-1))) & !subdown
                   & + dfxup1_1(kz) + dfxdn1_1(kz) & !detrainment (detupd+_detupk+detdo)
                   & + efxup1_1(kz) + efxdn1_1(kz)   !entrainment (entup, entupk, entdo, edtdoj)
              WRITE(KFJCLMSG,'(i5,7g15.6)') &
                   & kz,cfxup1_1(kz),cfxdn1_1(kz),dfxup1_1(kz),dfxdn1_1(kz), &
                   & efxup1_1(kz),efxdn1_1(kz),totmas1
           endif
           if (maxshall .ge. epsmax) then
              totmas2=cfxup2_1(kz) & !subin
                   & - cfxup2_1(kz-1) & !subdown
                   & + dfxup2_1(kz) & !detrainment
                   & + efxup2_1(kz)   !entrainment
              WRITE(KFJCLMSG,'(i5,7(15x),4g15.6)') &
                   & kz,cfxup2_1(kz),dfxup2_1(kz),efxup2_1(kz),totmas2
           endif
        enddo
     endif
     if (iprcnv .le. maxprint) then
        WRITE(KFJCLMSG,'(a/a5,11a15)') 'Output convective fluxes:', &
             & 'kl','cfu1','cfd1','dfu1','dfd1','efu1','efd1','totmas1', &
             &      'cfu2','dfu2','efu2','totmas2'
     endif
     do kl=2,nl
        if (maxdeep .ge. epsmax) then
           totmas1=cfxup1_2(kl)-abs(cfxdn1_2(kl)) & !subin
                & - (cfxup1_2(kl-1)-abs(cfxdn1_2(kl-1))) & !subdown
                & + dfxup1_2(kl) + dfxdn1_2(kl) & !detrainment (detupd+_detupk+detdo)
                & + efxup1_2(kl) + efxdn1_2(kl)   !entrainment (entup, entupk, entdo, edtdoj)
           if (iprcnv .le. maxprint) WRITE(KFJCLMSG,'(i5,7g15.6)') &
                & kl,cfxup1_2(kl),cfxdn1_2(kl),dfxup1_2(kl),dfxdn1_2(kl), &
                & efxup1_2(kl),efxdn1_2(kl),totmas1
           if (abs(totmas1) .gt. epstotm*maxdeep) then
              ! adjust de/entrainment to agree with convective mass-flux profile:               
              if (cfxup1_2(kl) .gt. cfxup1_2(kl-1)) then
                 efxup1_2(kl)=cfxup1_2(kl-1) - cfxup1_2(kl)
                 dfxup1_2(kl)=0.
              else
                 dfxup1_2(kl)=cfxup1_2(kl-1) - cfxup1_2(kl)
                 efxup1_2(kl)=0.
              endif
              if (abs(cfxdn1_2(kl-1)) .gt. abs(cfxdn1_2(kl))) then
                 efxdn1_2(kl)=abs(cfxdn1_2(kl)) - abs(cfxdn1_2(kl-1))
                 dfxdn1_2(kl)=0.
              else
                 dfxdn1_2(kl)=abs(cfxdn1_2(kl)) - abs(cfxdn1_2(kl-1))
                 efxdn1_2(kl)=0.
              endif
              totmas1=cfxup1_2(kl)-abs(cfxdn1_2(kl)) & !subin
                   & - (cfxup1_2(kl-1)-abs(cfxdn1_2(kl-1))) & !subdown
                   & + dfxup1_2(kl) + dfxdn1_2(kl) & !detrainment (detupd+_detupk+detdo)
                   & + efxup1_2(kl) + efxdn1_2(kl)   !entrainment (entup, entupk, entdo, edtdoj)
              if (iprcnv .le. maxprint) WRITE(KFJCLMSG,'(5x,7g15.6)') &
                   & cfxup1_2(kl),cfxdn1_2(kl),dfxup1_2(kl),dfxdn1_2(kl), &
                   & efxup1_2(kl),efxdn1_2(kl),totmas1
           endif
        endif
        if (maxshall .ge. epsmax) then
           totmas2=cfxup2_2(kl) & !subin
                & - cfxup2_2(kl-1) & !subdown
                & + dfxup2_2(kl) & !detrainment
                & + efxup2_2(kl)   !entrainment
           if (iprcnv .le. maxprint) WRITE(KFJCLMSG,'(i5,7(15x),4g15.6)') &
                & kl,cfxup2_2(kl),dfxup2_2(kl),efxup2_2(kl),totmas2
           if (abs(totmas2) .gt. epstotm*maxshall) then
              ! adjust de/entrainment to agree with convective mass-flux profile:               
              if (cfxup2_2(kl) .gt. cfxup2_2(kl-1)) then
                 efxup2_2(kl)=cfxup2_2(kl-1) - cfxup2_2(kl)
                 dfxup2_2(kl)=0.
              elseif (cfxup2_2(kl) .lt. cfxup2_2(kl-1)) then
                 dfxup2_2(kl)=cfxup2_2(kl)-cfxup2_2(kl-1)
                 efxup2_2(kl)=0.
              endif
              totmas2=cfxup2_2(kl) & !subin
                   & - cfxup2_2(kl-1) & !subdown
                   & + dfxup2_2(kl) & !detrainment
                   & + efxup2_2(kl)   !entrainment
              if (iprcnv .le. maxprint) WRITE(KFJCLMSG,'(5x,7(15x),4g15.6)') &
                   & cfxup2_2(kl),dfxup2_2(kl),efxup2_2(kl),totmas2
           endif
        endif
     enddo
  endif

end SUBROUTINE PRFwrf
    
   
