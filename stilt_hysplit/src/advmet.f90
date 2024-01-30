!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVMET           ADVection step returns local METeorology
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION STEP RETURNS LOCAL METEOROLOGY AT END POINT AFTER
!   INTEPOLATION OF METEOROLOGICAL INFORMATION FROM GRID TO POSITION
!   IN BOTH SPACE AND TIME.  NOT USED IN ADVECTION BUT BY OTHER
!   ROUTINES AS DISPERSION AND DEPOSITION.  THIS ROUTINE PROVIDES TH
!   ONLY INTERFACE OF METEOROLOGICAL INFORMATION TO NON-ADVECTION
!   SUBROUTINES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 18 Aug 1998 (RRD) - added tension stress to deformation
!                 21 Dec 1998 (RRD) - generalized vertical index and
!                                     simplified non-interpolation option
!                                     added rdep to argument list
!                                     compile option for spatial interpolation
!                 04 Mar 1999 (RRD) - eliminate negative precip
!                 06 Apr 1999 (RRD) - created pressure profile, pass to METO
!                 30 Apr 1999 (RRD) - added terrain and height computation
!
! USAGE:  CALL ADVMET(BACK,VMIX,CDEP,RDEP,TRAJ,XP,YP,JET,DM,KVEL,
!              KCYCLE,NXS,NYS,NLVL,FHOUR,IFHR,K1,K2,GD,Z0,LU,ZT,
!              T1,Q1,P1,D1,X1,H1,RT1,UF1,VF1,
!              T2,Q2,P2,D2,X2,H2,RT2,UF2,VF2)
!   INPUT ARGUMENT LIST:
!     BACK  - log backward direction integration flag
!     VMIX  - log flag to return mixing profile
!     CDEP  - log flag to return deposition variables
!     RDEP  - log flag to indicate resistance deposition
!     TRAJ  - log flag to return trajectory variables
!     XP,YP - real      horizontal particle position
!     JET   - int elapsed time
!     DM    - real      minutes between data time periods
!     KVEL  - int vertical velocity computation method
!     KCYCLE      - int precipitation accumulation cycle (min)
!     NXS,NYS   - int   horizontal grid dimensions
!     NLVL  - int number of vertical levels
!     FHOUR     - int   forecast hour at last and next times
!     IFHR      - int   forecast hour interpolated to current time
!     GD    - real      grid distance (m)
!     Z0    - real      aerodynamic roughness length (m)
!     LU    - int land-use category (1-11)
!     ZT    - real      terrain height elevation (m)
!     ... indicies 1 and 2 define data at last and next time period
!     T,Q,P,D     - real      input meteo temp, rh, press, density
!     X,H   - real      vertical and horizontal mixing
!     U,V   - real      input meteo u,v winds
!     RT    - real      input rainfall total accumulation (m)
!     UF    - real      friction velocity (m/s)
!     VF    - real      stability function
!   OUTPUT ARGUMENT LIST:
!     GBLMET COMMON WITH VARIABLES DEFINED IN DEFMETO%INC
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: advmet.f90,v 1.19 2009-03-05 08:26:21 gerbig Exp $
!
!$$$

! JCL:add two arrays of Lagrangian timescales (TL1,TL2) and
!         two arrays of stddev of vertical velocity (SIGW1,SIGW2)
! JCL:(5/9/01)added arrays of horizontal velocity (U1,U2) and (V1,V2)
! CHG:(11/20/01) added conv. precip. flag CFLG to arguments
! CHG:(11/20/01)added conv. precip. rates (RC1,RC2) to arguments
! CHG:(11/20/01) added tot.cld and radiation flag (TCLF,RADF)
! CHG:(11/20/01)added tot. cld and radiation (TC1,TC2,SW1,SW2)
! CHG:(11/20/01) added energy fluxes (LF1,LF2,HF1,HF2) to arguments
! CHG:(11/20/01)added low.cld flag LCLF
! CHG:(11/20/01) added low.cld LC1, LC2
! JCL:(4/3/02)added arrays of mass violation (DMASS1,DMASS2)
! CHG:(9/24/02) added TLAT&TLON to arguments, for DSWF
! CHG:(22/01/03)added soil moisture flag SLMF
! CHG:(22/01/03) added soil moisture SM1 SM2
! CHG(09/16/03)changed input position XP,YP to XPI,YPI
! CHG(09/16/03) added RAMSFLG
! JCL:(07/12/2004) added cyclic boundary condition flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
SUBROUTINE ADVMET (BACK,VMIX,CDEP,RDEP,TRAJ,XPI,YPI,TLON,TLAT,JET, &
         DM,KVEL,KCYCLE,NXS,NYS,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,   &
         CFLG,TCLF,LCLF,RADF,SLMF,TFLG,                                 &
         T1,Q1,P1,D1,X1,H1,RT1,UF1,VF1,U1,V1,DMASS1,                    &
         T2,Q2,P2,D2,X2,H2,RT2,UF2,VF2,U2,V2,DMASS2,                    &
         TL1,TL2,SIGW1,SIGW2,RC1,RC2,HF1,HF2,LF1,LF2,                   &
         TC1,TC2,LC1,LC2,SW1,SW2,SM1,SM2,T01,T02,                       &
         RAMSFLG,ECMFLG,GLOBAL,NXP,NYP, &
         awrfflg, fluxflg)

!     dimension information
!     contains meteorological summary at last advection point
      USE module_defmeto

      IMPLICIT REAL*8 (A-H,O-Z)


      LOGICAL GLOBAL
!     input meteorological data arrays at times 1 and 2
! CHG:(11/20/01) add conv. precip RC, tot.cld TC and radiation SW
! CHG:(11/20/01)added energy fluxes (LF1,LF2,HF1,HF2)
! CHG:(12/04/01) low.cld LC
! CHG:(22/01/03)added soil moisture SM1 SM2
      REAL*8 T1(NXS,NYS,NZM), Q1(NXS,NYS,NZM), D1(NXS,NYS,NZM),         &
           P1(NXS,NYS,NZM), X1(NXS,NYS,NZM), RT1(NXS,NYS),              &
           H1(NXS,NYS,NZM), UF1(NXS,NYS),    VF1(NXS,NYS),              &
           RC1(NXS,NYS),HF1(NXS,NYS),LF1(NXS,NYS),                      &
           TC1(NXS,NYS),LC1(NXS,NYS),SW1(NXS,NYS),                      &
           SM1(NXS,NYS),T01(NXS,NYS)

      REAL*8 T2(NXS,NYS,NZM), Q2(NXS,NYS,NZM), D2(NXS,NYS,NZM),         &
           P2(NXS,NYS,NZM), X2(NXS,NYS,NZM), RT2(NXS,NYS),              &
           H2(NXS,NYS,NZM), UF2(NXS,NYS),    VF2(NXS,NYS),              &
           RC2(NXS,NYS),HF2(NXS,NYS),LF2(NXS,NYS),                      &
           TC2(NXS,NYS),LC2(NXS,NYS),SW2(NXS,NYS),                      &
           SM2(NXS,NYS),T02(NXS,NYS)
! JCL:arrays of Lagrangian timescale
      REAL*8 TL1(NXS,NYS,NZM),TL2(NXS,NYS,NZM)
! JCL:arrays of stddev of vertical velocity
      REAL*8 SIGW1(NXS,NYS,NZM),SIGW2(NXS,NYS,NZM)
! JCL:(5/9/01)added arrays of horizontal velocity
      REAL*8 U1(NXS,NYS,NZM), U2(NXS,NYS,NZM),                          &
             V1(NXS,NYS,NZM), V2(NXS,NYS,NZM)
! JCL:(4/3/02)added arrays of mass violation
      REAL*8 DMASS1(NXS,NYS,NZM),DMASS2(NXS,NYS,NZM)
! JCL:(4/3/02)temp variables to hold density
      REAL*8 DENS1,DENS2

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!     grid spacing, roughness length, terrain elevation
!     REAL*8 GD(NXS,NYS), Z0(NXS,NYS), ZT(NXS,NYS)
      REAL*8 GX(NXS,NYS),GY(NXS,NYS),Z0(NXS,NYS), ZT(NXS,NYS)
!     land use category
      INTEGER LU(NXS,NYS)

!     forecast times
      INTEGER FHOUR(2)

!     time step met files, time at end of met interval
      INTEGER JETDEL, JETMAX

!     flags: direction, vertical mixing, deposition, trajectory
      LOGICAL BACK,     VMIX,            CDEP,RDEP,   TRAJ

! CHG(09/11/03):add flag specifying whether data from RAMS or not
      LOGICAL, INTENT(IN) :: RAMSFLG, ECMFLG
      LOGICAL awrfflg, fluxflg

! CHG:(11/20/01) flag conv. precip, tot. cld., radiation
! CHG:(12/04/01) flag low. cld.
! CHG:(22/01/03) flag soil moisture SLMF
      LOGICAL CFLG, TCLF, LCLF, RADF, SLMF, TFLG


!---------------------------------------------------------------------------------------------------
! CHG(09/16/03) use proper grid stagger, see comments in advpnt.f
      IF (RAMSFLG) THEN
        XP = DNINT(XPI)
        YP = DNINT(YPI)
      ELSE
        XP = XPI
        YP = YPI
        if (awrfflg) then
!analogous to treatment in advpnt:
           XX2 = Xpi+0.5
           YY2 = Ypi+0.5
        end if
      END IF
! JCL:
!      WRITE(*,*) 'IN ADVMET:',XP,YP
!     interpolation factor for current time
      TF = MOD(DBLE(JET),DM)/DM
      IF (BACK) THEN
         TF = 1.0-TF
      ELSE
         IF (TF == 0.0) TF = 1.0
      END IF

!     compute current forecast hour
      RFHR = (FHOUR(K2)-FHOUR(K1))*TF+FHOUR(K1)
      IFHR = NINT(RFHR)

!     save position indicies within meteo subgrid
      II = NINT(XP)
      JJ = NINT(YP)
      ZX = METO%ZNDX
      KK = NINT(ZX)

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!     grid distance at point
!     METO%GDIS=GD(II,JJ)
! JCL:(07/26/2004) impose global cyclic boundary condition; otherwise could exit subgrid
      IF (GLOBAL) THEN
         IF (II > NXP) II = 1
         IF (JJ > NYP) JJ = NYP
      END IF
      if (ii .gt. nxs .or. jj .gt. nys) then
         print *, 'advmet warning: resetting ii,jj: ',ii,jj,' to: ',nxs,nys
         ii = min(ii,nxs)
         jj = min(jj,nys)
      endif
      METO%GDISX = GX(II,JJ)
      METO%GDISY = GY(II,JJ)
!     terrain surface elevation
      METO%ZTER = ZT(II,JJ)

! CHG:(11/20/01) energy fluxes
!        interpolate sensible heat flux to point and time
         CALL ADVINT2D (HF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (HF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%SHTF = (VAR2-VAR1)*TF+VAR1
!        interpolate latent heat flux to point and time
         CALL ADVINT2D (LF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (LF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%LHTF = (VAR2-VAR1)*TF+VAR1

! near surface temperature (T02M) for WRF to drive VPRM
      IF (TFLG .and. awrfflg) THEN
!        interpolate T02M to point and time
         CALL ADVINT2D (T01,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (T02,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%T02M = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%T02M = -99.0
      END IF

! CHG:(11/20/01) total cloud cover and downw. radiation
      IF (TCLF) THEN
!        interpolate total cloud cover to point and time
         CALL ADVINT2D (TC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (TC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%TCLD = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%TCLD = -99.0
      END IF

! CHG:(22/01/03) soil moisture
      IF (SLMF) THEN
!        interpolate total cloud cover to point and time
         CALL ADVINT2D (SM1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (SM2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%SOLW = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%SOLW = -99.0
      END IF

! CHG:(12/04/01) low cloud cover
      IF (LCLF) THEN
!        interpolate low cloud cover to point and time
         CALL ADVINT2D (LC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (LC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%LCLD = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%LCLD = -99.0
      END IF

      IF (RADF) THEN
!        interpolate downw. rad. flux to point and time
         CALL ADVINT2D (SW1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (SW2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
! CHG:(9/24/02) better interpolate "cloudiness factor"
! i.e. DSWF/DSFW(clear)
! so first get clear sky DSWF
         IF (BACK) THEN
! CHG:(1/28/03) fixed bug in JET1 and JET2
!            JET1 = JET+(1-TF)*DM
!            JET2=JET-TF*DM
            JET1 = JET+TF*DM
            JET2 = JET-(1.0-TF)*DM
         ELSE
            JET1 = JET-TF*DM
            JET2 = JET+(1.0-TF)*DM
         END IF

!     interpolation factor for current time
!      TF = MOD(DBLE(JET),DM)/DM
!      IF(BACK)THEN
!         TF = 1.0-TF
!      ELSE
!         IF(TF == 0.0) TF = 1.0
!      END IF

         JETDEL=IABS(JET1-JET2)
         JETMAX=MAX0(JET1,JET2)
         CALL SUNAVE(JETMAX,JETDEL,TLAT,TLON,SWFC1)
         CALL SUNCLR (JET,TLAT,TLON,SWFC)
         !for other fields than ECMWF use instantaneous clear sky radiation,
         !since radiation is not time integrated over met. output time step
         if(.not. ecmflg) SWFC1=SWFC
! Ratio DSWF(cld) to DSWF(clear); VAR? is now "cloudiness factor"
         IF (VAR1 < 10d0*SWFC1 .AND. SWFC1 > 0d0) THEN
            VAR1 = VAR1/SWFC1
         ELSE
            VAR1 = -1d0
         END IF
         IF (VAR2 < 10d0*SWFC1 .AND. SWFC1 > 0d0) THEN
            VAR2=VAR2/SWFC1
         ELSE
            VAR2 = -1d0
         END IF
         IF(JET1 > JET2) THEN
            IF(VAR1 < 0d0) THEN
               METO%DSWF=0d0
            ELSE
               METO%DSWF=VAR1*SWFC
            ENDIF
         ELSE
            IF(VAR2 < 0d0) THEN
              METO%DSWF=0d0
            ELSE
               METO%DSWF=VAR2*SWFC
            ENDIF
         ENDIF
!         WRITE(45,*)"advmet JET 1 2,TF:",JET,JET1,JET2,TF,DM
!         WRITE(45,*) "advmet VAR1,VAR2:",VAR1,VAR2,METO%DSWF
!         WRITE(45,*)"advmet SWFC,SWFC1,SWFC2:",SWFC,SWFC1,SWFC2
!         WRITE(45,*) "TLAT,TLON:",TLAT,TLON
      ELSE
         METO%DSWF = -999.0
      END IF

!==>surface level parameters

      IF (RDEP) THEN
!        roughness length and land-use category at point
         METO%AERO = Z0(II,JJ)
         METO%LAND = LU(II,JJ)

!        interpolate friction velocity to point and time
         CALL ADVINT2D (UF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (UF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%USTR = (VAR2-VAR1)*TF+VAR1
!        optional statement to simplify previous code
!        METO%USTR = (UF2(II,JJ)-UF1(II,JJ))*TF+UF1(II,JJ)

!        interpolate stability function to point and time
         CALL ADVINT2D (VF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT2D (VF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         METO%PSI = (VAR2-VAR1)*TF+VAR1
!        optional statement to simplify previous code
!        METO%PSI = (VF2(II,JJ)-VF1(II,JJ))*TF+VF1(II,JJ)
      END IF

! CHG:(11/20/01) get rain anyway
!      IF(CDEP) THEN
!        interpolate 2D precipitation total (m) to position
         IF (KCYCLE > 0) THEN
            CALL ADVINT2D (RT1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT2D (RT2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
!           optional statement to simplify previous code
!           VAR2 = RT2(II,JJ)
!           VAR1=RT1(II,JJ)

!           check if previous time requires zero accumulation
!           it is an even multiple rain accumulation initialization
            IF (MOD(FHOUR(K1)*60,KCYCLE) == 0) VAR1 = 0.0
!           convert accumulation to precipitation rate (m/min)
            METO%RAIN = MAX(DBLE(0.0),(VAR2-VAR1)/DM)
         ELSE
!           accumlation cycle not set implies no precip field
            METO%RAIN = 0.0
         END IF
! CHG:(11/20/01) get rain anyway
!      END IF


      IF (CFLG) THEN
         !  interpolate 2D precipitation total (m) to position
         IF (KCYCLE > 0) THEN
            CALL ADVINT2D (RC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT2D (RC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
            ! optional statement to simplify previous code
            !VAR2 = RC2(II,JJ)
            !VAR1=RC1(II,JJ)

            ! check if previous time requires zero accumulation
            ! it is an even multiple rain accumulation initialization
            IF (MOD(FHOUR(K1)*60,KCYCLE) == 0) VAR1 = 0.0
!           convert accumulation to precipitation rate (m/min)
            METO%CRAI = MAX(DBLE(0.0),(VAR2-VAR1)/DM)
         ELSE IF (ECMFLG) THEN                  ! instantaneous fields (averaged of 3 hours before)
            ! interpolate convective rain (in m/s) to point and time
            CALL ADVINT2D (RC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT2D (RC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
            METO%CRAI = ((VAR2-VAR1)*TF+VAR1)*60d0 ! output is in m/min
         ELSE
!           accumlation cycle not set implies no precip field
            METO%CRAI = -0.9
         END IF
      ELSE
! CHG:(11/20/01) if no conv. precip avail., set to -999
         METO%CRAI = -0.9
      END IF


!==>parameter value at vertical position

      IF (VMIX) THEN
!        to obtain the horizontal diffusivity in m2/sec
! JCL:(4/26/02) since ZX could be less than 1.0, prevent ADVINT from accessing 0th element in array
         ZTMP = ZX
         IF (ZX < 1.0) ZTMP = 1.0
! CHG(09/16/03) use next higher level for RAMS vertical (see adviec comments CHG)
         IF (RAMSFLG) ZTMP = DINT(ZX+1.0)
         ZTMP = DMIN1(DMAX1(DBLE(1.0),ZTMP),DBLE(NLVL))
         CALL ADVINT (H1,NXS,NYS,NZM,XP,YP,ZTMP,GLOBAL,NXP,NYP,VAR1)
         CALL ADVINT (H2,NXS,NYS,NZM,XP,YP,ZTMP,GLOBAL,NXP,NYP,VAR2)
!        CALL ADVINT(H1,NXS,NYS,NZM,XP,YP,ZX,VAR1)
!        CALL ADVINT(H2,NXS,NYS,NZM,XP,YP,ZX,VAR2)
         METO%HMIX = (VAR2-VAR1)*TF+VAR1
!        optional statement to simplify previous code
!        METO%HMIX = (H2(II,JJ,KK)-H1(II,JJ,KK))*TF+H1(II,JJ,KK)

!==>parameter profiles used in deposition and mixing calculations

         DO KL=1,NLVL
!           set vertical interpolation point to index position
            ZK = KL

! JCL:      interpolate Lagrangian timescale to position
! JCL:(10/30/02) no longer conduct linear interpolation, since screws up turbulence profiles
! JCL:(10/31/02)move extraction of TL to ADVPNT
!            XX = DNINT(XP)
!            YY=DNINT(YP)
!            CALL ADVINT(TL1,NXS,NYS,NZM,XX,YY,ZK,VAR1)
!            CALL ADVINT(TL2,NXS,NYS,NZM,XX,YY,ZK,VAR2)
!            METO%TL(KL) = (VAR2-VAR1)*TF+VAR1

! JCL:      interpolate std dev of vertical velocity to position
! JCL:(10/30/02) no longer conduct linear interpolation, since screws up turbulence profiles
! JCL:(10/31/02)move extraction of SIGMAW to ADVPNT
!            XX = DNINT(XP)
!            YY=DNINT(YP)
!            CALL ADVINT(SIGW1,NXS,NYS,NZM,XX,YY,ZK,VAR1)
!            CALL ADVINT(SIGW2,NXS,NYS,NZM,XX,YY,ZK,VAR2)
!            METO%SIGW(KL) = (VAR2-VAR1)*TF+VAR1

! JCL:(5/9/01)interpolate horizontal winds to position
! CHG(09/16/03) don't use this for RAMS
            IF (.NOT.RAMSFLG) THEN
               if (.not. awrfflg) then
                  CALL ADVINT (U1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                  CALL ADVINT (U2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                  METO%UUNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                  CALL ADVINT (V1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                  CALL ADVINT (V2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                  METO%VVNEXT(KL) = (VAR2-VAR1)*TF+VAR1
               else
! for awrf:
                  if (.not. fluxflg) then
! instantantenous velocities, on staggered C-grid
                     CALL ADVINT (U1,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                     CALL ADVINT (U2,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                     METO%UUNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                     CALL ADVINT (V1,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR1)
                     CALL ADVINT (V2,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                     METO%VVNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                  else
! time-averaged coupled u,v (decoupled in prfcom): (on C-grid)
                     IF (BACK) THEN
                        CALL ADVINT (U1,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                        CALL ADVINT (V1,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                     else
                        CALL ADVINT (U2,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                        CALL ADVINT (V2,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                     end IF !back
                     METO%UUNEXT(KL) = VAR1
                     METO%VVNEXT(KL) = VAR2
                  end if !fluxflg
               end if !awrfflg
                 !of IF(.NOT.RAMSFLG)THEN
            ELSE
! CHG(09/16/03) use this for RAMS
! use later time for fluxes
              IF (BACK) THEN
                CALL ADVINT (U1,NXS,NYS,NZM,XPI,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                CALL ADVINT (V1,NXS,NYS,NZM,XP,YPI,ZK,GLOBAL,NXP,NYP,VAR2)
              ELSE
                CALL ADVINT (U2,NXS,NYS,NZM,XPI,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                CALL ADVINT (V2,NXS,NYS,NZM,XP,YPI,ZK,GLOBAL,NXP,NYP,VAR2)
              END IF
              METO%UUNEXT(KL) = VAR1
              METO%VVNEXT(KL) = VAR2
                   !of IF(.NOT.RAMSFLG)THEN ... ELSE ...
            END IF

!           interpolate mixing to position
            CALL ADVINT (X1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (X2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            METO%VMIX(KL) = (VAR2-VAR1)*TF+VAR1
!           optional statement to simplify previous code
!           METO%VMIX(KL) = (X2(II,JJ,KL)-X1(II,JJ,KL))*TF+X1(II,JJ,KL)

! JCL:(2/23/01) need density for correction in PARDSP, so should extract
!                 density profile even when CDEP = F
!           interpolate density to position
            CALL ADVINT (D1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,DENS1)
            CALL ADVINT (D2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,DENS2)
            METO%DENS(KL) = (DENS2-DENS1)*TF+DENS1

! JCL:(4/3/02)should NOT have spatial interpolation of mass violation b/c mass violation
!           was calculated for a gridcell, so should keep same value
! CHG(10/10/03) special for RAMS
            IF (.NOT.RAMSFLG) THEN
              XX = DINT(XP)
              YY = DINT(YP)
              ZZ = DINT(ZK)
                 !for RAMS
            ELSE
              XX = XP
              YY = YP
              ZZ = DBLE(ZK)
            ENDIF
            CALL ADVINT (DMASS1,NXS,NYS,NZM,XX,YY,ZZ,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (DMASS2,NXS,NYS,NZM,XX,YY,ZZ,GLOBAL,NXP,NYP,VAR2)
            IF (.NOT. (RAMSFLG .OR. fluxflg) .OR. ECMFLG) THEN
               METO%DMASSNEXT(KL) = (VAR2-VAR1)*TF+VAR1
            ELSE IF (BACK) THEN
                  METO%DMASSNEXT(KL) = VAR1
               ELSE
                  METO%DMASSNEXT(KL) = VAR2
            END IF
! JCL:(4/3/02)the density change term is extremely small, so not worry
!          also need the density change term (kg/m3/min) as part of mass budget!
            METO%DMASSNEXT(KL) = METO%DMASSNEXT(KL)*2.0/(DENS2+DENS1)     &
                        +((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
!        take into account the sign associated with direction in TIME
            IF (BACK) METO%DMASSNEXT(KL) = -1.0*METO%DMASSNEXT(KL)


! JCL:(2/28/01) want to output surface temperature, which uses pressure, so should
!                extract pressure profile even when CDEP = F
!              interpolate pressure to position (currently not saved)
            CALL ADVINT (P1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (P2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            METO%PRES(KL) = (VAR2-VAR1)*TF+VAR1

! JCL:(2/28/01) want to output surface temperature, so should extract
!                temperature profile even when CDEP = F
!           compute temperature from pressure and density
            METO%TEMP(KL) = 100.0*METO%PRES(KL)/METO%DENS(KL)/287.04

! JCL:(2/28/01) relative humidity profile needed to calculate incident
!                solar radiation in SUNFLX, so extract even when CDEP = F
!           interpolate humidity to position
            CALL ADVINT (Q1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (Q2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            METO%RHFR(KL) = (VAR2-VAR1)*TF+VAR1

            IF (CDEP) THEN
!              interpolate pressure to position (currently not saved)
               CALL ADVINT (P1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
               CALL ADVINT (P2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
               METO%PRES(KL) = (VAR2-VAR1)*TF+VAR1
!              optional statement to simplify previous code
!              METO%PRES(KL)=
!    :             (P2(II,JJ,KL)-P1(II,JJ,KL))*TF+P1(II,JJ,KL)

!              interpolate humidyty to position
               CALL ADVINT (Q1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
               CALL ADVINT (Q2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
               METO%RHFR(KL) = (VAR2-VAR1)*TF+VAR1
!              optional statement to simplify previous code
!              METO%RHFR(KL)=
!    :             (Q2(II,JJ,KL)-Q1(II,JJ,KL))*TF+Q1(II,JJ,KL)

!              interpolate density to position
               CALL ADVINT (D1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
               CALL ADVINT (D2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
               METO%DENS(KL) = (VAR2-VAR1)*TF+VAR1
!              optional statement to simplify previous code
!              METO%DENS(KL)=
!    :             (D2(II,JJ,KL)-D1(II,JJ,KL))*TF+D1(II,JJ,KL)

!              compute temperature from pressure and density
               METO%TEMP(KL) = 100.0*METO%PRES(KL)/METO%DENS(KL)/287.04
            END IF

         END DO
      END IF

!==>trajectory option only returns simple marker variable

      IF (TRAJ) THEN
!        set marker variable according to vertical motion option
         IF (KVEL == 0 .OR. KVEL == 1 .OR. KVEL == 3 .OR. KVEL == 4) THEN
!           data, isobaric, density, isosigma -> save pressure
            CALL ADVINT (P1,NXS,NYS,NZM,XP,YP,ZX,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (P2,NXS,NYS,NZM,XP,YP,ZX,GLOBAL,NXP,NYP,VAR2)

         ELSEIF(KVEL == 2) THEN
!           isentropic -> save potential temperature
            CALL ADVINT (T1,NXS,NYS,NZM,XP,YP,ZX,GLOBAL,NXP,NYP,VAR1)
            CALL ADVINT (T2,NXS,NYS,NZM,XP,YP,ZX,GLOBAL,NXP,NYP,VAR2)

         END IF
         METO%TMRK = (VAR2-VAR1)*TF+VAR1
      END IF

END SUBROUTINE ADVMET
