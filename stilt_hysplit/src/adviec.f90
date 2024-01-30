!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVIEC           ADVection via Improved Euler-Cauchy
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION VIA IMPROVED EULER-CAUCHY COMPUTES THE 3D ADVECTION OF
!   POINT IN SPACE AND TIME USING A TWO-STEP PROCESS.  THE ADVECTION
!   THE RESULTS OF AN AVERAGE VELOCITY AT THE INITIAL POINT IN SPACE
!   TIME AND THE VELOCITY AT THE FIRST GUESS POSITION (USING THE INITIAL
!   VELOCITY).  ALL INTERPOLATION IS LINEAR.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 26 Aug 1999 (RRD) - common block to define vertical grid
!
! USAGE:  CALL ADVIEC (U1,V1,W1,U2,V2,W2,NXS,NYS,NZM,NLVL,DM,JET,
!              ZMDL,X1,Y1,Z1,X2,Y2,Z2,ZX,DT,BACK)
!   INPUT ARGUMENT LIST:
!     U,V   - real      component winds in meters per second
!     W           - real      vertical wind in grid/minute
!     ...indicies 1&2 for data at the last and next time
!     NXS,NYS     - int horizontal subgrid dimensions
!     NZM   - int vertical grid dimension
!     NLVL  - int number of levels to process
!     DM    - real      minutes between winds 1 and 2
!     JET   - int current elapsed time (minutes)
!     ZMDL  - real      internal grid vertical model domain top (m)
!     X1,Y1,Z1    - real      old position at time t
!     DT    - real      integration step (minutes)
!     BACK      - log   logical flag to indicate direction
!   OUTPUT ARGUMENT LIST:
!     X2,Y2,Z2    - real      new position at time t+dt
!     ZX    - real      last estimate of vertical index
!     DEAD  - log   logical flag to indicate if particle off grid
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: adviec.f90,v 1.11 2006-11-23 14:27:00 gerbig Exp $
!
!$$$

! JCL:(8/13/01) add grd ht (ZTER) & 1st sigma level (ZSG1) as input to calculate AGL for accurate interpolation
! JCL:(3/1/01) add WWOUT used to output the Wbar [sigma/min]
! JCL: add DUDZ, DVDZ to argument list--will be OUTPUT--CHG needs it to parameterize box
! JCL: 'GD' is array of grid size (m)--needed to convert windspeed from (grid/min)=>(m/s)
! JCL:(09/01/03) pass on wind error flag, UPRIME/UERR and VPRIME/VERR from previous timestep
! JCL:(09/01/03) pass on random seed 'RSEED' for modeling transport error as stochastic process
! JCL:(09/01/03) DXERR and DYERR are the horizontal displacements resulting from transport error
! CHG(09/11/03) Pass on RAMSFLG
! CHG(09/11/03) Pass on density fields D1,D2
! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
      SUBROUTINE ADVIEC (U1,V1,W1,U2,V2,W2,D1,D2,NXS,NYS,NZM,NLVL,DM,JET, &
                         ZMDL,ZTER,ZSG1,X1,Y1,Z1,X2,Y2,Z2,ZX,DT,BACK,     &
                         GX,GY,DUDZ,DVDZ,WWOUT,awrfflg,fluxflg, zsg,      &
                         RAMSFLG,ECMFLG,GLOBAL,NXP,NYP,DEAD)
!                        WINDERRTF,UERRPREV,VERRPREV,
!                        RSEED,DXERR,DYERR)

      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BACK,GLOBAL,DEAD
      REAL*8   U1(NXS,NYS,NZM), V1(NXS,NYS,NZM), W1(NXS,NYS,NZM)
      REAL*8   U2(NXS,NYS,NZM), V2(NXS,NYS,NZM), W2(NXS,NYS,NZM)
! CHG(09/18/03) add density fields
      REAL*8   D1(NXS,NYS,NZM), D2(NXS,NYS,NZM)

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
! JCL:grid size array
!     REAL*8 GD(NXS,NYS)
      REAL*8 GX(NXS,NYS),GY(NXS,NYS)

! JCL:variables used to calculate wind shear
      REAL*8 ZLO, ZHI, ULO, VLO, UHI, VHI, DUDZ, DVDZ
! JCL:(3/1/01) variable used to output the Wbar (vertical velocity) [grid/min]
      REAL*8 WWOUT

! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! JCL:(09/01/03) WINDERRTF specifies whether to include wind errs as Markov process (=1) or not (=0)
!      INTEGER WINDERRTF
! JCL:(09/01/03) UPRIME/UERR and VPRIME/VERR from previous timestep
!      REAL*8 UERRPREV,VERRPREV
! JCL:(09/01/03) add random seed 'RSEED' for modeling transport error as stochastic process
!      INTEGER RSEED
! JCL:(09/01/03) stddev of wind errors [m/s] and its Lagrangian correlation timescale [min]
!      REAL*8 SIGUVERR,TLUVERR
!      LOGICAL FTEST

! CHG(09/16/03):add flag specifying whether data from RAMS or not
      LOGICAL, INTENT(IN) :: RAMSFLG, ECMFLG
      LOGICAL awrfflg, fluxflg

      real*8 zsg(nlvl)

!     common block to pass parameters for vertical grid
      COMMON /ZZTOKK/ AA,BB,CC

!---------------------------------------------------------------------------------------------------
!     set default to on grid
      dead=.FALSE.

!     need to save initial position in two-pass integration
      XT = X1
      YT = Y1
      ZT = Z1
!     determine interpolation time factor
      TF = DMOD(DBLE(JET),DM)/DM
      IF (BACK) THEN
         IF (TF == 0.0) TF = 1.0
         TF = 1.0-TF
      END IF

!     two pass to obtain average from velocity at old and new
      DO IPASS = 1,2

         !        compute vertical interpolation factor as fraction
         !        from integer index based upon quadratic relation between
         !        height (or sigma) and array index position

         !        height agl when zsfc = 0
         !        relation between index and sigma when zsfc = 0
         ! jcl:     ZX is an INDEX;see RUNSET: quadratic relationship defines K<=>ht levels=>sigmas
         ! CHG(09/11/03) don't use this with RAMS
         ! use the ind_zsg routine instead to allow for externally specified zsg levels:
         !!$         ZZ = ZMDL*(1.0-DMIN1(DBLE(1.0),ZT))
         !!$         IF(.NOT.RAMSFLG) ZX = (-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
         IF (.NOT.RAMSFLG) CALL ind_zsg (zmdl,zsg,nlvl,zt,zx,aa,bb,cc)
         ! CHG(09/11/03) with RAMS, use value assigned in advpnt, so use as input and output rather than output only
         ! so nothing to do here

         ! JCL:(8/13/01) enable vertical index to be < 1.0 for interpolation
         !        between 1st level and ground surface
         !        ZX = DMIN1(DMAX1(DBLE(1.0),ZX),DBLE(NLVL))
         ZX = DMIN1(DMAX1(DBLE(0.0),ZX),DBLE(NLVL))

         ! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
         ! JCL:(8/13/01) calculate ht of first sigma level and ht of current particle [m AGL]
         ! CHG(09/10/03) correct transformation between sigma and agl
         !        ZAGL1 = (1.0-ZSG1)*ZMDL*(ZMDL/(ZMDL-ZTER))
         ZAGL1 = (1.0-ZSG1)*(ZMDL-ZTER)
         !        ZAGL = (1.0-ZT)*ZMDL*(ZMDL/(ZMDL-ZTER))
         ZAGL = (1.0-ZT)*(ZMDL-ZTER)

         ! JCL:(8/13/01) have new interpolation subroutine (ADVINTWIND instead of ADVINT) to interpolate winds
         ! CHG(09/11/03) don't use this for RAMS

         IF (RAMSFLG) THEN

            ! CHG(09/11/03) use this for RAMS
            XXT = XT-0.5d0
            YYT = YT-0.5d0
            ! check if off grid since sometimes OFFG is not set
            IF (XXT < 1d0 .OR. YYT < 1d0 .OR. INT(XT) > NXS .OR. INT(YT) > NYS) THEN
               dead = .TRUE.
               RETURN
            END IF

            XX = DNINT(XT)                            !DNINT for scalars; for x-fluxes use X1, not XX
            YY = DNINT(YT)                            !DNINT for scalars; for y-fluxes use YT, not YY

            ! vertical grid: 1st level = 1st flux level;
            ! so assume ZX = 0.3, means need scalars at 1, and fluxes between 0 and 1
            ! so assume ZX = 0.7, means need scalars at 1, and fluxes between 0 and 1
                             !to get scalar from proper level
            ZZZ = DINT(ZX+1.0)
            ZZZ = DMIN1(DMAX1(DBLE(1.0),ZZZ),DBLE(NLVL))
            IF (XXT < 1d0  .OR.  YYT < 1d0  .OR.  XX > NXS  .OR.  YY > NYS) THEN
              dead = .TRUE.
              RETURN
            END IF

            !  interpolate to position at last time
            CALL ADVINT (U1,NXS,NYS,NZM,XXT,YY ,ZZZ,GLOBAL,NXP,NYP,UU1)
            CALL ADVINT (V1,NXS,NYS,NZM,XX ,YYT,ZZZ,GLOBAL,NXP,NYP,VV1)
            !  below 1st level using 'zero wind' boundary-condition at ground
            CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XX,YY,ZX,GLOBAL,NXP,NYP,WW1)
            CALL ADVINT (D1,NXS,NYS,NZM,XX,YY,ZZZ,GLOBAL,NXP,NYP,DD1)

            !  interpolate to position at next time
            CALL ADVINT (U2,NXS,NYS,NZM,XXT,YY ,ZZZ,GLOBAL,NXP,NYP,UU2)
            CALL ADVINT (V2,NXS,NYS,NZM,XX ,YYT,ZZZ,GLOBAL,NXP,NYP,VV2)
            CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XX,YY,ZX,GLOBAL,NXP,NYP,WW2)
            CALL ADVINT (D2,NXS,NYS,NZM,XX ,YY ,ZZZ,GLOBAL,NXP,NYP,DD2)

            ! CHG divide by density
            UU1 = UU1/DD1
            VV1 = VV1/DD1
            WW1 = WW1/DD1
            UU2 = UU2/DD2
            VV2 = VV2/DD2
            WW2 = WW2/DD2
            !  interpolate to current time
            IF (BACK) THEN
              UU = UU1
              VV = VV1
              WW = WW1
            ELSE
              UU = UU2
              VV = VV2
              WW = WW2
            END IF

            !WRITE(*,*)'adviec: TF',TF,JET
            !WRITE(45,*) ZX,WW

         ELSE IF (awrfflg) THEN

            !U(V) is staggered in x(y)-direction (C-grid) in WRF
            ! Note that x1 corresponds to position on mass grid, but staggered direction
            ! is off 0.5 (first staggered gridpoint is at position -0.5, so x1 on mass
            ! grid corresponds to x1+0.5 on staggered grid - this is different from RAMS)
            XXT = XT+0.5d0
            YYT = YT+0.5d0
            ! check if off grid since sometimes OFFG is not set
            IF (XT < 1d0  .OR.  YT < 1d0  .OR.  CEILING(XXT) > NXS  .OR.  CEILING(YYT) > NYS) THEN
               dead = .TRUE.
               RETURN
            END IF
            zxw = min(zx+dble(0.5),dble(nlvl))
            CALL ADVINTWIND (ZAGL,ZAGL1,U1,NXS,NYS,NZM,XXT,YT ,ZX ,GLOBAL,NXP,NYP,UU1)
            CALL ADVINTWIND (ZAGL,ZAGL1,V1,NXS,NYS,NZM,XT ,YYT,ZX ,GLOBAL,NXP,NYP,VV1)
            CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XT ,YT ,ZXW,GLOBAL,NXP,NYP,WW1)

            !  interpolate to position at next time
            CALL ADVINTWIND (ZAGL,ZAGL1,U2,NXS,NYS,NZM,XXT,YT ,ZX ,GLOBAL,NXP,NYP,UU2)
            CALL ADVINTWIND (ZAGL,ZAGL1,V2,NXS,NYS,NZM,XT ,YYT,ZX ,GLOBAL,NXP,NYP,VV2)
            CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XT ,YT ,ZXW,GLOBAL,NXP,NYP,WW2)

            IF (fluxflg) THEN
               !  use the time-averaged values, do not time-interpolate
               IF (BACK) THEN
                  UU = UU1
                  VV = VV1
                  WW = WW1
               ELSE
                  UU = UU2
                  VV = VV2
                  WW = WW2
               END IF
            ELSE
               !  interpolate to current time
               UU = (UU2-UU1)*TF+UU1
               VV = (VV2-VV1)*TF+VV1
               WW = (WW2-WW1)*TF+WW1
            END IF
            
         ELSE

            ! check if off grid since sometimes OFFG is not set
            IF (XT < 1d0  .OR.  YT < 1d0  .OR.  CEILING(XT) > NXS  .OR.  CEILING(YT) > NYS) THEN
               dead = .TRUE.
               RETURN
            END IF
            CALL ADVINTWIND (ZAGL,ZAGL1,U1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,UU1)
            CALL ADVINTWIND (ZAGL,ZAGL1,V1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,VV1)
            CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,WW1)

            !  interpolate to position at next time
            CALL ADVINTWIND (ZAGL,ZAGL1,U2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,UU2)
            CALL ADVINTWIND (ZAGL,ZAGL1,V2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,VV2)
            CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,WW2)

            !  interpolate to current time
            UU = (UU2-UU1)*TF+UU1
            VV = (VV2-VV1)*TF+VV1
            WW = (WW2-WW1)*TF+WW1
            
         END IF


         ! JCL:   if want to turn vertical advection (Wbar) OFF, set WW = 0.0
         !         UU = 0.0
         !         VV = 0.0
         !         WW = 0.0


         IF (IPASS == 1) THEN
            !           first pass position simple integration
            XT = X1+UU*DT
            YT = Y1+VV*DT
            ZT = Z1+WW*DT
            ! JCL:(3/1/01) store the Wbar from the 1st pass
            WWOUT = WW
         ELSE
            !           final pass average of first and second positions
            X2 = 0.5*(XT+X1+UU*DT)
            Y2 = 0.5*(YT+Y1+VV*DT)
            Z2 = 0.5*(ZT+Z1+WW*DT)
            !           final pass average of first and second ds_dz

            ! JCL:(3/1/01) average the Wbar from the two passes
            WWOUT = 0.5*(WWOUT+WW)
         END IF

      END DO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! JCL:(09/01/03) incoporate transport error by representing errors in winds as Markov process CCCCCCCCC
!      IF(WINDERRTF == 1) THEN
!
! JCL:(09/01/03) read in statistics of wind errors
!      INQUIRE(FILE='WINDERR',EXIST = FTEST)
!      IF(FTEST) THEN
!         OPEN(46,FILE='WINDERR')
!         READ(46,*) SIGUVERR !stdev of U&V errs [m/s]
!         READ(46,*) TLUVERR  !correlation Lagr timescale of U&V errs [min]
!         READ(46,*) COR.Z    !correlation lengthscale of U&V errs in vertical dir [m]
!         CLOSE(46)
!      ELSE
!         WRITE(*,*)'File "WINDERR" not found!'
!         STOP
!      END IF
!
! JCL:(09/01/03) call GASDEV for calculating Gaussian random fluctuations--represent wind errors
!      UU = GASDEV(RSEED,SIGUVERR)
!      VV = GASDEV(RSEED,SIGUVERR)
! JCL:(09/01/03) retrieve fluctuations from last timestep
!      UPRIME = SIGUVERR*UERRPREV
!      VPRIME = SIGUVERR*VERRPREV
! JCL:(09/01/03) autocorrelation function based upon time step and TL
!      RAUTO = DEXP(-1.0*DABS(DT)/TLUVERR)
!      UPRIME = RAUTO*UPRIME+DSQRT(1.0-RAUTO*RAUTO)*UU
!      VPRIME = RAUTO*VPRIME+DSQRT(1.0-RAUTO*RAUTO)*VV
! JCL:(09/01/03) conversion factor for (grid/min)=>(m/s)
!      CONVFAC = GD(INT(X2),INT(Y2))/60.0
! JCL:(09/01/03) move particle with stochastic wind error
!      X2 = X2 + (UPRIME/CONVFAC)*DT
!      Y2 = Y2 + (VPRIME/CONVFAC)*DT
! JCL:(09/01/03) store displacement resulting from stochastic wind error; will be applied in HYMODELC
!      DXERR = (UPRIME/CONVFAC)*DT
!      DYERR = (VPRIME/CONVFAC)*DT
! JCL:(09/01/03) update the stored wind fluctuations
!      UERRPREV = UPRIME/SIGUVERR
!      VERRPREV = VPRIME/SIGUVERR
!
! JCL:(09/01/03) IF(WINDERRTF == 1) THEN
!      END IF
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! JCL:(3/1/01) turn off calculations of DUDZ & DVDZ b/c not outputing these variables
!CCCCCCCCC
! JCL:get vertical position at the new time
!     height agl when zsfc = 0
!      ZZ = ZMDL*(1.0-DMIN1(DBLE(1.0),Z2))
!     relation between index and sigma when zsfc = 0
!      ZX = (-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
! jcl comment: makes sure that 1.0<=ZX<=NLVL
!      ZX = DMIN1(DMAX1(DBLE(1.0),ZX),DBLE(NLVL))

! JCL:calculate the altitudes [m] at the two height indices bracketing particle location
!      ZLO = (1/(4.0*AA))*(((AINT(ZX)*2.0*AA+BB)**2)-BB**2)+CC
!      ZHI = (1/(4.0*AA))*((((AINT(ZX)+1)*2.0*AA+BB)**2)-BB**2)+CC

! JCL:get horizontal velocities @ hts just below&above new location
!      CALL ADVINT (U2,NXS,NYS,NZM,X2,Y2,DINT(ZX),ULO)
!      CALL ADVINT (V2,NXS,NYS,NZM,X2,Y2,DINT(ZX),VLO)
!      CALL ADVINT (U2,NXS,NYS,NZM,X2,Y2,DINT(ZX)+1.0,UHI)
!      CALL ADVINT (V2,NXS,NYS,NZM,X2,Y2,DINT(ZX)+1.0,VHI)

! JCL:convert velocities from (grid/min)=>(m/s)
! JCL:conversion factor for (grid/min)=>(m/s)
!      CONVFAC = GD(INT(X2),INT(Y2))/60.0
!      ULO = ULO*CONVFAC
!      UHI = UHI*CONVFAC
!      VLO = VLO*CONVFAC
!      VHI = VHI*CONVFAC
! JCL:calculate the wind shears
!      DUDZ = (UHI-ULO)/(ZHI-ZLO)
!      DVDZ = (VHI-VLO)/(ZHI-ZLO)
!CCCCCCCCC

      RETURN
      END
