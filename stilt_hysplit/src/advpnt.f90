!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVPNT           ADVection of one PoiNT in space
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION OF ONE POINT IN SPACE IS THE PRIMARY ROUTINE THAT IS U
!   TRAJECTORY AND DISPERSION SIMULATIONS.  IT IS CALLED EACH TIME A
!   IS REQUIRED TO BE ADVECTED OVER A DURATION OF ONE TIME STEP.
!   ROUTINE CHECKS METEO DATA IN ARRAY, IF POINT FITS WITHIN THE TIM
!   SPACE LIMITS CALCULATION CONTINUES, OTHERWISE NEW DATA ARE INPUT
!   IN ADDITION METEO VARIABLES AT THE END OF THE STEP ARE PLACED IN
!   METO%STRUCTURE.  THOSE VALUES ARE USED IN SUBSEQUENT DISPERSION
!   TRAJECTORY OUTPUT ROUTINES.  NOTE THAT ARRAY DIMENSIONS ARE AT T
!   COMPILED MAXIMUM IN THIS UPPER LEVEL ROUTINE.  ALL ROUTINES CALL
!   HERE REQUIRE SUB-GRID DIMENSIONS, EXCEPT IN THE VERTICAL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 10 Apr 1998 (RRD)
!                 18 Aug 1998 (RRD) - added isotropic turbulence option
!                 27 Oct 1998 (RRD) - Clean return (KG=0) on no met data
!                 21 Dec 1998 (RRD) - added endrec to metinp call
!                                     new subroutine to compute horiz mixing
!                                     added rdep to argument list
!                 25 Feb 1999 (RRD) - corrected ftime initialization
!                 04 Mar 1999 (RRD) - changed argument list to metpos
!                 19 Apr 1999 (RRD) - added terrain height array
!
! USAGE:  CALL ADVPNT(BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,ISOT,
!              XP,YP,ZP,JET,DT,KG,KGRID,NGRD,ZSG,NLVL,ZMDL,KVEL,
!              UBAR,IFHR,*)
!   INPUT ARGUMENT LIST:
!     BACK      - log   backward direction integration flag
!     VMIX      - log   flag to return mixing profile
!     CDEP      - log   flag to return deposition variables
!     RDEP      - log   flag to indicate that resistance deposition is on
!     TRAJ      - log   flag to return trajectory variables
!     KSFC      - int   index for the height of the sfc layer top
!     ISOT      - int   flag to set isotropic turbulence (0-no; 1-yes)
!     XP,YP,ZP  - real  position before advection
!     JET       - int   elapsed time (min)
!     DT        - real  advection time step (min)
!     KG        - int   grid selection indicator for calculation
!     KGRID     - int   marker index to indicate current open file
!     NGRD      - int   number of meteo grids
!     ZSG       - real  vertical sigma levels
!     NLVL      - int   number of vertical levels
!     ZMDL      - real  model top m agl
!     KVEL      - int   vertical velocity remapping option
!   OUTPUT ARGUMENT LIST:
!     XP,YP,ZP  - real  updated position after advection
!     UBAR      - real  advection velocity for component
!     IFHR      - int   current forecast hour
!     *         - alternate return for point off grid
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT 30 - diagnostic messages
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: advpnt.f90,v 1.22 2009-03-05 08:26:22 gerbig Exp $
!
!$$$

! JCL:(3/1/01)add WWOUT as output argument
! JCL:add DUDZ & DVDZ as output arguments--CHG needs them for box
! JCL:(9/16/02) NHRSZI and ZIPRESC are used to scale mixed-layer heights
! CHG:(12/05/01) add 'CONVDUR' (reset to 0 in PARDSP after conv. redistribution was done)
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(09/01/03) add wind error flag, UPRIME/UERR and VPRIME/VERR from previous timestep
! JCL:(09/01/03) add random seed 'RSEED' for modeling transport error as stochastic process
! JCL:(09/01/03) DXERR and DYERR are the horizontal displacements resulting from transport error
! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! CHG(09/18/03) pass on RAMSFLG
      SUBROUTINE ADVPNT(BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,ISOT,             &
     &   XP,YP,ZP,JET,DT,KG,KGRID,NGRD,ZSG,NLVL,ZMDL,KVEL,              &
     &   UBAR,IFHR,DUDZ,DVDZ,WWOUT,ZICONTROLTF,NHRSZI,ZIPRESC,          &
     &   CONVDUR,ICONVECT,RAMSFLG,ECMFLG,is_off_grid)

      USE module_defgrid
      USE module_defmeto
      USE module_defconc

      IMPLICIT REAL*8 (A-H,O-Z)


      logical is_off_grid !replaces alternate return mechanism

!     meteo variables defined by (x,y,z,t), where t varies from
!     1 (last observation time) to 2 (next observation time)

!     forecast hour associated with data
      INTEGER FHOUR(2)
!     horizontal wind components (input m/s; output grid/min)
      REAL*8 U(NXM,NYM,NZM,2), V(NXM,NYM,NZM,2)
!     vertical velocity (input dp/dt; output sigma/min)
      REAL*8 W(NXM,NYM,NZM,2)
!     ambient temperature (input deg-K; output pot-K)
      REAL*8 T(NXM,NYM,NZM,2)
!     moisture (input rh or kg/kg; output rh/100)
      REAL*8 Q(NXM,NYM,NZM,2)
!     local air density (kg/m3),local air pressure (mb)
      REAL*8 D(NXM,NYM,NZM,2), P(NXM,NYM,NZM,2)
!     vertical and horizontal mixing coefficients (m2/s)
      REAL*8 X(NXM,NYM,NZM,2), H(NXM,NYM,NZM,2)

! JCL(03/27/03): add arrays to store grids of TL & SIGW from RAMS
! note: the TLRAMS array is also used for awrfflg=T and fluxflg=T
      REAL*8 TLRAMS(NXM,NYM,NZM,2), SIGWRAMS(NXM,NYM,NZM,2)

      !  arrays to store convective fluxes, zero initialization for fluxes not in ECMWF input
      !  Radius up/downdraft computed later
      REAL*8, SAVE :: CFXUP1(NXM,NYM,NZM,2), CFXUP2(NXM,NYM,NZM,2), &
                      CFXDN1(NXM,NYM,NZM,2),                            &
                      DFXUP1(NXM,NYM,NZM,2), DFXUP2(NXM,NYM,NZM,2), &
                      DFXDN1(NXM,NYM,NZM,2),                            &
                      EFXUP1(NXM,NYM,NZM,2), EFXUP2(NXM,NYM,NZM,2), &
                      EFXDN1(NXM,NYM,NZM,2)

! CHG(09/25/03) add RAMS turb. kin. energy TKEN
      REAL*8 TKEN(NXM,NYM,NZM,2)

! JCL:Lagrangian timescale (s)
      REAL*8 TL(NXM,NYM,NZM,2)
! JCL:stddev of vertical velocity (m/s)
      REAL*8 SIGW(NXM,NYM,NZM,2)

!     model sfc press (mb), rainfall total (m)
! CHG:(11/20/01) add rainfall convective
      REAL*8 P0(NXM,NYM,2), RT(NXM,NYM,2), RC(NXM,NYM,2)
! CHG:(11/20/01) add cloud cover and swhortw. flux
      REAL*8 TC(NXM,NYM,2), SW(NXM,NYM,2)
! CHG:(12/04/01) add low cover
      REAL*8 LC(NXM,NYM,2)
! CHG:(22/01/03) add soil moisture
      REAL*8 SM(NXM,NYM,2)
!     low-level horizontal wind (m/s), ambient temp (deg-K)
      REAL*8 U0(NXM,NYM,2), V0(NXM,NYM,2), T0(NXM,NYM,2), W0(NXM,NYM,2), &
           & HPBL(NXM,NYM,2)
! additional variables for un-coupling WRF momentum fluxes
! (NOTE: we will reuse the TLRAMS array for the dry inverse density ALT)
      real*8 muu(NXM,NYM,2),muv(NXM,NYM,2),mu(NXM,NYM,2)
      real*8 msfu(NXM,NYM,2),msfv(NXM,NYM,2),msft(NXM,NYM,2)
!     sfc fluxes: u,v momentum (n/m2) or (kg/m2-s), heat (w/m2)
!     u momentum flux gets replaced by U* in prfcom
      REAL*8 UF(NXM,NYM,2), VF(NXM,NYM,2), SF(NXM,NYM,2)
! CHG:(11/20/01) add sfc flux for water (w/m2)
      REAL*8 LF(NXM,NYM,2)

!     grid dist (m), roughness length (m), terrain elevation (m)
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!     REAL*8 GD(NXM,NYM), Z0(NXM,NYM), ZT(NXM,NYM)
      REAL*8 GX(NXM,NYM), GY(NXM,NYM), Z0(NXM,NYM), ZT(NXM,NYM)

! JCL:(4/28/00)add mixed-layer height (m)
      REAL*8 ZML(NXM,NYM,2)

! CHG:(12/04/01)add limit of convection height (m)
      REAL*8 ZLOC(NXM,NYM,2)

! JCL:(4/3/02)mass violation grid (fraction/min)
      REAL*8 DMASS(NXM,NYM,NZM,2)
! JCL:(4/3/02)temp variables to hold density
      REAL*8 DENS1,DENS2

!     land use category
      INTEGER LU(NXM,NYM)

!     internal model sigma levels
      REAL*8 ZSG(NLVL)

! JCL:(10/31/02)variables to store profiles of turbulence variables for the two analysis times
!     REAL*8 TLK1(NLVL),TLK2(NLVL),SIGWK1(NLVL),SIGWK2(NLVL)
      REAL*8 TLK1(NZM),TLK2(NZM),SIGWK1(NZM),SIGWK2(NZM)

! JCL:(10/31/02)heights at different sigma levels [m AGL]
      REAL*8 ZLVLS(NZM)

! JCL: add DUDZ, DVDZ--the wind shear--CHG needs it to parameterize box
      REAL*8 DUDZ, DVDZ

! CHG:(12/04/01) add duration for conv. redistribution
      INTEGER CONVDUR

! JCL:(9/16/02) variables needed to prescribe ZI
      REAL*8  ZIPRESC(150)
      INTEGER ZICONTROLTF,GRIDHR

! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! JCL:(09/01/03) WINDERRTF specifies whether to include wind errs as Markov process (=1) or not (=0)
!      INTEGER WINDERRTF
! JCL:(09/01/03) UPRIME/UERR and VPRIME/VERR from previous timestep
!      REAL*8 UERRPREV,VERRPREV
! JCL:(09/01/03) add random seed 'RSEED' for modeling transport error as stochastic process
!      INTEGER RSEED

!     input data times: requested and read from file
      INTEGER MTIME(2), FTIME(2)

!     integration direction, off-grid, sub-grid, mixing, deposit, traj, particle dead (special case of off-grid)
      LOGICAL BACK, OFFG, NEWX, VMIX, CDEP, RDEP, TRAJ, DEAD

! CHG(09/11/03):add flag specifying whether data from RAMS or not
      LOGICAL, INTENT(IN) :: RAMSFLG, ECMFLG
      LOGICAL awrfFLG, fluxflg, deepflg, shallflg

! JCL:add the two new variables TL & SIGW to be saved
! CHG:(11/20/01) add rainfall convective
! CHG:(11/20/01) add cloud cover and shortw. flux
! CHG:(11/20/01) add sfc flux for water (w/m2)
! CHG:(12/04/01) add low cloud cover and ZLOC
! JCL:(4/7/02) add array of mass violation to be saved
! CHG:(22/01/03) add soil moisture SM
! JCL(03/27/03): add RAMS turbulence variables
! CHG(09/23/03) add convective fluxes from RAMS
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
!     save definitions and meteo data between calls
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
      integer :: icall
      SAVE NXS,NYS,MTIME,FTIME,FHOUR,JETP,                              &
     &     LU,Z0,GX,GY,UF,VF,SF,P0,RT,U0,V0,W0,T0,hpbl,U,V,W,T,Q,D,P,X,H,  &
     &     TL,SIGW,ZML,RC,LF,TC,LC,SW,ZLOC,DMASS,SM,TLRAMS,SIGWRAMS,    &
           TKEN,icall,ZT
! needed for WRF model input
      save muu,muv,mu,msfu,msfv,msft,fluxflg, deepflg, shallflg

!     initialize special diagnostic  counters
      DATA JETP, NOFFH, NOFFV, NOFFT, icall / 0, 0, 0, 0, 0 /

!---------------------------------------------------------------------------------------------------
!      WRITE(*,*)'Entered ADVPNT'

      if (icall .eq. 0) then
         CFXUP2(:,:,:,:)=0d0
         DFXUP2(:,:,:,:)=0d0
         EFXUP2(:,:,:,:)=0d0
      endif

      icall=icall+1

      is_off_grid = .FALSE.

!     initialize array data time from file read to no data
      IF(KGRID.EQ.0)THEN
         FTIME(1)=0
         FTIME(2)=0
      END IF

!     diagnostic output when time changes
      IF(JETP.NE.JET)THEN
!        IF(NOFFH.GT.0)WRITE(30,*)'NOTICE advpnt - off grid:',NOFFH
!        IF(NOFFV.GT.0)WRITE(30,*)'NOTICE advpnt - out top :',NOFFV
!        IF(NOFFT.GT.0)WRITE(30,*)'NOTICE advpnt - out time:',NOFFT
         NOFFH=0
         NOFFV=0
         NOFFT=0
         JETP=JET
      END IF

!     determine meteo file record numbers and sub-grid location
! JEL(03/06/2007) Save values in case we end up ignoring grid expansion
      NXSB=NXS  
      NYSB=NYS  
      LX1B=METO%LX1
      LY1B=METO%LY1
      KGRIDB=KGRID
      CALL METPOS(BACK,XP,YP,JET,KG,NGRD,NXS,NYS,                       &
     &   METO%LX1,METO%LY1,METO%LXC,METO%LYC,METO%LXR,METO%LYR,         &
     &   KREC1,KREC2,OFFG,NEWX,MTIME,FTIME,KGRID)

! JCL:
!      WRITE(45,*)'####',JET+INT(DT)-(CONC(1)%START%MACC)
!     &           ,JET+INT(DT)
!      WRITE(45,*)MTIME(1),MTIME(2)
!      WRITE(45,*)FTIME(1),FTIME(2)


!     location not on current meteo grid
!     set default for KGOLD
      KGOLD = 0
      DO WHILE (OFFG)
!        when more grids available determine why calculation offgrid
         IF(KG.LT.NGRD)THEN
            IF(MTIME(2).LT.FILE(KG+1,1)%FIRST%MACC.OR.                  &
     &         MTIME(2).GT.FILE(KG+1,1)%LAST%MACC)THEN
!              terminate calculation for point, but stay on same grid
!              because there is no new data at this time
! CHG(09/18/03) keep info about old KG
               KGOLD=KG
               KG=0
               NOFFT=NOFFT+1
               is_off_grid = .TRUE.
! JEL(03/06/2007) Ignore grid expansion
               NXS=NXSB
               NYS=NYSB
               METO%LX1=LX1B
               METO%LY1=LY1B
               KGRID=KGRIDB
               RETURN
            END IF

! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
!           convert position from old grid to true
            IF(GRID(KG)%LATLON)THEN
              CALL GBL2LL(KG,XP,YP,TLAT,TLON)
            ELSE
              CALL CXY2LL(GRID(KG)%GBASE,XP,YP,TLAT,TLON, GRID(KG)%proj)
            END IF

!           grid number increases - no going back (although possible)
            KG=KG+1

! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
!           convert true back to new grid coordinates
            IF(GRID(KG)%LATLON)THEN
              CALL GBL2XY(KG,TLAT,TLON,XP,YP)
            ELSE
              CALL CLL2XY(GRID(KG)%GBASE,TLAT,TLON,XP,YP,GRID(KG)%proj)
            END IF

! JCL:(02/18/2005) added Draxler's fix in HYSPLIT Vers. 4.5
!           when switching grids redefine subgrid center (RRD: 1/3/00)
            METO%LXC=XP
            METO%LYC=YP

!           keep checking grids until no new grids left
            CALL METPOS(BACK,XP,YP,JET,KG,NGRD,NXS,NYS,                 &
     &         METO%LX1,METO%LY1,METO%LXC,METO%LYC,METO%LXR,METO%LYR,   &
     &         KREC1,KREC2,OFFG,NEWX,MTIME,FTIME,KGRID)

!           update current active grid index
            KGRID=KG

         ELSE
!           when on last grid terminate calculation at this point
            NOFFH=NOFFH+1
! CHG(09/18/03) keep info about old KG
            KGOLD=KG
            KG=0
            IF(KREC1.EQ.0.AND.KREC2.EQ.0)THEN
               is_off_grid = .TRUE.
               RETURN
            END IF
!           if more data are required continue loading data
            OFFG=.FALSE.
         END IF
      END DO


!     test for postion off final meteo grid (terminate)
      IF(ZP.LT.ZSG(NLVL))THEN
         NOFFV=NOFFV+1
! CHG(09/18/03) keep info about old KG
         KGOLD=KG
         KG=0
         IF(KREC1.EQ.0.AND.KREC2.EQ.0)THEN
            is_off_grid = .TRUE.
            RETURN
         END IF
      END IF

!     test for end of meteorology file (no more data)
      IF((FILE(KGRID,2)%PERIOD).EQ.0)THEN
!        correction 10/27/98
! CHG(09/18/03) keep info about old KG
         KGOLD=KG
         KG=0
         WRITE(30,*)'WARNING advpnt - no more meteo data'
         is_off_grid = .TRUE.
         RETURN
      END IF

!     confirm # data levels doesn't exceed compiled array limit
!     NZS defines the vertical input data subgrid
      NZS=MIN0(MLVL,NZM,GRID(KGRID)%NZ-1)

      !  use different NXY value for RAMS/ECMWF(X) to allow 16 bit instead of 8
      NXY = GRID(KGRID)%NX*GRID(KGRID)%NY
      IF (RAMSFLG .OR. GRID(KGRID)%model_id(1:3) == 'ECX' &
           & .OR. GRID(KGRID)%model_id == 'DWRF' ) NXY = NXY*2
      awrfflg = GRID(KGRID)%model_id(2:4) .eq. 'WRF'

!     set meteo indices
      K1=FILE(KGRID,1)%PERIOD
      K2=FILE(KGRID,2)%PERIOD
!     unit numbers for data
      KUNIT1=FILE(KGRID,1)%KUNIT
      KUNIT2=FILE(KGRID,2)%KUNIT
!     last valid record in each file
      KEND1=FILE(KGRID,1)%ENDREC
      KEND2=FILE(KGRID,2)%ENDREC


!     shift in sub-grid position requires variable reloading
      IF (NEWX) CALL METGRD(KGRID,METO%LX1,METO%LY1,NXS,NYS,GX,GY,Z0,LU)


!     test if new data required at the last time (k1)
      IF(KREC1.GT.0)THEN

!      WRITE(45,*)'KREC1:',(JET+INT(DT)-(CONC(1)%START%MACC))/60
!     &           ,KREC1,K1
!        load meteorological data according to positioning
!        reads data for NZS levels
! CHG:(11/20/01)added conv. precip. rates (RC) to arguments
! CHG:(11/20/01)add tot. cloud cover and shortw. flux
! CHG:(11/20/01) add sfc flux for water (w/m2)
! CHG:(12/04/01) add low cloud cover
! CHG:(22/01/03) add soil moisture
! JCL(03/27/03): add arrays to store grids of TL & SIGW from RAMS
! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
         CALL METINP(BACK,KGRID,KUNIT1,KREC1,NXY,METO%LX1,METO%LY1,     &
     &      NXS,NYS,NZS,FTIME(1),KEND1,FHOUR(K1),ZT,                    &
     &      P0(:,:,K1),T0(:,:,K1),U0(:,:,K1),V0(:,:,K1),W0(:,:,K1),hpbl(:,:,K1),     &
     &      muu(:,:,K1),muv(:,:,K1),mu(:,:,K1),                         &
     &      msfu(:,:,K1),msfv(:,:,K1),msft(:,:,K1),fluxflg, deepflg, shallflg,             &
     &      UF(:,:,K1),VF(:,:,K1),SF(:,:,K1),RT(:,:,K1),                &
     &      U(:,:,:,K1),V(:,:,:,K1),W(:,:,:,K1),T(:,:,:,K1),            &
     &      Q(:,:,:,K1),P(:,:,:,K1),RC(:,:,K1),LF(:,:,K1),              &
     &      TC(:,:,K1),LC(:,:,K1),SW(:,:,K1),SM(:,:,K1),                &
     &      TLRAMS(:,:,:,K1),SIGWRAMS(:,:,:,K1),CFXUP1(:,:,:,K1),       &
     &      CFXUP2(:,:,:,K1),CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),         &
     &      DFXUP2(:,:,:,K1),EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),         &
     &      DFXDN1(:,:,:,K1),EFXDN1(:,:,:,K1),TKEN(:,:,:,K1))
!        WRITE(45,*)'advpnt, K1,K2:',K1,K2
         ZISCALE=1.0
! JCL:(9/16/02) if PBL height is prescribed
         IF(ZICONTROLTF.EQ.1)THEN
            if (nhrszi .eq. 1 .and. zipresc(nhrszi) .lt. 0.0) then
! Special case to denote using hpbl from met file
               ziscale = -1.
            else
! JCL:      time of loaded met grid (# of hrs elapsed since model starting time)
               GRIDHR=ABS((MTIME(1)-(CONC(1)%START%MACC))/60)
               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(GRIDHR+1)
            endif
         END IF

!        vertical interpolation to terrain following coordinates
!        data from NZS levels processed to NLVL of internal grid
! JCL:   add 'TL' & 'SIGW' as arguments to PRFCOM
! JCL:(4/28/00)add ZML (mixed-layer ht) as output from PRFCOM
! CHG:(12/04/04)add ZLOC (lim. of convection ht) as output from PRFCOM
! JCL:(4/3/02)add mass violation grid as output from PRFCOM
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer height
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        CALL PRFCOM(ISOT,KGRID,KSFC,GD,Z0,ZT,
         CALL PRFCOM(ISOT,KGRID,KSFC,GX,GY,Z0,ZT,                       &
     &      NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,ICONVECT,                    &
     &      P0(:,:,K1),T0(:,:,K1),U0(:,:,K1),V0(:,:,K1),W0(:,:,K1),hpbl(:,:,K1),UF(:,:,K1),     &
     &      VF(:,:,K1),SF(:,:,K1),U(:,:,:,K1),V(:,:,:,K1),              &
     &      W(:,:,:,K1),T(:,:,:,K1),Q(:,:,:,K1),                        &
     &      P(:,:,:,K1),D(:,:,:,K1),X(:,:,:,K1),TL(:,:,:,K1),           &
     &      SIGW(:,:,:,K1),ZML(:,:,K1),ZLOC(:,:,K1),DMASS(:,:,:,K1),    &
     &      ZISCALE,TLRAMS(:,:,:,K1),mu(:,:,K1),                        &
     &      muu(:,:,K1),muv(:,:,K1),msfu(:,:,K1),msfv(:,:,K1),msft(:,:,K1), &
     &      fluxflg, deepflg, shallflg, &
     &      CFXUP1(:,:,:,K1),       &
     &      CFXUP2(:,:,:,K1),CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),         &
     &      DFXUP2(:,:,:,K1),EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),         &
     &      DFXDN1(:,:,:,K1),EFXDN1(:,:,:,K1),TKEN(:,:,:,K1))
! JCL:(03/27/03) assign SIGW & TL outputted directly from RAMS
!         IF(GRID(KG)%MODEL_ID.EQ.'RAMS')THEN
!            SIGW(1:NXM,1:NYM,1:NZM,K1)=SIGWRAMS(1:NXM,1:NYM,1:NZM,K1)
!            TL(1:NXM,1:NYM,1:NZM,K1)=TLRAMS(1:NXM,1:NYM,1:NZM,K1)
!         END IF
! JCL:
!     loop over each grid position
!      DO K=1,NLVL
!      DO J=1,NYS-1
!      DO I=1,NXS-1
!         WRITE(45,*)I,J,K,SIGWRAMS(I,J,K),TLRAMS(I,J,K)
!      END DO
!      END DO
!      END DO

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        compute horizontal mixing coefficient
         IF(VMIX) CALL STBHOR(NXS,NYS,NLVL,GX,GY,                       &
     &            U(:,:,:,K1),V(:,:,:,K1),H(:,:,:,K1))

!        reset meteo time position if data are missing
         IF(FTIME(1).NE.MTIME(1))THEN
            WRITE(30,*)'WARNING advpnt: meteo time-1 off'
            WRITE(30,*)'Internal time:  ',JET
            WRITE(30,*)'Input data time:',FTIME(1)
            WRITE(30,*)'Expected time:  ',MTIME(1)
            MTIME(1)=FTIME(1)
         END IF
      END IF


!     test if new data required at the next time (k2)
      IF(KREC2.GT.0)THEN

!      WRITE(*,*)'KREC2:',(JET+INT(DT)-(CONC(1)%START%MACC))/60
!     &           ,KREC2,K2

!        load data for subgrid
! CHG:(11/20/01)added conv. precip. rates (RC) to arguments
! CHG:(11/20/01)add tot. cloud cover and shortw. flux
! CHG:(11/20/01) add sfc flux for water (w/m2)
! CHG:(12/04/01) add low cloud cover
! CHG:(22/01/03) add soil moisture
! JCL(03/27/03): add arrays to store grids of TL & SIGW from RAMS
! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
         CALL METINP(BACK,KGRID,KUNIT2,KREC2,NXY,METO%LX1,METO%LY1,     &
     &      NXS,NYS,NZS,FTIME(2),KEND2,FHOUR(K2),ZT,                    &
     &      P0(:,:,K2),T0(:,:,K2),U0(:,:,K2),V0(:,:,K2),W0(:,:,K2),hpbl(:,:,K2),                &
     &      muu(:,:,K2),muv(:,:,K2),mu(:,:,K2),                         &
     &      msfu(:,:,K2),msfv(:,:,K2),msft(:,:,K2),fluxflg, deepflg, shallflg,             &
     &      UF(:,:,K2),VF(:,:,K2),SF(:,:,K2),RT(:,:,K2),                &
     &      U(:,:,:,K2),V(:,:,:,K2),W(:,:,:,K2),T(:,:,:,K2),            &
     &      Q(:,:,:,K2),P(:,:,:,K2),RC(:,:,K2),LF(:,:,K2),              &
     &      TC(:,:,K2),LC(:,:,K2),SW(:,:,K2),SM(:,:,K2),                &
     &      TLRAMS(:,:,:,K2),SIGWRAMS(:,:,:,K2),CFXUP1(:,:,:,K2),       &
     &      CFXUP2(:,:,:,K2),CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),         &
     &      DFXUP2(:,:,:,K2),EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),         &
     &      DFXDN1(:,:,:,K2),EFXDN1(:,:,:,K2),TKEN(:,:,:,K2))

         ZISCALE=1.0
! JCL:(9/16/02) if PBL height is prescribed
         IF(ZICONTROLTF.EQ.1)THEN
            if (nhrszi .eq. 1 .and. zipresc(nhrszi) .lt. 0.0) then
! Special case to denote using hpbl from met file
               ziscale = -1.
            else
! JCL:      time of loaded met grid (# of hrs elapsed since model starting time)
               GRIDHR=ABS((MTIME(2)-(CONC(1)%START%MACC))/60)
               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(GRIDHR+1)
            endif
         END IF

! JCL:   add 'TL' & 'SIGW' as arguments to PRFCOM
! JCL:(4/28/00)add ZML (mixed-layer ht) as output from PRFCOM
!        vertical interpolation of profile
! CHG:(12/04/02)add ZLOC (lim. of convection ht) as output from PRFCOM
! JCL:(4/3/02)add mass violation grid as output from PRFCOM
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer height
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        CALL PRFCOM(ISOT,KGRID,KSFC,GD,Z0,ZT,
         CALL PRFCOM(ISOT,KGRID,KSFC,GX,GY,Z0,ZT,                       &
     &      NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,ICONVECT,                    &
     &      P0(:,:,K2),T0(:,:,K2),U0(:,:,K2),V0(:,:,K2),W0(:,:,K2),hpbl(:,:,K2),UF(:,:,K2),     &
     &      VF(:,:,K2),SF(:,:,K2),U(:,:,:,K2),                          &
     &      V(:,:,:,K2),                                                &
     &      W(:,:,:,K2),T(:,:,:,K2),Q(:,:,:,K2),                        &
     &      P(:,:,:,K2),D(:,:,:,K2),X(:,:,:,K2),TL(:,:,:,K2),           &
     &      SIGW(:,:,:,K2),ZML(:,:,K2),ZLOC(:,:,K2),DMASS(:,:,:,K2),    &
     &      ZISCALE,TLRAMS(:,:,:,K2),mu(:,:,K2),                        &
     &      muu(:,:,K2),muv(:,:,K2),msfu(:,:,K2),msfv(:,:,K2),msft(:,:,K1), &
     &      fluxflg, deepflg, shallflg, &
     &      CFXUP1(:,:,:,K2),       &
     &      CFXUP2(:,:,:,K2),CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),         &
     &      DFXUP2(:,:,:,K2),EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),         &
     &      DFXDN1(:,:,:,K2),EFXDN1(:,:,:,K2),TKEN(:,:,:,K2))

! JCL:(03/27/03) assign SIGW & TL outputted directly from RAMS
!         IF(GRID(KG)%MODEL_ID.EQ.'RAMS')THEN
!            SIGW(1:NXM,1:NYM,1:NZM,K2)=SIGWRAMS(1:NXM,1:NYM,1:NZM,K2)
!            TL(1:NXM,1:NYM,1:NZM,K2)=TLRAMS(1:NXM,1:NYM,1:NZM,K2)
!         END IF
!      DO KKKK=1,NZM
!      DO JJJJ=1,NYS
!      DO IIII=1,NXS
!       IF(iiii.eq.25.and.jjjj.eq.25)then
!         WRITE(*,*)'grids DN TL TLR:',IIII,JJJJ,KKKK,
!     &   D(IIII,JJJJ,KKKK,K1),TL(IIII,JJJJ,KKKK,K1),
!     &   TLRAMS(IIII,JJJJ,KKKK,K1)
!       end if
!      END DO
!      END DO
!      END DO
!      WRITE(*,*)'K1:',K1

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        compute horizontal mixing coefficient
         IF(VMIX) CALL STBHOR(NXS,NYS,NLVL,GX,GY,                       &
     &            U(:,:,:,K2),V(:,:,:,K2),H(:,:,:,K2))


!        missing data metpos will wait longer before reading
         IF(FTIME(2).NE.MTIME(2))THEN
            WRITE(30,*)'WARNING advpnt: meteo time-2 off'
            WRITE(30,*)'Internal time:  ',JET
            WRITE(30,*)'Input data time:',FTIME(2)
            WRITE(30,*)'Expected time:  ',MTIME(2)
            MTIME(2)=FTIME(2)
         END IF
      END IF

!     determine minutes between data time periods for interpolation
      DM=IABS(MTIME(2)-MTIME(1))
      IF(DM.GT.720)THEN
         WRITE(30,*)'ERROR advpnt: excessive interpolation'
         WRITE(30,*)'Greater than 720 min: ',DM
         STOP
      END IF

!     optional vertical velocity remapping
      IF(KREC1.GT.0.AND.KVEL.GT.0)                                      &
     &   CALL METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,                      &
     &   U(:,:,:,K1),V(:,:,:,K1),W(:,:,:,K1),P(:,:,:,K1),P(:,:,:,K2),   &
     &   T(:,:,:,K1),T(:,:,:,K2),D(:,:,:,K1),D(:,:,:,K2))

      IF(KREC2.GT.0.AND.KVEL.GT.0)                                      &
     &   CALL METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,                      &
     &   U(:,:,:,K2),V(:,:,:,K2),W(:,:,:,K2),P(:,:,:,K1),P(:,:,:,K2),   &
     &   T(:,:,:,K1),T(:,:,:,K2),D(:,:,:,K1),D(:,:,:,K2))

!     map position from meteo grid to sub-grid
      X1=XP-METO%LX1+1
      Y1=YP-METO%LY1+1
      Z1=ZP

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! CHG(09/11/03) flag specifying whether data from RAMS or not
!      RAMSFLG=GRID(KG)%MODEL_ID.EQ.'RAMS'


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! JCL4/13/00: following lines calculate the vertical profiles of TL & SIGMAW
!     at location BEFORE mean advection--will be used in PARDSP, along with
!     vertical profiles AFTER mean advection, for vertical turbulent transport
! JCL: a lot of the code below are similar to lines in ADVMET

! JCL:(4/25/00)have lines diff. from those in ADVMET b/c interpolation
!        factor for before advection should have value of 0 when subgrid is reloaded,
!        instead of value of 1 for the case after advection
!     interpolation factor for current time
      IF (ABS(DM) <= TINY(DM)) STOP 'ADVPNT: DM == 0. STOP.'
      TF = MOD(DBLE(JET),DM)/DM
      IF (BACK .AND. (ABS(TF) >= TINY(TF))) TF=1d0-TF

! CHG:(12/04/01) set duration for conv. redistribution,
!     to be done once for every analysis time
         IF (ABS(TF) <= TINY(TF)) CONVDUR=DM
         IF(JET.EQ.CONC(1)%START%MACC)CONVDUR=IDNINT(DM*(1d0-TF))

! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
      DO KL=1,NLVL
!        set vertical interpolation point to index position
         ZK=DBLE(KL)

! JCL:   interpolate Lagrangian timescale to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
         XX=DNINT(X1)
         YY=DNINT(Y1)
         CALL ADVINT(TL(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADVINT(TL(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
! JCL:(10/31/02) more sophisticated temporal interpolation that involves scaling by local zi
!         METO%TLPREV(KL)=(VAR2-VAR1)*TF+VAR1
         TLK1(KL)=VAR1
         TLK2(KL)=VAR2

! JCL:   interpolate std dev of vertical velocity to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
         XX=DNINT(X1)
         YY=DNINT(Y1)
         CALL ADVINT(SIGW(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,               &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADVINT(SIGW(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,               &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
! JCL:(10/31/02) more sophisticated temporal interpolation that involves scaling by local zi
!         METO%SIGWPREV(KL)=(VAR2-VAR1)*TF+VAR1
         SIGWK1(KL)=VAR1
         SIGWK2(KL)=VAR2
!         WRITE(45,*)"TLSIGW",KL,TLK1(KL),TLK2(KL),SIGWK1(KL),SIGWK2(KL)

! JCL:(5/9/01)interpolate horizontal density profiles to position
         if (awrfflg) then
                   !U(V) is staggered in x(y)-direction (C-grid) in WRF
! Note that x1 corresponds to position on mass grid, but staggered direction
! is off 0.5 (first staggered gridpoint is at position -0.5, so x1 on mass
! grid corresponds to x1+0.5 on staggered grid - this is different from RAMS)
            XX1=X1+0.5
            YY1=Y1+0.5
            if (.not. fluxflg) then
! instantantenous velocities
               CALL ADVINT(U(:,:,:,K1),NXS,NYS,NZM,XX1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
               CALL ADVINT(U(:,:,:,K2),NXS,NYS,NZM,XX1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               METO%UUPREV(KL)=(VAR2-VAR1)*TF+VAR1
               CALL ADVINT(V(:,:,:,K1),NXS,NYS,NZM,X1,YY1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
               CALL ADVINT(V(:,:,:,K2),NXS,NYS,NZM,X1,YY1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               METO%VVPREV(KL)=(VAR2-VAR1)*TF+VAR1
            else
! time-averaged coupled u,v (decoupled in prfcom):
               IF(BACK)KLATE=K1
               IF(.NOT.BACK)KLATE=K2
               CALL ADVINT(U(:,:,:,KLATE),NXS,NYS,NZM,XX1,Y1,ZK,             &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               METO%UUPREV(KL)=VAR2
               CALL ADVINT(V(:,:,:,KLATE),NXS,NYS,NZM,X1,YY1,ZK,             &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               METO%VVPREV(KL)=VAR2
            endif !fluxflg
         else !not awrfflg
            CALL ADVINT(U(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
            CALL ADVINT(U(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
            METO%UUPREV(KL)=(VAR2-VAR1)*TF+VAR1
            CALL ADVINT(V(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
            CALL ADVINT(V(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
            METO%VVPREV(KL)=(VAR2-VAR1)*TF+VAR1
         endif

! JCL:   interpolate air density to position
         CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
         CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
         METO%DENSPREV(KL)=(DENS2-DENS1)*TF+DENS1

! CHG:(03/17/2004) interpolate temperature to position (need for dry air density)
         CALL ADVINT(T(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP1)
         CALL ADVINT(T(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP2)
         METO%TEMPPREV(KL)=(TEMP2-TEMP1)*TF+TEMP1

! CHG:(03/17/2004) interpolate rel. humidity to position (need for dry air density)
         CALL ADVINT(Q(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR1)
         CALL ADVINT(Q(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR2)
         METO%RHFRPREV(KL)=(RHFR2-RHFR1)*TF+RHFR1

! JCL:(4/7/02)should NOT have spatial interpolation of mass violation b/c mass violation
!         was calculated for a gridcell, so should keep same value
! JCL:(4/3/02)interpolate mass violation (fraction of mass/min) profile to position
!         CALL ADVINT(DMASS(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,VAR1)
!         CALL ADVINT(DMASS(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,VAR2)
!         So decided to use ADVINT instead and feed it REAL values that are integers
         XX=DINT(X1)
         YY=DINT(Y1)
         ZZ=DINT(ZK)
         CALL ADVINT(DMASS(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZZ,              &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADVINT(DMASS(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZZ,              &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
! JCL:   temporal interpolation
         METO%DMASSPREV(KL)=(VAR2-VAR1)*TF+VAR1
! JCL:(4/3/02)the density change term is extremely small, so not worry
!        also need the density change term (kg/m3/min) as part of mass budget!
         METO%DMASSPREV(KL)=METO%DMASSPREV(KL)+                         &
     &                     ((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
!        take into account the sign associated with direction in TIME
         IF(BACK)METO%DMASSPREV(KL)=-1.0*METO%DMASSPREV(KL)
      END DO

!CCCCCCCCCCCCCCCRAMS specialCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG(09/11/03) Do following steps for RAMS to "interpolate"
! Allways use K2 for fluxes (later time, they contain averages for period K1-K2
! Use nearest time for scalars
      ELSE IF(RAMSFLG)THEN
! CHG In Rams met files fluxes are staggered 1/2 Dx aggainst scalars
! location of (i,j)=(1,1) is scalar location (center of grid cell)
! all dummies should have been skipped in the rams2arl job file
! so assume i=2.3, means need scalars at 2, and fluxes at 1.8
! so assume i=2.7, means need scalars at 3, and fluxes at 2.2
                   !staggered
        XX1=X1-0.5
                   !staggered
        YY1=Y1-0.5
                     !DNINT for scalars; for x-fluxes use XX1, not XX
        XX=DNINT(X1)
                     !DNINT for scalars; for y-fluxes use YX1, not YY
        YY=DNINT(Y1)
!        WRITE(*,*)'advpnt: XX1,YY1',XX1,YY1
        IF(X1.GE.DBLE(NXS))WRITE(*,*)'advpnt: X1 off grid, >NXS',XX1
        IF(Y1.GE.DBLE(NYS))WRITE(*,*)'advpnt: Y1 off grid, >NYS',YY1
        IF(XX1.LE.1.0)WRITE(*,*)'advpnt: XX1 off grid, <1.0',XX1
        IF(YY1.LE.1.0)WRITE(*,*)'advpnt: YY1 off grid, <1.0',YY1
!        IF(X1.GE.DBLE(NXS-1))WRITE(*,*)'advpnt: X1 off grid, >NXS',XX1
!        IF(Y1.GE.DBLE(NYS-1))WRITE(*,*)'advpnt: Y1 off grid, >NYS',YY1
!        IF(XX1.LE.2.0)WRITE(*,*)'advpnt: XX1 off grid, <1.0',XX1
!        IF(YY1.LE.2.0)WRITE(*,*)'advpnt: YY1 off grid, <1.0',YY1
!        WRITE(*,*)'X1,Y1',X1,Y1
        IF(TF.GE.0.5)THEN
          KNEAR=K2
          KFAR=K1
        ELSE
          KNEAR=K1
          KFAR=K2
        END IF
        IF(BACK)KLATE=K1
        IF(.NOT.BACK)KLATE=K2
        DO KL=1,NLVL
!        set vertical interpolation point to index position
          ZK=DBLE(KL)
          CALL ADVINT(TLRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,         &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
          CALL ADVINT(SIGWRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,       &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)

! DMM: 11/9/2004 Negative values should not be present in TLPREV; it causes pardsp to hang
        IF(VAR1.LT.0.0) THEN
            METO%TLPREV(KL)=0.0
        ELSE
            METO%TLPREV(KL)=VAR1
          ENDIF
          METO%SIGWPREV(KL)=VAR2

          CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
          CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
          METO%DENSPREV(KL)=(DENS2-DENS1)*TF+DENS1

! CHG:(03/17/2004) interpolate temperature to position (need for dry air density)
          CALL ADVINT(T(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP1)
          CALL ADVINT(T(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP2)
          METO%TEMPPREV(KL)=(TEMP2-TEMP1)*TF+TEMP1

! CHG:(03/17/2004) interpolate rel. humidity to position (need for dry air density)
          CALL ADVINT(Q(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR1)
          CALL ADVINT(Q(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR2)
          METO%RHFRPREV(KL)=(RHFR2-RHFR1)*TF+RHFR1

          CALL ADVINT(U(:,:,:,KLATE),NXS,NYS,NZM,XX1,YY,ZK,             &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          METO%UUPREV(KL)=VAR2/DENS2
          CALL ADVINT(V(:,:,:,KLATE),NXS,NYS,NZM,XX,YY1,ZK,             &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          METO%VVPREV(KL)=VAR2/DENS2

          CALL ADVINT(DMASS(:,:,:,KLATE),NXS,NYS,NZM,XX,YY,ZK,          &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
          METO%DMASSPREV(KL)=VAR1
! JCL:(4/3/02)the density change term is extremely small, so not worry
!        also need the density change term (kg/m3/min) as part of mass budget!
!        WRITE(45,*)'advpnt,divu:',METO%DMASSPREV(KL)/(DENS2+DENS1),KL
          METO%DMASSPREV(KL)=METO%DMASSPREV(KL)*2.0/(DENS2+DENS1)       &
     &                     +((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
          IF(BACK)METO%DMASSPREV(KL)=-1.0*METO%DMASSPREV(KL)
!        WRITE(45,*)'advpnt,drodt:',
!     &    ((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM,KL
!        WRITE(45,*)'advpnt, sum:',
!     &    METO%DMASSPREV(KL),KL
!        WRITE(45,*)'K1, K2:',K1,K2
        END DO
      END IF
!        WRITE(45,*)'NXS,NYS:',NXS,NYS
!        WRITE(45,*)'sub grid low left c:',METO%LX1,METO%LY1,xx,yy
!        WRITE(45,*)'new','UUpr',METO%UUPREV
!        WRITE(45,*)'new','UUpr',METO%UUPREV
!        WRITE(*,*)'new','DENSpr',METO%DENSPREV
!        WRITE(*,*)'new','TLprev',METO%TLPREV
!        WRITE(*,*)'new','SIGWprev',METO%SIGWPREV
!        WRITE(45,*)'new','VVpr',METO%VVPREV
!        WRITE(45,*)'new','DMpr',METO%DMASSPREV
! CHG(09/11/03) END Rams special
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! JCL:(4/28/00)interpolate horizontally to get mixed-layer height
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
! JCL:(10/30/02)no longer conduct horizontal spatial interpolation, since want ZML to be consistent w/ turbulence profiles
      XX=DNINT(X1)
      YY=DNINT(Y1)
      ZZZ=1.0
      CALL ADVINTZML(ZML(:,:,K1:K1),NXS,NYS,XX,YY,ZZZ,                   &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADVINTZML(ZML(:,:,K2:K2),NXS,NYS,XX,YY,ZZZ,                   &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      ZMLK1=VAR1
      ZMLK2=VAR2
      METO%ZMLPREV=(VAR2-VAR1)*TF+VAR1

! CHG:(12/04/01)get NEAREST gridpoints value for ZLOC (limit of convection
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
      ZZZ=1.0
      CALL ADVINTZLOC(ZLOC(:,:,K1:K1),NXS,NYS,X1,Y1,ZZZ,                 &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADVINTZLOC(ZLOC(:,:,K2:K2),NXS,NYS,X1,Y1,ZZZ,                 &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      IF(TF.LE.0.5)METO%ZLOCPREV=VAR1
      IF(TF.GT.0.5)METO%ZLOCPREV=VAR2


!CCCCCCCCCCC EXTRACT TEMPORALLY-INTERPOLATED SIGMAW & TL ***BEFORE*** ADVECTION
! JCL:(10/31/02) New way to interpolate vertical profiles of sigmaw & TL between analysis times by scaling
!         the two analysis profiles with the temporally interpolated mixed-layer height
! 1) Find nearest model level to move Zi to
      XX=DNINT(X1)
      YY=DNINT(Y1)
! JCL:(10/31/02) extract ground height with new interpolation routine that deals with 2D arrays
      CALL ADVINT2D(ZT,NXS,NYS,XX,YY,                                   &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      ZGRD=VAR1
! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
!    calculate altitudes at different sigma levels
! CHG(09/10/03) correct transformation between sigma and agl
!     ZLVLS(1)=(1.0-ZSG(1))*(ZMDL-ZGRD)
!      ZLVLS(1)=(1.0-ZSG(1))*ZMDL*ZMDL/(ZMDL-ZGRD)
      ZLVLS(1)=(1.0-ZSG(1))*(ZMDL-ZGRD)
!     retrieve height of each model level [m AGL]
!     ZMIN=(1.0-ZSG(2))*(ZMDL-ZGRD)
!      ZMIN=(1.0-ZSG(2))*ZMDL*ZMDL/(ZMDL-ZGRD)
      ZMIN=(1.0-ZSG(2))*(ZMDL-ZGRD)
      ZLVLS(2)=ZMIN
      ZMIXNEW=ZMIN
      DO KK=3,NLVL
! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
!        ZAGL=(1.0-ZSG(KK))*(ZMDL-ZGRD)
!         ZAGL=(1.0-ZSG(KK))*ZMDL*ZMDL/(ZMDL-ZGRD)
         ZAGL=(1.0-ZSG(KK))*(ZMDL-ZGRD)
         ZLVLS(KK)=ZAGL
         IF(DABS(ZAGL-METO%ZMLPREV).LT.DABS(ZMIXNEW-METO%ZMLPREV))THEN
            ZMIXNEW=ZAGL
         END IF
      END DO
! JCL:(10/31/02) min mixed layer height set to second model level (~75 m)
      METO%ZMLPREV=DMAX1(ZMIXNEW,ZMIN)

! 2)  Find CLOSEST vertical level based on coordinate that's SCALED BY zi
! CHG(09/11/03) don't do this for RAMS
      IF(.NOT.RAMSFLG)THEN
         ZMIX=METO%ZMLPREV
         KK1=1
         KK2=1
         DO KK=1,NLVL
           ! altitude of current timestep (scaled by local zi)
           ZNOWSC=ZLVLS(KK)/ZMIX
           IF (KK1 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK1)/ZMLK1-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK1+1)/ZMLK1-ZNOWSC))
                 KK1=KK1+1
                 IF (KK1 == NLVL) EXIT
              END DO
           END IF
           IF (KK2 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK2)/ZMLK2-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK2+1)/ZMLK2-ZNOWSC))
                 KK2=KK2+1
                 IF (KK2 == NLVL) EXIT
              END DO
           END IF
           METO%SIGWPREV(KK)=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
           METO%TLPREV(KK)=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
         END DO
      ELSE
         ! CHG(09/11/03) assign vertical index here for RAMS, not in ADVIEC
         ! Do this consistent with other runs: ZNDX can be < 1 for below 1st level
         ! For RAMS the 1st level is a flux level...
         ! Here use DREC()%HEIGHT rather than internal heights ZSG(), since need to compare
         ! heights in magl; note there is a shift in indices: HEIGHT(2)~ZSG(1); see hymodelc
         KGNOW = MAX(KG,KGOLD)
         KZ=1
         Z1Z=ZMDL*(1.0-DMIN1(DBLE(1.0),Z1))
         DO WHILE(DREC(KGNOW)%HEIGHT(KZ).LT.Z1Z)
           KZ=KZ+1
         END DO
         KZ=KZ-1
         IF (KZ.GT.0) FRAC=(Z1Z-DREC(KGNOW)%HEIGHT(KZ))/                   &
                              (DREC(KGNOW)%HEIGHT(KZ+1)-DREC(KGNOW)%HEIGHT(KZ))
         IF(KZ.EQ.0)FRAC=Z1Z/DREC(KGNOW)%HEIGHT(KZ+1)
         METO%ZNDX=DBLE(KZ)+FRAC-1.0
      END IF

! JCL:(3/1/01)add WWOUT to output mean vertical velocity [sigma/min]
! JCL:add GD, DUDZ, DVDZ to argument list--used by CHG for box size
!     compute new position
! JCL:(8/13/01)add grdht(METO%ZTER) & 1st sigma level to calculate AGL accurately--to enable accurate interpolation
! JCL:(09/01/03)pass on wind error flag, UPRIME/UERR and VPRIME/VERR from previous timestep
! JCL:(09/01/03) pass on random seed 'RSEED' for modeling transport error as stochastic process
! JCL:(09/01/03) DXERR and DYERR are the horizontal displacements resulting from transport error
! CHG(09/11/03) Pass on RAMSFLG
! CHG(09/11/03) Pass on Density fields for RAMS
! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
      CALL ADVIEC                                                       &
     &   (U(:,:,:,K1), V(:,:,:,K1), W(:,:,:,K1),                        &
     &    U(:,:,:,K2), V(:,:,:,K2), W(:,:,:,K2),                        &
     &    D(:,:,:,K1), D(:,:,:,K2),                                     &
     &    NXS,NYS,NZM,NLVL,DM,JET,ZMDL,METO%ZTER,ZSG(1),                &
     &    X1,Y1,Z1,X2,Y2,Z2,METO%ZNDX,DT,BACK,GX,GY,DUDZ,DVDZ,WWOUT,    &
     &    awrfflg, fluxflg, zsg, &
     &    RAMSFLG,ECMFLG,GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DEAD)
!    :    WINDERRTF,UERRPREV,VERRPREV,RSEED,DXERR,DYERR,RAMSFLG)

      if (dead) then
        WRITE(45,*)'advpnt: particle off grid',X1,Y1
        is_off_grid = .TRUE.
        RETURN
      end if
!     save advection distance as a wind speed (grid pts / min)
      UBAR=DMAX1(ABS(X2-X1),ABS(Y2-Y1))/ABS(DT)
!     map position back to meteo grid
      XP=X2+METO%LX1-1
      YP=Y2+METO%LY1-1
      ZP=DMIN1(DBLE(1.0),Z2)

! JCL: add TL & SIGW to arguments of ADVMET
! JCL:(5/9/01)also added arrays of horizontal velocity (U,V) to arguments
! CHG:(11/20/01)added conv. precip. flag to arguments
! CHG:(11/20/01)added conv. precip. rates (RC) to arguments
! CHG:(11/20/01)added tot. cloud and radiation flag to arguments
! CHG:(11/20/01)added tot. cloud cover and radiation flux to arguments
! CHG:(11/20/01)added surf. energy fluxes (LF,SF) to arguments
! CHG:(11/20/01)added low cloud cover flag to arguments
! CHG:(11/20/01)added low cloud cover to arguments
! JCL:(4/3/02)added mass violation grids (DMASS) to arguments
! CHG:(22/01/03)added soil moisture flag to arguments
! CHG:(22/01/03)added soil moisture to arguments
! CHG:(9/24/02)added LAT&LON to arguments, for DSWF
!      meteo variables interpolated to last advection point
      JTIME=JET+NINT(DT)

! CHG(9/24/02) get updated lat and lon
      IF(KG.NE.0)THEN
! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
         IF(GRID(KG)%LATLON)THEN
           CALL GBL2LL(KG,XP,YP,TLAT,TLON)
         ELSE
           CALL CXY2LL(GRID(KG)%GBASE,XP,YP,TLAT,TLON, GRID(KG)%proj)
                !IF(GRID(KG)%LATLON)THEN
         END IF
             !IF(KG.NE.0)THEN
      END IF
      IF(KG.EQ.0)THEN
! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
         IF(GRID(KGOLD)%LATLON)THEN
           CALL GBL2LL(KGOLD,XP,YP,TLAT,TLON)
         ELSE
           CALL CXY2LL(GRID(KGOLD)%GBASE,XP,YP,TLAT,TLON, GRID(KGOLD)%proj)
                !IF(GRID(KGOLD)%LATLON)THEN
         END IF
             !IF(KG.EQ.0)THEN
      END IF

! JCL:
! CHG(09/16/03)pass on RAMSFLG
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
      CALL ADVMET(BACK,VMIX,CDEP,RDEP,TRAJ,X2,Y2,TLON,TLAT,JTIME,DM,      &
     &   KVEL,DREC(KGRID)%ACYCLE,NXS,NYS,NLVL,FHOUR,IFHR,K1,K2,           &
     &   GX,GY,Z0,LU,ZT,DREC(KGRID)%CFLG,DREC(KGRID)%TCLF,                &
     &   DREC(KGRID)%LCLF,DREC(KGRID)%RADF,DREC(KGRID)%SLMF,              &
     &   DREC(KGRID)%TFLG,                                                &
     &   T(:,:,:,K1),Q(:,:,:,K1),P(:,:,:,K1),D(:,:,:,K1),                 &
     &   X(:,:,:,K1),H(:,:,:,K1),RT(:,:,K1),UF(:,:,K1),VF(:,:,K1),        &
     &   U(:,:,:,K1),V(:,:,:,K1),DMASS(:,:,:,K1),                         &
     &   T(:,:,:,K2),Q(:,:,:,K2),P(:,:,:,K2),D(:,:,:,K2),                 &
     &   X(:,:,:,K2),H(:,:,:,K2),RT(:,:,K2),UF(:,:,K2),VF(:,:,K2),        &
     &   U(:,:,:,K2),V(:,:,:,K2),DMASS(:,:,:,K2),                         &
     &   TL(:,:,:,K1),TL(:,:,:,K2),SIGW(:,:,:,K1),SIGW(:,:,:,K2),         &
     &   RC(:,:,K1),RC(:,:,K2),                                           &
     &   SF(:,:,K1),SF(:,:,K2),LF(:,:,K1),LF(:,:,K2),                     &
     &   TC(:,:,K1),TC(:,:,K2),LC(:,:,K1),LC(:,:,K2),                     &
     &   SW(:,:,K1),SW(:,:,K2),SM(:,:,K1),SM(:,:,K2),                     &
     &   T0(:,:,K1),T0(:,:,K2),                                           &
     &   RAMSFLG,ECMFLG,GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY, &
     &   awrfflg, fluxflg)

! CHG(09/23/03) call ADVMETGRELL to prepare profiles of convective fluxes
      IF (RAMSFLG .OR. GRID(KGRID)%model_id(1:3) == 'ECX' .or. (awrfflg .and. (deepflg .or. shallflg))) then
         IF (BACK) THEN
! CHG(09/23/03) for back, use K1 (absolute later time)
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
            CALL ADVMETGRELL (X2,Y2,CFXUP1(:,:,:,K1),CFXUP2(:,:,:,K1),            &
                              CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),DFXUP2(:,:,:,K1), &
                              EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),DFXDN1(:,:,:,K1), &
                              EFXDN1(:,:,:,K1),TKEN(:,:,:,K1),NXS,NYS,NLVL,       &
                              GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY)
         ELSE
! CHG(09/23/03) for forward, use K2 (absolute later time)
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
            CALL ADVMETGRELL (X2,Y2,CFXUP1(:,:,:,K2),CFXUP2(:,:,:,K2),            &
                              CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),DFXUP2(:,:,:,K2), &
                              EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),DFXDN1(:,:,:,K2), &
                              EFXDN1(:,:,:,K2),TKEN(:,:,:,K2),NXS,NYS,NLVL,       &
                              GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY)
         END IF
      END IF

! JCL:(4/28/00)interpolate horizontally to get mixed-layer height
!     interpolation factor for time after advection
      TF=DMOD(DBLE(JTIME),DM)/DM
      IF(BACK)THEN
         TF=1.0-TF
           !make sure after advection use K2 when TF=0.0 (JTIME updated)
      ELSE
         IF(TF.EQ.0.0)TF=1.0
      END IF

! JCL:(10/31/02)extract profiles of TL & SIGMAW at the advected particle location for the two analysis times
! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
      DO KL=1,NLVL
!        set vertical interpolation point to index position
         ZK=KL

! JCL:   interpolate Lagrangian timescale to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
         XX=DNINT(X2)
         YY=DNINT(Y2)
         CALL ADVINT(TL(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADVINT(TL(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         TLK1(KL)=VAR1
         TLK2(KL)=VAR2

! JCL:   interpolate std dev of vertical velocity to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
         XX=DNINT(X2)
         YY=DNINT(Y2)
         CALL ADVINT(SIGW(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,               &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADVINT(SIGW(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,               &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         SIGWK1(KL)=VAR1
         SIGWK2(KL)=VAR2

      END DO
             !of IF(.NOT.RAMSFLG)THEN
      END IF

! JCL:spatial interpolation to get local ML ht with 'ADVINTZML'
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
! JCL:(10/30/02)no longer conduct horizontal spatial interpolation, since want ZML to be consistent w/ turbulence profiles
      XX=DNINT(X2)
      YY=DNINT(Y2)
      ZZZ=1.0
      CALL ADVINTZML(ZML(:,:,K1:K1),NXS,NYS,XX,YY,ZZZ,                   &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADVINTZML(ZML(:,:,K2:K2),NXS,NYS,XX,YY,ZZZ,                   &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      ZMLK1=VAR1
      ZMLK2=VAR2
      METO%ZMLNEXT=(VAR2-VAR1)*TF+VAR1

! CHG:(12/04/01)get nearest gridpoints value for ZLOC (limit of convection
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
      ZZZ=1.0
      CALL ADVINTZLOC(ZLOC(:,:,K1:K1),NXS,NYS,X2,Y2,ZZZ,                 &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADVINTZLOC(ZLOC(:,:,K2:K2),NXS,NYS,X2,Y2,ZZZ,                 &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      IF(TF.LE.0.5)METO%ZLOCNEXT=VAR1
      IF(TF.GT.0.5)METO%ZLOCNEXT=VAR2

! JCL:(6/28/02)interpolate air density to particle position
      ZTMP=METO%ZNDX
      IF(ZTMP.LT.1.0)ZTMP=1.0
                           !see comments in advmet CHG(09/16/03)
      ZTMPR=DINT(ZTMP+1.0)
      ZTMPR=DMIN1(DMAX1(DBLE(1.0),ZTMPR),DBLE(NLVL))
! CHG(09/16/03) for rams no spatial interpolation in density
      IF (.NOT. RAMSFLG) THEN
        CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,X2,Y2,ZTMP,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
        CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,X2,Y2,ZTMP,                 &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
      ELSE
        CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZTMPR,                &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
        CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZTMPR,                &
     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
      END IF
      METO%DENSLOCAL=(DENS2-DENS1)*TF+DENS1

!CCCCCCCCCCC EXTRACT TEMPORALLY-INTERPOLATED SIGMAW & TL ***AFTER*** ADVECTION
! JCL:(10/31/02) New way to interpolate vertical profiles of sigmaw & TL between analysis times by scaling
!         the two analysis profiles with the temporally interpolated mixed-layer height
! 1) Find nearest model level to move Zi to
      XX=DNINT(X2)
      YY=DNINT(Y2)
! JCL:(10/31/02) extract ground height with new interpolation routine that deals with 2D arrays
      CALL ADVINT2D(ZT,NXS,NYS,XX,YY,                                   &
     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      ZGRD=VAR1
! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
!    calculate altitudes at different sigma levels
! CHG(09/10/03) correct transformation between sigma and agl
!     ZLVLS(1)=(1.0-ZSG(1))*(ZMDL-ZGRD)
!      ZLVLS(1)=(1.0-ZSG(1))*ZMDL*ZMDL/(ZMDL-ZGRD)
      ZLVLS(1)=(1.0-ZSG(1))*(ZMDL-ZGRD)
!     retrieve height of each model level [m AGL]
!     ZMIN=(1.0-ZSG(2))*(ZMDL-ZGRD)
!      ZMIN=(1.0-ZSG(2))*ZMDL*ZMDL/(ZMDL-ZGRD)
      ZMIN=(1.0-ZSG(2))*(ZMDL-ZGRD)
      ZLVLS(2)=ZMIN
      ZMIXNEW=ZMIN
      DO KK=3,NLVL
! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
!        ZAGL=(1.0-ZSG(KK))*(ZMDL-ZGRD)
!         ZAGL=(1.0-ZSG(KK))*ZMDL*ZMDL/(ZMDL-ZGRD)
         ZAGL=(1.0-ZSG(KK))*(ZMDL-ZGRD)
         ZLVLS(KK)=ZAGL
         IF(DABS(ZAGL-METO%ZMLPREV).LT.DABS(ZMIXNEW-METO%ZMLPREV))THEN
            ZMIXNEW=ZAGL
         END IF
      END DO
! JCL:(10/31/02) min mixed layer height set to second model level (~75 m)
      METO%ZMLNEXT=DMAX1(ZMIXNEW,ZMIN)

! 2)  Find CLOSEST vertical level based on coordinate that's SCALED BY zi
! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
         ZMIX=METO%ZMLNEXT
         KK1=1
         KK2=1
         DO KK=1,NLVL
           ! altitude of current timestep (scaled by local zi)
           ZNOWSC=ZLVLS(KK)/ZMIX
           IF (KK1 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK1)/ZMLK1-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK1+1)/ZMLK1-ZNOWSC))
                 KK1=KK1+1
                 IF (KK1 == NLVL) EXIT
              END DO
           END IF
           IF (KK2 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK2)/ZMLK2-ZNOWSC).GT.                      &
                           ABS(ZLVLS(KK2+1)/ZMLK2-ZNOWSC))
                 KK2=KK2+1
                 IF (KK2 == NLVL) EXIT
              END DO
           END IF
           METO%SIGW(KK)=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
           METO%TL(KK)=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
         END DO
           !do this for RAMS
      ELSE
                     !DNINT for scalars
        XX=DNINT(X2)
                     !DNINT for scalars
        YY=DNINT(Y2)
!        IF(TF.GE.0.5)THEN
!          KNEAR=K2
!        ELSE
!          KNEAR=K1
!        END IF
! CHG(3/3/2004) Now get time averaged sigw and tl from RAMS, stored at later time
        KNEAR=K2
        DO KK=1,NLVL
          ZK=DBLE(KK)
          CALL ADVINT(TL(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,             &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
          CALL ADVINT(SIGWRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,       &
     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          METO%TL(KK)=VAR1
          METO%SIGW(KK)=VAR2
        END DO
            !of IF(.NOT.RAMSFLG)THEN ... ELSE
      ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! JCL:
!      WRITE(*,*)"ADVPNT:ZMIX1 = ",METO%ZMLPREV,"   ZMIX2=",METO%ZMLNEXT
!      DO KL=1,NLVL
!         WRITE(*,*)KL,METO%SIGWPREV(KL),METO%SIGW(KL),METO%TLPREV(KL),
!     &                                                 METO%TL(KL)
!      END DO

      RETURN
      END
