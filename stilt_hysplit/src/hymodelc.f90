!$$$  SUBROUTINE DOCUMENTATION BLOCK
!
! MAIN PROGRAM: HYMODELC     MAIN HYSPLIT PROGRAM FOR AIR CONCENTRATIONS
!   PRGMMR: DRAXLER          ORG: R/E/AR      DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   IS THE MAIN PROGRAM FOR THE TRANSPORT AND DISPERSION MODEL HYSPLIT. GIVEN
!   A TIME, SOURCE LOCATION, AND RELEASE AMOUNT, THE MODEL COMPUTES AIR
!   CONCENTRATIONS OVER PRESPECIFIED SAMPLING PERIODS, ON PRESELECTED
!   LATITUDE-LONGITUDE GRIDS AT VARIOUS HEIGHTS ABOVE GROUND.  REQUIRED
!   METEOROLOGICAL DATA ARE PRESUMED TO HAVE ALREADY BEEN PREPARED FOR MODEL
!   INPUT. SEE ROUTINES AVN2ARL FOR A DISCUSSION OF METEO DATA PREPARATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 May 1998 (RRD)
!                 17 Aug 1998 (RRD) - added delt & isot option to namelist
!                                   - particle dump and initialization
!                 22 Dec 1998 (RRD) - included area source emission file
!                 12 Feb 1999 (RRD) - pardump option fixed for grid change
!                 04 Mar 1999 (RRD) - subgrid optimiztion
!                 20 Apr 1999 (RRD) - terrain compression factor
!                 14 May 1999 (RRD) - pgrd=0 test prior to advpnt
!                 09 Jun 1999 (RRD) - changed structure of PARDUMP file
!                 15 Jun 1999 (RRD) - added tratio to namelist
!                 26 Aug 1999 (RRD) - pass parameters for vertical grid
!                                   - puff to particle conversion
!                 08 Nov 1999 (RRD) - correction to PARDUMP read
!
! USAGE:  HYMODELC
!
!   INPUT PARAMETERS:
!     PROMPTED ON STANDARD INPUT UNLESS FILE NAMED "CONTROL" EXISTS
!   INPUT FILES:
!     CONTROL (unit 30) - CONTAINS REQUIRED MODEL INPUT PARAMETERS
!     METEOROLOGICAL DATA (units 10,12) - AS NAMED IN CONTROL FILE
!     ROUGLEN.ASC (unit 40) - OPTIONAL AERODYNAMIC ROUGHNESS LENGTH
!     LANDUSE.ASC (unit 41) - OPTIONAL LAND-USE
!   OUTPUT FILES:
!     BINARY CONCENTRATIONS (unit 21,22) as named in CONTROL file
!     MESSAGE (unit 30) for diagnostic ouputs
!     STARTUP (unit 31) for a copy of standard input
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: hymodelc.f90,v 1.40 2009-03-05 08:26:22 gerbig Exp $
!
!$$$

PROGRAM HYMODELC

!     array sizes
      USE module_defgrid                              ! meteorology grid and file
      USE module_defconc                              ! pollutants and concentration grid
      USE module_defmeto                              ! meteo variables returned after advection
      USE module_defspot                              ! multiple source structure

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER, PARAMETER :: dp=KIND(1d0)

! JCL:(11/03/03) store profile of wind errors
      TYPE WINDPROF
         REAL*8 U(NZM)
         REAL*8 V(NZM)
      END TYPE WINDPROF
      TYPE (WINDPROF), ALLOCATABLE :: UVERR(:)
      REAL*8 UERR(NZM),VERR(NZM),UERR2(NZM),VERR2(NZM)
      REAL*8 UUERR_T(NZM),VVERR_T(NZM)
      REAL*8 SIGUVERR,TLUVERR,ZCORLEN
      REAL*8 SIGZIERR,TLZIERR,HORCORZILEN,RELZIERR,ZIERR(MAXPAR)

!     pollutant particle/puff array elements
      REAL*8 XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR),SIGH(MAXPAR),       &
     &     SIGV(MAXPAR),SIGX(MAXPAR),MASS(MAXDIM,MAXPAR)
      INTEGER PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)

!     index variable to contain position sorted elements
      INTEGER NSORT(MAXPAR)

!     master concentration array (x,y,z,species,grids)
      REAL*8 CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!     deposition accumulation variable
      REAL*8 DEPT(MAXTYP)

! JCL:(5/18/00)array to hold WPRIME/SIGMAW for each particle exiting PARDSP
      REAL*8 WWPREV(MAXPAR)
! JCL:(09/01/03)array to store UPRIME/UERR & VPRIME/VERR for each particle--to store wind err fluctuations
!      REAL*8 UERRPREV(MAXPAR),VERRPREV(MAXPAR)
! CHG(09/24/03)array to store rel. area coverage (updraft/downdraft) at particle location
      REAL*8 AREAPRU(MAXPAR),AREAPRD(MAXPAR)
! CHG(09/24/03)array to store cloud index (1=updraft, 2=env., 3=downdraft)
      INTEGER ICNDX(MAXPAR)

! JCL:(6/29/00)array to hold total SAMPTT (calculated in PARDSP) for each particle before results are written out
      REAL*8 SAMPTTCUM(MAXPAR)
      REAL*8 FOOTCUM(MAXPAR)

! CHG(24/06/04)array to store vertical displacement due to convection
      REAL*8 ZFXCUM(MAXPAR)

! CHG:(12/05/01)array to hold 'CONVDUR' (calc. in ADVPNT, reset to 0 after conv. redistribution was done in PARDSP)
      INTEGER CONVDUR(MAXPAR), CONVTMP
      logical is_off_grid !replaces alternate return mechanism

! JCL:(4/3/02)weighting of particles due to mass violation (calculated in ADVPNT and tallied in PARDSP) for each particle
      REAL*8 DMASSWT(MAXPAR)
! JCL:(4/3/02)the total mass violation experienced by particle in PARDSP during each timestep
      REAL*8 DMASSTT

!     model sigma levels, mass diagnostic information
      REAL*8 ZSG(NZM),      ZMASS(NZM)
      real*8 zprofm(nzm) !flux-level zsg expressed as height AGL (without terrain compression)
      real*8 zsigw, zsigwdn !flux-level zsg used in zprofm computation

      ! index record extended header,  file information
      CHARACTER(10000) :: HEADER
      CHARACTER(256)   :: FNAME

! JCL:(05/12/2004) ALLOCATE probably didn't work originally b/c didn't have proper f90 compiler
! JCL:new variables to correct 'mistake' in original code
! JCL:set aside 100 locations for XARG & YARG, since dynamic allocation of
!     arrays doesn't seem to work. XARG&YARG are the X&Y coords of starting
!     locations, so 100 is a pretty large number--unlikely that will exceed
!      DIMENSION XARG(100),YARG(100)
      REAL*8, ALLOCATABLE :: XARG(:), YARG(:)

! JCL:variables to denote particles' positions that would be written out
      REAL*8 XOUT,YOUT,ZOUT

! JCL:standard deviation of vertical & horizontal velocity that is calculated in PARDSP
      REAL*8 SIGMAW
      REAL*8 SIGMAU

! JCL:amount of time [min.] that particle 'sees' the ground--calculated in PARDSP:
!      SAMPTT=(# particle touchdowns)*(deltat); deltat is timestep of subloop in PARDSP
      REAL*8 SAMPTT
      REAL*8 FOOT

! CHG:used to store Lagrangian timescale--calculated in PARDSP
      REAL*8 SAMPTT2

! JCL:variables used to store wind shear [(m/s)/m]
      REAL*8 DUDZ, DVDZ

      INTEGER SMIN

!     integration direction, deposition flags, file test, area source
      LOGICAL BACK, CDEP, RDEP, SDEP, FTEST, QFILE, PRESET

! JCL:logical flag to specify whether data should be written out to PARTICLE.DAT or not
      LOGICAL WRITEOUT

! JCL:flag to tell whether particle has 'seen' vegetation on ground
      INTEGER SEEVEG

! JCL:random seed for random number generator--would always be
!     '5' if the flag 'RANDOM' is selected to be '0'; otherwise
!     the random seed itself is 'randomized' using the current time & RANDOM_NUMBER
      INTEGER              :: RSEED, dtt(8), i
      REAL                 :: harvest

! JCL:the number of seconds since midnight--used to generate random seed
      REAL*8 SECS

! JCL:declre the variables read in as NAMELIST
      INTEGER INITD,KHMAX,NUMPAR,QCYCLE,KRND,ISOT,NDUMP,RANDOM
      REAL*8 FRMR,DELT

! JCL:(2/28/2004) enable user to choose DYNAMICALLY the output variables
      CHARACTER(LEN=4), DIMENSION(100) :: VARSIWANT
      INTEGER                        :: IVMAX

! JCL:(9/16/02)flag to say whether PBL height was prescribed or not
      INTEGER ZICONTROLTF
! JCL:(9/16/02)vector to store the scaling factors to prescribe PBL height
      REAL*8 ZIPRESC(150)

! JCL:(09/01/03)flag specifying whether to include wind errors as Markov process
      INTEGER WINDERRTF

! JCL:counter to keep track of # of [min] since last time output written to PARTICLE.DAT
      REAL*8 COUNTOUT

!     assumed maximum wind speed in km/min for calculation of subgrid size
      REAL(dp) :: UMAXI=1.2

! JCL:(2/13/2001)variable to store Wbar (vertical velocity)
      REAL*8 WWOUT

! JCL:(4/27/01)counter to keep track of how many particles have left model area
      INTEGER COUNTNPAROUT

      REAL(dp) :: TLON=-HUGE(TLON), TLAT=-HUGE(TLON)

! CHG(09/16/03):add flag specifying whether data from RAMS or not
      LOGICAL RAMSFLG, awrfflg, ECMFLG
      INTEGER :: dummy

!     parameters that define vertical grid polynomial
      COMMON /ZZTOKK/ AA,BB,CC

!     special simulation setup parameters read from namelist
! JCL:add an extra variable into SETUP--'RANDOM', which tells random number generator
!     whether to have a different random sequence each time model is run; if set FALSE,
!     then generates same random sequence each time=>useful for debugging purposes
! JCL:add 'TLFRAC'--determines the internal timestep in PARDSP as fraction of TL
! JCL:add 'OUTDT'--the interval to write out data to 'PARTICLE.DAT' [min]
! CHG:add 'NTURB'--flag for no turbulence (to calc. mean wind trajectories)
! JCL:(6/29/00)'VEGHT' is ht [m] below which a particle would be counted as 'seeing' grd vegetation
! JCL:(4/27/01)'OUTFRAC' determines fraction of particle # that leave model area that if exceeded,
!       then model stops
! JCL:(9/16/02) ZICONTROLTF specifies whether prescribe mixed-layer height (=1) or not (=0)
! CHG:(9/17/02) add 'ICONVECT'--flag for convection
! JCL:(09/01/03) WINDERRTF specifies whether to include wind errs as Markov process (=1) or not (=0)
! JCL:(02/28/2004) additional vars to dynamically determine output variables
      NAMELIST/SETUP/INITD,KHMAX,NUMPAR,QCYCLE,FRME,FRMR,KRND,          &
     &               DELT,ISOT,NDUMP,RANDOM,TRATIO,AA,BB,CC,TLFRAC,     &
     &               OUTDT,NTURB,VEGHT,OUTFRAC,ZICONTROLTF,ICONVECT,    &
     &               WINDERRTF,IVMAX,UMAXI,VARSIWANT

!=>default model configuration options

!     constants that control puff rounding and merging
!     standard fractional constants used each hour
!                    horizontal, vertical, temporal
      DATA FRHS,FRVS,FRTS /1.00,   0.5,   0.1/
!     enhanced rounding parameters used at KRND intervals
!                    horizontal, vertical, temporal
      DATA FRHE,FRVE,FRTE /1.75,   2.0,   0.2/


!---------------------------------------------------------------------------------------------------
!     title line to standard output
      WRITE(*,*)'HYSPLIT4 (Nov 99) - Initialization'

!=>required for NCEP operational implementation
!     CALL W3TAGB('HYMODELC',1998,0238,0068,'R/E/AR ')

!     unique input file information added as suffix to CONTROL.{FNAME}
!     standard output files get same extension

! jcl:the following two lines were commented out in PC source file
!    so go ahead and do so as well
!      CALL GETARG(1,FNAME)
!      KLEN=INDEX(FNAME,' ')-1
! jcl:'KLEN=0' was found in PC source file, so add it as well
       KLEN=0

!     title line to standard output
      WRITE(*,*)'HYSPLIT4 (Nov 99) - Initialization'


! JCL:initalize output flag to FALSE
      WRITEOUT=.FALSE.

!=>additional model setup parameters read from namelist

!     initial distribution: 0-3dPart 1-Gauss 2-TopHat 3-1/0 4-2/0
      INITD=4
!     maximum duration in hours of any particle/puff
      KHMAX=9999
!     number of puffs in simulation or particles to release
      NUMPAR=100
!     optional cycling of emissions (hours) anatex=60
      QCYCLE = 0
!     mass rounding fraction for enhanced merging
      FRME=0.10
!     mass removal fraction during enhanced merging
      FRMR=0.0
!     enhanced merging interval (hours)
      KRND=6
!     integration time step (0 - autoset; >0 - constant minutes)
      DELT=0.0
!     isotropic turbulence option (0 - off; 1 - on)
      ISOT=0
!     initialize and dump particles to/from file
!     0-none    1-read/write    2-read only     3-write only
      NDUMP=0
!     advection stability ratio
      TRATIO=0.75
! JCL:default random seed is 5
      RSEED=5
! JCL:default of RANDOM is TRUE  (0-FALSE; 1-TRUE)
      RANDOM=1
! JCL:fraction of TL (Lagrangian timescale) to set as timestep in dispersion subroutine
      TLFRAC=0.1
! JCL:default value of OUTDT is 0.0--means output is written EVERY timestep
      OUTDT=0.0
! JCL:ht [fraction of PBL ht or m] below which a particle's time spent is tallied;
!        If <=1.0, then specifies fraction of PBL ht
      VEGHT=0.5
! CHG:default value of NTURB is FALSE means TURBULENCE is on
      NTURB=0
! JCL:(4/27/01)default value of OUTFRAC--specifies fraction of particle number over which
!        model would stop if fraction of particles that leaves model area exceeds this number
!        e.g., if OUTFRAC=0.9, once over 90% of particles leave model area, then model stops
      OUTFRAC=0.9
! JCL:(9/16/02) initialize PBL height prescription flag to FALSE (0-FALSE; 1-TRUE)
      ZICONTROLTF=0
! CHG:(9/17/02)default value of ICONVECT is FALSE ~means CONVECTION is off
      ICONVECT=0
      ICONVECTRAMS=0
! JCL:(09/01/03) wind err flag to FALSE (0-FALSE; 1-TRUE)
! JCL:(09/01/03) wind err flag to FALSE (0-FALSE; 1-TRUE, 2: zi error, 3: zi and wind error)
      WINDERRTF=0
! JCL:(02/28/2004) default output variables to be written to 'PARTICLE.DAT'
!     number of variables to output
      IVMAX=5
!     4-lettered codes for the different variables
      VARSIWANT(1:5)=(/'TIME','INDX','LONG','LATI','ZAGL'/)
!     vertical grid polynomial defaults Z= aa*k^2 + bb*k + cc
!     mechanism: if set in SETUP.CFG - no further change, otherwise set in sbr runset depending
!        on meteo model ID
      AA = -HUGE(AA)
      BB = AA
      CC = AA

      INQUIRE(FILE='SETUP.CFG',EXIST=FTEST)
      IF(FTEST)THEN
         OPEN(30,FILE='SETUP.CFG',STATUS='OLD',ACTION='READ')
         READ(30,NML=SETUP)
         CLOSE(30)
         WRITE(*,*)'NOTICE: using SETUP.CFG namelist file'
      END IF

! JCL:(02/28/2004) make sure that not select too many output variables
      IF(IVMAX.GT.100) STOP 'Too many output vars (IVMAX) in SETUP.CFG!'

!     optional conversion of puffs to particles flag
      PRESET=.FALSE.
      IF(INITD.GT.4)THEN
         PRESET=.TRUE.
         INITD=INITD-4
         WRITE(*,*)'NOTICE: puff to particle conversion enabled'
      END IF

!=>check some initial model parameters for consistency

!     MAss removal fraction set to zero for particle simulations
      IF(INITD.EQ.0)FRMR=0.0

!     simulation number should be less than dimensioned value
      NUMPAR=MIN0(NUMPAR,MAXPAR)
!     check compilation dimensions for special applications
      IF (MAXDIM > MAXTYP) WRITE(*,*) 'Warning module_defconc/module_defsize: maxtyp >= maxdim'

!     special limit for horizontal splitting and vertical particle
!     numpar gives the particles in the initial release, need to leave
!     room in the array for further horizontal splitting
      NUMSPL=NUMPAR
      IF(INITD.EQ.3.OR.INITD.EQ.4)NUMSPL=MAXPAR/2


!=>standard meteorological model initialization

!     set the basic simulation control parameters (unit 30)
      CALL RUNSET(NLOC,NHRS,NGRD,NLVL,KSFC,ZSG,ZMDL,SFCL,               &
                  BACK,IUNIT,KVEL,KLEN,FNAME, AA, BB, CC)

!     initialize information about each meteo data grid
      CALL METINI(HEADER,NGRD,SPOT(1)%OLAT,SPOT(1)%IBYR,                &
     &   SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,BACK,KVEL)

! CHG(09/18/03) moved from after setting JET=MC0
!     set default meteorological grid number for starting particles
      KG=1
!     set current active meteorological grid ( when no particles = 0 )
      KGRID=0

! CHG(09/09/03) Need to use excact flux levels from RAMS
! Assume no change in vertical grid between nests
! i.e. set NLVL (number of internal levels),
! ZSG (pseudo sigma internal coordinate),
! ZMDL (to corresponding value from RAMS at (NZRAMS-1)),
! KSFC (vertical index for top of sfc layer, was 2, now 1), SFCL (=HGT(KSFC)),
      IF(GRID(KG)%MODEL_ID.EQ.'RAMS')THEN
! CHG&JCL (03/10/2004) set convection flag for RAMS
!     so that other routines are not affected
        ICONVECTRAMS=ICONVECT
! CHG&JCL (03/10/2004) no excess convection for RAMS
        IF(ICONVECT.EQ.1)ICONVECT=0
        RAMSFLG = .TRUE.
        NLVL=GRID(KG)%NZ-2
        ZMDL=DREC(KG)%HEIGHT(NLVL+1)
        ZSG(1:NLVL)=1.0-DREC(KG)%HEIGHT(2:(NLVL+1))/ZMDL
        KSFC=1
        SFCL=DREC(KG)%HEIGHT(KSFC+1)
      ELSE
        RAMSFLG = .FALSE.
      ENDIF
      awrfflg = GRID(KG)%MODEL_ID(2:4) .EQ. 'WRF'
      ECMFLG  = GRID(KG)%MODEL_ID(1:2) == 'EC'
      if (awrfflg .OR. ECMFLG) then
         if (iconvect .eq. 1) then !turn on cgrell, turn off extreme conv. (as for RAMS)
            iconvectrams = 1
            iconvect = 0
         elseif (iconvect .eq. -1) then !turn off cgrell, turn on extreme conv.
            iconvectrams = 0
            iconvect = 1
         endif
      end if

! Optionally read in AGL heights of model layers from external file:
      if (.not. RAMSFLG) then
         call read_zsg ('ZSG_LEVS.IN',45,zmdl,zsg,nzm,nlvl,ksfc,sfcl,aa,bb,cc)
      endif

! CHG(09/16/03) moved from runset
! check starting height limit against model top, using ZSG and ZMDL
      max_hgt = zmdl * (1.-zsg(nlvl))
!!$      IF(RAMSFLG)THEN
      DO N=1,NLOC
!         starting height limit
!!$        IF(SPOT(N)%OLVL.GT.ZMDL)THEN
        IF(SPOT(N)%OLVL.GT.max_hgt)THEN
          WRITE(*,*)'ERROR hymodelc: Start height above mdl domain'
          WRITE(*,*)'   Source - ',N,'    Height - ',SPOT(N)%OLVL
!!$          WRITE(*,*)'   Model domain - ',INT(DREC(KG)%HEIGHT(NLVL+1))
          WRITE(*,*)'   Model domain - ',max_hgt
          WRITE(*,*)'   Model top ht - ',ZMDL
          STOP
        END IF
      END DO
!!$      ELSE
!!$        DO N=1,NLOC
!!$!          starting height limit
!!$           IF(SPOT(N)%OLVL.GT.(AA*NLVL*NLVL+BB*NLVL+CC))THEN
!!$              WRITE(*,*)'ERROR runset: Start height above mdl domain'
!!$              WRITE(*,*)'  Source - ',N,'    Height - ',SPOT(N)%OLVL
!!$              WRITE(*,*)'  Model domain - ',INT(AA*NLVL*NLVL+BB*NLVL+CC)
!!$              WRITE(*,*)'  Model top ht - ',ZMDL
!!$              STOP
!!$           END IF
!!$        END DO
!!$      END IF

!     convert starting time to accumulated minutes
      CALL TM2MIN(SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,               &
     &   SPOT(1)%IBHR,0,MC0)


!     elapsed time
      JET=MC0

!     default initial meteo grid (km) for time-step calculation
!     the first grid defined always has the finest resolution
      CGSIZE=GRID(KG)%SIZE

! JCL:
      OPEN(45,FILE='JCLmessage',FORM='FORMATTED')

! JCL:open file to write particle information
      OPEN(50,FILE='PARTICLE.DAT',FORM='FORMATTED')

! get a random number RSEED to initialize RAN3, if RANDOM=1
! JCL(050803): implemented new random seed initialization suggested by Stefan Koerner
      IF(RANDOM.EQ.1)THEN
        CALL RANDOM_SEED                       ! initializes the build-in random number generator
        CALL DATE_AND_TIME (values=dtt)
        DO i=1,MAX(dtt(8),1)                   ! dtt(8) contains the milliseconds of the wall-clock
         CALL RANDOM_NUMBER (harvest)          ! the build-in random number generator for uniform
        END DO                                 !    distribution between 0 and 1
        RSEED = ABS(NINT((harvest-0.5)*1e4))
      END IF


! JCL:
      WRITE(45,*) 'RSEED: ',RSEED

! JCL:(3/16/01) Instead of setting initial WWPREV (WPRIME normalized by stddev) to 0.0,
!               set it to a GAUSSIAN DISTRIBUTION
! JCL:(6/29/00) Initialized the cumulative SAMPTT to 0.0
! JCL:(4/3/02)  Initialized the weights of particles from mass violation [fraction of gridcell]
! CHG(09/24/03) Initialize rel. area coverage up/downdraft at particle location
! CHG(09/24/03) Initialize cloud index to environment
      DO I=1,MAXPAR
         WWPREV(I)=GASDEV(RSEED,DBLE(1.0))
         SAMPTTCUM(I)=0.0
         FOOTCUM(I)=0.0
         DMASSWT(I)=1.0
         AREAPRU(I)=0.0
         AREAPRD(I)=0.0
                    !start in environment
         ICNDX(I)=2
! CHG(24/06/04) add ZFXCUM as vertical displacement due to convective flux
!  (deep or shallow, up or downdraft) along trajectory [m]
         ZFXCUM(I)=0.0
      END DO

!     set starting location
      DO N=1,NLOC
! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
         IF(GRID(KG)%LATLON)THEN
            CALL GBL2XY(KG,SPOT(N)%OLAT,SPOT(N)%OLON,                   &
     &                     SPOT(N)%XP,SPOT(N)%YP)
         ELSE
            CALL CLL2XY(GRID(KG)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,       &
     &                  SPOT(N)%XP, SPOT(N)%YP, GRID(KG)%proj)
         END IF


         IF(SPOT(N)%XP.LT.1.OR.SPOT(N)%XP.GT.GRID(KG)%NX.OR.            &
     &      SPOT(N)%YP.LT.1.OR.SPOT(N)%YP.GT.GRID(KG)%NY)THEN
            WRITE(*,*)'ERROR main: source point off grid'
            WRITE(*,*)'Position:',SPOT(N)%OLAT,SPOT(N)%OLON
            WRITE(*,*)'Grid Loc:',SPOT(N)%XP,  SPOT(N)%YP
            STOP
         END IF
!        default without terrain adjustment

         SPOT(N)%ZP=(ZMDL-SPOT(N)%OLVL)/ZMDL
      END DO

! JCL:write starting location to PARTICLE.DAT
!     use any old variable--ZOUT in this case, to get formatting correct
!      ZOUT=0.0
!      WRITE(50,*)ZOUT,ZOUT,SPOT(1)%OLAT,SPOT(1)%OLON,SPOT(1)%OLVL
! JCL:(11/6/02) not needed anymore
      WRITE(50,*)SPOT(1)%OLAT,SPOT(1)%OLON,SPOT(1)%OLVL

! JCL:(9/16/02) read in prescribed scaling factors for mixed-layer height
      IF(ZICONTROLTF.EQ.1)THEN
      INQUIRE(FILE='ZICONTROL',EXIST=FTEST)
      IF(FTEST)THEN
         OPEN(46,FILE='ZICONTROL',ACTION='READ')
         READ(46,*)NHRSZI
!         WRITE(45,*)'NHRSZI=',NHRSZI
         DO KNZI=1,NHRSZI
            READ(46,*)ZIPRESC(KNZI)
!            WRITE(45,*)KNZI,' ZIPRESC=',ZIPRESC(KNZI)
         END DO
         WRITE(45,'(a,i4,a,g15.6)') ' ZI prescribed from ZICONTROL: nhrszi= ', &
              & nhrszi,' zipresc(1)= ',ZIPRESC(1)
         CLOSE(46)
      ELSE
         WRITE(*,*)'File ZICONTROL not found!'
         STOP
      END IF
      END IF


!=>concentration calculation initialization

!     emission (source term) initializaton
      CALL EMSSET(NUMTYP,IUNIT,                                         &
     &   SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR)
! JCL: added BACK as argument to be able to assign concentration
!      sampling starting & stopping times accurately
!     sampling grid initialization (update CGSIZE if required)


      CALL CONSET(ZMDL,NUMGRD,IUNIT,CGSIZE,                             &
     &   SPOT(1)%OLAT,SPOT(1)%OLON,                                     &
     &   SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,BACK)
!     open concentration output files (units 21,...), zero variables
      CALL CONINI(NLOC,NUMGRD,NUMTYP,CSUM)
!     set the deposition parameters for each pollutant type
      CALL DEPSET(NUMTYP,IUNIT,CDEP,RDEP,SDEP)

!=>diagnostic file initialization

!     close all initialization files
      IF(IUNIT.EQ.5)THEN
         CLOSE(31)
      ELSE
         CLOSE(30)
      END IF

!     define unique name for message file if required
      IF(KLEN.LE.0)THEN
         OPEN(30,FILE='MESSAGE')
      ELSE
         OPEN(30,FILE='MESSAGE.'//FNAME(1:KLEN))
      END IF
      WRITE(*,*)'Calculation Started ... please be patient'

!=>diagnostic message on status of deposition flags

      IF(DIRT(1)%DOGAS)WRITE(30,*)'Gas pollutant  : ',DIRT(1)%DOGAS
      IF(DIRT(1)%DODRY)WRITE(30,*)'Dry deposition : ',DIRT(1)%DODRY
      IF(DIRT(1)%DOWET)WRITE(30,*)'Wet deposition : ',DIRT(1)%DOWET
      IF(DIRT(1)%DORES)WRITE(30,*)'Dynamic dry    : ',DIRT(1)%DORES
      IF(DIRT(1)%DOGRV)WRITE(30,*)'Grav settling  : ',DIRT(1)%DOGRV
      IF(DIRT(1)%DOSUS)WRITE(30,*)'Resuspension   : ',DIRT(1)%DOSUS
      IF(DIRT(1)%DORAD)WRITE(30,*)'Radioactive    : ',DIRT(1)%DORAD

!=>initialize gridded emissions array, use with emsgrd

!     enables gridded emissions from an ascii file (emission.txt)
!     used in conjunction with call emsgrd
      CALL EMSDAT(NLOC,NUMTYP,QFILE)

!     initial value of particle/puff counter and time step velocity
      KPM=0
!     initial wind speed max (1.2 km/min = 20 m/s) for dt
      UMAX=UMAXI

!=>optional initialization from previous simulation dump file

      IF(NDUMP.GT.0)THEN
!        all combinations require file to be opened
         INQUIRE(FILE='PARDUMP',EXIST=FTEST)
         OPEN(32,FILE='PARDUMP',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

!        update array for read options (ndump = 1 [rw], 2 [r], 3 [w] )
         IF(FTEST.AND.(NDUMP.EQ.1.OR.NDUMP.EQ.2))THEN

            READ (32,ERR=100,END=100)KPM,NUMPOL
            DO J=1,KPM
               READ (32) (MASS(I,J),I=1,NUMPOL)
               READ (32) TLAT,TLON,ZPOS(J),SIGH(J),SIGV(J),SIGX(J)
               READ (32) PAGE(J),HDWP(J),PTYP(J),PGRD(J),NSORT(J)

! CHG(11/05/03) also read DMASSWT(KP)
               READ (32)DMASSWT(J)

! CHG&JCL (03/10/2004) avoid using PGRD(KG) as index when eqal to 0 (i.e. off grid)
!              reset meteo grid index back to default start grid
               IF(PGRD(J).NE.0)THEN
!              convert positions from lat/lon to x,y after reading dump
! CHG (06/28/2004) changed from GRID(PGRD(KG)) to GRID(PGRD(J)) (was a bug fixed in Vers 4.5)
!                  not missing then reset meteo grid index back to default start grid
                 PGRD(J)=KG
! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
                 IF(GRID(PGRD(J))%LATLON)THEN
                   CALL GBL2XY(PGRD(J),TLAT,TLON,XPOS(J),YPOS(J))
                 ELSE
!                  convert positions to lat/lon before dump
                   CALL CLL2XY(GRID(PGRD(J))%GBASE,TLAT,TLON,           &
     &                        XPOS(J),YPOS(J), GRID(PGRD(J))%proj)
                         !IF(GRID(PGRD(J))%LATLON)THEN
                 END IF

               ELSE
                 XPOS(J)=0.0
                 YPOS(J)=0.0
               END IF

!              convert height to sigma ASSUMING TERRAIN HEIGHT = 0!!!!
! jcl          simply the reverse transformation of output;as long as consistent, doesn't matter
               ZPOS(J)=(ZMDL-ZPOS(J))/ZMDL

!              flag any puffs that might be offgrid, previous
!              may have created particles on a different grid
               IF(XPOS(J).LE.1.OR.XPOS(J).GE.GRID(KG)%NX.OR.            &
     &            YPOS(J).LE.1.OR.YPOS(J).GE.GRID(KG)%NY)               &
     &            PGRD(J)=0
            END DO

!           delete any particles that may not be on new grid
            XXX=0.0
            CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,         &
     &           PAGE,PTYP,PGRD,NSORT,KHMAX,XXX)
            WRITE(30,*)'NOTICE: hymodelc particle initialization',KPM
            GOTO 110

  100       WRITE(30,*)'WARNING: hymodelc - error in reading particle'
            WRITE(30,*)'         initialization file PARDUMP'
            KPM=0

  110       CONTINUE
            REWIND(32)

!        end read-write test
         END IF
!     end ndump exist test
      END IF

!=>starting sigma height initialization

!     convert to gp/min
      VMAX=UMAX/GRID(MAX0(KG,KGRID))%SIZE

! JCL:(05/12/2004) ALLOCATE probably didn't work originally b/c didn't have proper f90 compiler
! jcl:instead of giving ADVRNG SPOT.XP & SPOT.YP as arguments,
!     give it XARG & YARG as arguments
!     comment out ALLOCATE b/c in UNIX version, has already allocated 100 elements
      ALLOCATE( XARG(NLOC))
      ALLOCATE( YARG(NLOC))
! jcl:add loop to would let XARG & YARG take on values of
!     SPOT.XP and SPOT.YP
! CHG(09/08/03) replace loop
      XARG(1:NLOC)=SPOT(1:NLOC)%XP
      YARG(1:NLOC)=SPOT(1:NLOC)%YP

!     redefine subgrid parameters
!     CALL ADVRNG(MAX0(KG,KGRID),VMAX,NLOC,SPOT.XP,SPOT.YP)
      CALL ADVRNG(MAX0(KG,KGRID),VMAX,NLOC,XARG,YARG)

      DT=10.0
      DO 150 N=1,NLOC
         XPOS(N)=SPOT(N)%XP
         YPOS(N)=SPOT(N)%YP
         ZPOS(N)=SPOT(N)%ZP
         KMET=MAX0(KG,KGRID)

! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! CHG(09/18/03) pass on RAMSFLG
!        dummy call to advpnt to return terrain height
         CALL ADVPNT                                                    &
     &      (BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,ISOT,                   &
     &       XPOS(N),YPOS(N),ZPOS(N),JET,DT,KMET,KGRID,                 &
     &       NGRD,ZSG,NLVL,ZMDL,KVEL,UBAR,IFHR,DUDZ,DVDZ,               &
     &       WWOUT,ZICONTROLTF,NHRSZI,ZIPRESC,CONVDUR(1),ICONVECT,      &
     &       RAMSFLG,ECMFLG,is_off_grid)
         if (is_off_grid) goto 150
!        terrain adjusted sigma

         SPOT(N)%ZP=1.0-SPOT(N)%OLVL/(ZMDL-METO%ZTER)
  150 END DO

! JCL:initialize counter: # of [min] since last time output written to PARTICLE.DAT
      COUNTOUT=0.0

!=>main loop over number of simulation hours

      DO KH=1,NHRS

!     time step minutes function of gridsize and max wind speed
!     such that (uvel dt <= TRATIO * delx) and dt(max)=60, where uvel
!     is defined as grid points per minute and delx is the grid
!     spacing in km
      IF (UMAX <= 0d0) UMAX=UMAXI
      MAXDT=MIN0(60,DREC(MAX0(KG,KGRID))%DELTA)
      DT=AMAX0(1, MIN0(MAXDT, NINT(TRATIO*CGSIZE/UMAX)))
      DO WHILE (MOD(MAXDT,INT(DT)).NE.0.AND.INT(DT).GT.1)
         DT=DT-1.0
      END DO

!     namelist over-ride for fixed time step when delt<>0
      IF(DELT.GT.0.0)DT=DELT

! CHG(09/23/03) don't use fixed time step for RAMS
      IF(RAMSFLG.AND.DELT.GT.0.0)THEN
        WRITE(*,*)"ERROR: don't use fixed time step for RAMS"
        STOP
      END IF

! JCL:
! the following line was found in the trajectory model HYMODELT
! but was missing from the distributed source code for HYMODELC
! this line is needed to run dispersion model BACKWARDS
      IF (BACK) THEN
        DT=-DT
      END IF

!     convert to gp/min before call to advrng
      UMAX=UMAX/GRID(MAX0(KG,KGRID))%SIZE

!     redefine subgrid parameters
      IF(KPM.LE.0)THEN
!        no particels or first time initialize with source
!          ALLOCATE( XARG(NLOC))
!          ALLOCATE( YARG(NLOC))
!     jcl:add loop to would let XARG & YARG take on values of
!     SPOT.XP and SPOT.YP
! CHG(09/08/03) replaced loop
          XARG(1:NLOC)=SPOT(1:NLOC)%XP
          YARG(1:NLOC)=SPOT(1:NLOC)%YP

!         CALL ADVRNG(MAX0(KG,KGRID),UMAX,NLOC,SPOT%XP,SPOT%YP)
!     jcl:change erroneous statement above to:
          CALL ADVRNG(MAX0(KG,KGRID),UMAX,NLOC,XARG,YARG)
      ELSE
!        subsequent time use mean particle position to define grid
         CALL ADVRNG(MAX0(KG,KGRID),UMAX,KPM,XPOS,YPOS)
      END IF

!     reset time step variables
      UMAX=0.0
      CGSIZE=GRID(1)%SIZE


!CCCCCCCCCCCCCCCCCCCC
!=>sub-loop for number of time steps per hour
! JCL:IABS part takes care of DT<0
      NSTEP=IABS(60/INT(DT))

      DO KS=1,NSTEP

!        update grid on which new particles start if required
         IF(KGRID.GT.KG)THEN
            KG=KGRID
! CHG(09/22/03) update RAMSFLG
            RAMSFLG=GRID(KG)%MODEL_ID.EQ.'RAMS'
            ECMFLG =GRID(KG)%MODEL_ID(1:2) == 'EC'
            AWRFFLG=GRID(KG)%MODEL_ID(2:4).EQ.'WRF'

!           remap starting positions to new grid
            DO N=1,NLOC

! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
              IF(GRID(KG)%LATLON)THEN
                 CALL GBL2XY(KG,SPOT(N)%OLAT,SPOT(N)%OLON,              &
     &                          SPOT(N)%XP, SPOT(N)%YP)
              ELSE
                 CALL CLL2XY(GRID(KG)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,  &
     &                    SPOT(N)%XP,SPOT(N)%YP, GRID(KG)%proj)
                     !IF(GRID(KG)%LATLON)THEN
              END IF

                     !DO N=1,NLOC
            END DO
         END IF

!=>emissions each time step are from all grid cells or from a point

         IF(QFILE)THEN
! JCL:  add BACK as argument
!           start any new particles from grid cells (from emsdat)
            CALL EMSGRD(NUMTYP,KPM,INITD,DT,JET,KG,                     &
     &         NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,      &
     &         PTYP,PGRD,REAL(QCYCLE,dp),NUMPAR,TBAR,BACK)

         ELSE
! JCL:  add BACK as argument
!           start any new particles from a point (from emsset)
! CHG: HERE IS THE PROBLEM!
            CALL EMSPNT(NLOC,NUMTYP,KPM,INITD,DT,JET,KG,                &
     &         NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,      &
     &         PTYP,PGRD,REAL(QCYCLE,dp),NUMPAR,BACK)
         END IF

!        diagnostic variables initialization
         TMASS=0.0
! CHG(09/08/03) replaced loop
         ZMASS(1:NZM)=0.0

!        optional position sortting for small meteo subgrids
!        this routine requires further optimization
         IF(NGRD.EQ.1.AND.NXM.LT.GRID(KG)%NX.AND.NYM.LT.GRID(KG)%NY)    &
     &      CALL ADVSRT(KPM,METO%LX1,METO%LY1,XPOS,YPOS,NSORT)

!        computation of solar angle only required for gaseous dry deposition
!        or certain specialized chemistry applications, for simplicity
!        assume that same angle applies over the entire grid at this time
         IF(RDEP)                                                       &
     &      CALL SUNANG                                                 &
     &      (SPOT(1)%IBYR,JET,SPOT(1)%OLAT,SPOT(1)%OLON,EA,SEA)

!        simplified resuspension uses last advection point meteorological
!        data for computation - deposition grid not yet linked to meteo
         IF(SDEP)                                                       &
     &      CALL DEPSUS(INITD,KG,NUMGRD,NUMTYP,DT,KPM,CSUM,MASS,        &
     &      XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,NSORT)

! JCL:   increment # of [min] since last time output written to PARTICLE.DAT
         COUNTOUT=COUNTOUT+ABS(INT(DT))

! JCL:   write out results to PARTICLE.DAT only in specified intervals
!             and OUTDT=0 means that results would be output EACH timestep
         WRITEOUT=(OUTDT.EQ.0.0).OR.(COUNTOUT.GE.OUTDT)
! JCL:   ensures results always written for LAST timestep
         WRITEOUT=WRITEOUT.OR.(KH.EQ.NHRS.AND.KS.EQ.NSTEP)

! JCL:(4/27/01)initialize counter for number of particles that left grid area
         COUNTNPAROUT=0

! CHG:(9/17/02) initialize temporary convtmp variable to 0
         CONVTMP=0

!CCCCCCCCC
!=>loop through all particles/puffs

!         WRITE(*,*)'(2) KPM=',KPM,KS,WINDERRTF

!-----------------------------------------------------------
! JCL:(11/03/03) store profile of wind errors
         IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3))THEN
            WRITE(*,*)'ENTERED INITIALIZATION!!!!'

!           read in statistics of wind errors
            INQUIRE(FILE='WINDERR',EXIST=FTEST)
            IF(FTEST)THEN
              OPEN(46,FILE='WINDERR')
              !stdev of U&V errs [m/s]
              READ(46,*)SIGUVERR
              !correlation Lagr timescale of U&V errs [min]
              READ(46,*)TLUVERR
              !correlation lengthscale of U&V errs in vertical dir [m]
              READ(46,*)ZCORLEN
              !correlation lengthscale of U&V errs in hor. dir [km]
              READ(46,*)HORCORLEN
              CLOSE(46)
            ELSE
              WRITE(*,*)'============File "WINDERR" not found=========='
              STOP
                   !IF(FTEST)THEN
            END IF

!           initialize error profiles with appropriate vertical correlation
           !KPM is not initialized at this point, so not use KPM!!!
           ALLOCATE(UVERR(KPM))
!            ALLOCATE(UVERR(MAXPAR))
            DO KPP=1,KPM
               UVERR(KPP)%U(1)=GASDEV(RSEED,DBLE(1.0))
               UVERR(KPP)%V(1)=GASDEV(RSEED,DBLE(1.0))
            DO KZZ=2,NZM
               DELTAZ=(ZMDL-METO%ZTER)*(ZSG(KZZ-1)-ZSG(KZZ))
               RAUTO=DEXP(-1.0*DELTAZ/ZCORLEN)
               UU=GASDEV(RSEED,DBLE(1.0))
               UVERR(KPP)%U(KZZ)=RAUTO*UVERR(KPP)%U(KZZ-1)+             &
     &                           DSQRT(1.0-RAUTO*RAUTO)*UU
               VV=GASDEV(RSEED,DBLE(1.0))
               UVERR(KPP)%V(KZZ)=RAUTO*UVERR(KPP)%V(KZZ-1)+             &
     &                           DSQRT(1.0-RAUTO*RAUTO)*VV
!               WRITE(*,*)"initial:",KZZ,UVERR(KPP)%U(KZZ)
                   !DO KZZ=2,NZM
            END DO
                   !DO KPP=1,KPM
            END DO

!***********test******************
!            U0=GASDEV(RSEED,DBLE(2.5))
!            RAUTO=DEXP(-30.0/240.0)
!            DO JJJJ=1,1000
!              UU=RAUTO*U0+DSQRT(1.0-RAUTO*RAUTO)*GASDEV(RSEED,DBLE(2.5))
!               U0=UU
!               WRITE(45,*)UU
!            END DO
!***********test******************
                   !IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3))THEN
         END IF
!-----------------------------------------------------------
! CHG:(27/04/06) initialize zi errors
         IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3))THEN
            WRITE(*,*)'ENTERED ZI ERROR INITIALIZATION!!!!'

!           read in statistics of ZI errors
            INQUIRE(FILE='ZIERR',EXIST=FTEST)
            IF(FTEST)THEN
              OPEN(47,FILE='ZIERR')
              !stdev of ZI errs [%]
              READ(47,*)SIGZIERR
              !correlation Lagr timescale of U&V errs [min]
              READ(47,*)TLZIERR
              !correlation lengthscale of U&V errs in hor. dir [km]
              READ(47,*)HORCORZILEN
              CLOSE(47)
              WRITE(*,*)'SIGZIERR:   ',SIGZIERR
              WRITE(*,*)'TLZIERR:    ',TLZIERR
              WRITE(*,*)'HORCORZILEN:',HORCORZILEN
            ELSE
              WRITE(*,*)'============File "ZIERR" not found=========='
              STOP
                   !IF(FTEST)THEN
            END IF
            DO KPP=1,KPM
               ZIERR(KPP)=GASDEV(RSEED,SIGZIERR)
            END DO
                   !IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3))THEN
         END IF
!-----------------------------------------------------------

         DO KPT=1,KPM

!           current computation grid index number
            KP=NSORT(KPT)

!           skip terminated particles
            IF(PGRD(KP).EQ.0)GO TO 200

! JCL:(5/9/01)store horizontal position before mean advection--to pass onto PARDSP
            XPOSPREV=XPOS(KP)
            YPOSPREV=YPOS(KP)

! JCL:(3/1/01)add WWOUT to output the mean vertical velocity [sigma/min]
! JCL:      add DUDZ&DVDZ as output arguments--CHG uses them to
!                   parameterize trajectory box size
!           advects a single point for one time step
! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
! CHG:(12/05/01) add 'CONVDUR' (calc. in ADVPNT, reset to 0 after conv. redistribution was done in PARDSP)
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! CHG(09/18/03) pass on RAMSFLG
            CALL ADVPNT                                                 &
     &       (BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,ISOT,                  &
     &       XPOS(KP),YPOS(KP),ZPOS(KP),JET,DT,PGRD(KP),KGRID,          &
     &       NGRD,ZSG,NLVL,ZMDL,KVEL,UBAR,IFHR,DUDZ,DVDZ,               &
     &       WWOUT,ZICONTROLTF,NHRSZI,ZIPRESC,CONVDUR(KP),ICONVECT,     &
     &       RAMSFLG,ECMFLG,is_off_grid)
            if (is_off_grid) goto 200

!CCCCCCCCCCCCCCCCCCCC  CONVECTION USING RAMS FLUXES CCCCCCCCCCCCCCCCCCCC
! CHG(09/23/03) call CGRELL for convective redistribution w/ RAMS fluxes
            IF ((RAMSFLG .OR. ECMFLG .or. awrfflg) .AND. NTURB .EQ. 0 .AND. ICONVECTRAMS .EQ. 1)THEN
                                      !convert to normalized height
              Z1=(1.0-ZPOS(KP))*ZMDL
              if (awrfflg) then
! this uses the same logic as the dzsig computation in prfcom:
                 zsigwdn=1.
                 zsigw=2.*zsg(1)-zsigwdn
                 zprofm(1) = (1.0-zsigw)*zmdl
                 zsigwdn=zsigw
                 do kk = 2,nlvl
                    zsigw=2.*zsg(kk)-zsigwdn
                    zprofm(kk) = (1.-zsigw)*zmdl
                    zsigwdn=zsigw
                 enddo
              elseif (ramsflg) then
                 zprofm(1:nlvl) = DREC(KG)%HEIGHT(2:(NLVL+1))
              ELSE IF (ECMFLG) THEN                         ! arithm. mean of mid-level z values
                 FORALL (kk=1:nlvl-1) zprofm(kk) = 0.5d0*zmdl*(2d0-zsg(kk)-zsg(kk+1))
                 zprofm(nlvl) = zmdl
              else
                 stop 'hymodelc cgrell if-struct: internal logic error'
              endif
! CHG(09/24/03) relative area coverage up/downdraft not yet available from RAMS
! area should be near zero if no flux
! CHG(09/25/03) try to use RAUP1 ~ CFU1(cloud base) / (DENS(CB)*(1+sqrt(TKE*2) ))
!             Get cloud base: first level w/ non-zero CFXUP
              KCB1=1
              KCB2=1
              KDB1=1
              !cloud base deep
              DO WHILE(METO%CFXUP1(KCB1).EQ.0.0.AND.KCB1.LT.NLVL)
                KCB1=KCB1+1
              END DO
              !cloud base shallow
              DO WHILE(METO%CFXUP2(KCB2).EQ.0.0.AND.KCB2.LT.NLVL)
                KCB2=KCB2+1
              END DO
              !downdraft (should allways be to ground)
              DO WHILE(METO%CFXDN1(KDB1).EQ.0.0.AND.KDB1.LT.NLVL)
                KDB1=KDB1+1
              END DO
!             get top of convection
              KCT1=NLVL
              KCT2=NLVL
              KDT1=NLVL
               !cloud top deep
              DO WHILE(METO%CFXUP1(KCT1).EQ.0.0.AND.KCT1.GT.1)
                KCT1=KCT1-1
              END DO
              !cloud top shallow
              DO WHILE(METO%CFXUP2(KCT2).EQ.0.0.AND.KCT2.GT.1)
                KCT2=KCT2-1
              END DO
              !downdraft top
              DO WHILE(METO%CFXDN1(KDT1).EQ.0.0.AND.KDT1.GT.1)
                KDT1=KDT1-1
              END DO
!             Extract TKEN and relative area coverage
              RAUP=0.0
              IF (DREC(KG)%TKEF) THEN
                IF (KCB2.LT.NLVL)RAUP=METO%CFXUP2(KCB2)/                   &
                   (METO%DENS(KCB2)*(1+DSQRT(METO%TKEN(KCB2)*2)))
                IF (KCB1.LT.NLVL)RAUP=METO%CFXUP1(KCB1)/                   &
                   (METO%DENS(KCB1)*(1+DSQRT(METO%TKEN(KCB1)*2)))
              ELSE
                IF (KCB2.LT.NLVL)RAUP=METO%CFXUP2(KCB2)/                   &
                   (METO%DENS(KCB2)*(1+METO%SIGW(KCB2)*2))
                IF (KCB1.LT.NLVL)RAUP=METO%CFXUP1(KCB1)/                   &
                   (METO%DENS(KCB1)*(1+METO%SIGW(KCB1)*2))
              END IF
                                         !reset to environment
              IF(RAUP.EQ.0.0)ICNDX2=2
              !reset area coverage here since not entering CGRELL
              IF(KCT1.LT.KCB1)AREAPRU(KP)=0.0
              !reset area coverage here since not entering CGRELL
              IF(KCT2.LT.KCB2)AREAPRU(KP)=0.0
              !reset area coverage here since not entering CGRELL
              IF(KDT1.LT.KDB1)AREAPRD(KP)=0.0

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!             AREA=METO%GDIS*METO%GDIS
              AREA=METO%GDISX*METO%GDISY
!        IF(KP.EQ.1)WRITE(*,*)'Z before',Z1,ICNDX(KP)
! Only call CGRELL if convection
              IF(RAUP.GT.0.0)THEN
! assign relative area coverage at proper levels
                !initialize with small number
                METO%RAUP1(1:NZM)=1.0E-15
                METO%RAUP2(1:NZM)=1.0E-15
                METO%RADN1(1:NZM)=1.0E-15
                IF(KCT1.GE.KCB1)METO%RAUP1(KCB1:KCT1)=RAUP
                IF(KCT2.GE.KCB2)METO%RAUP2(KCB2:KCT2)=RAUP
                IF(KDT1.GE.KDB1)METO%RADN1(KDB1:KDT1)=RAUP

! diagnostics...
!       IF(KP.EQ.23)THEN
!       DO KKK=1,NLVL
!         WRITE(45,*)'hymodelc',KKK,METO%CFXUP1(KKK),METO%EFXUP1(KKK),
!     &     METO%DFXUP1(KKK),METO%RAUP1(KKK),AREAPRU(KP),Z1,ICNDX(KP),
!     &     DREC(KG)%HEIGHT(1+KKK)
!       END DO
!       ENDIF
!                WRITE(*,*)'hymodelc, before CGRELL',Z1,ICNDX(KP)
                CALL CGRELL(METO%CFXUP1,METO%CFXUP2,METO%CFXDN1,        &
     &           METO%DFXUP1,METO%DFXUP2,METO%EFXUP1,METO%EFXUP2,       &
     &           METO%DFXDN1,METO%EFXDN1,                               &
     &           METO%RAUP1,METO%RAUP2,METO%RADN1,                      &
     &           AREA,                                                  &
     &           AREAPRU(KP),AREAPRD(KP),METO%DENS,NLVL,                &
     &           zprofm(1:nlvl),Z1,Z2,                     &
     &           ICNDX(KP),ICNDX2,dummy,DT,BACK)
                                     !convert back to sigma
                ZPOS(KP)=1.0-Z2/ZMDL
!                WRITE(*,*)'Z after',Z2,ICNDX2
! CHG(24/06/04) add ZFXCUM as vertical displacement due to convective flux
!  (deep or shallow, up or downdraft) along trajectory [m]
                ZFXCUM(KP)=ZFXCUM(KP)+DABS(Z2-Z1)
                ! WRITE(*,'(a,i3,a,f0.1)') 'ZFXCUM(', KP, '): ',ZFXCUM(KP)

              END IF !of IF(RAUP.GT.0.0)

              !assign changed cloud index (moved by CHG(3/4/2004))
              ICNDX(KP)=ICNDX2

            END IF !of IF(RAMSFLG.AND.NTURB.EQ.0)THEN
!CCCCCCCCCCCCCCCCC END CONVECTION USING RAMS FLUXES CCCCCCCCCCCCCCCCCCCC

! CHG:need to set CONVDUR only for first call of ADVPNT during current particle loop
            IF(CONVDUR(KP).LT.CONVTMP)CONVDUR(KP)=CONVTMP
            CONVTMP=CONVDUR(KP)
! CHG(11/29/02) NO CONV. REDIST. WHEN NO CRAIN (OR RAIN IF CRAIN<0, i.e. NA)
!            IF(METO%CRAI.EQ.0.0.OR.
!     &         METO%CRAI.LT.0.0.AND.METO%RAIN.EQ.0.0)CONVDUR(KP)=0
!           skip terminated particles
            IF(PGRD(KP).EQ.0)GO TO 200

!           surface terrain height for this position
            ZSFC=METO%ZTER

!           convert advection u (gp/min) to km/min [minimum=1m/s]
            UMAX=MAX(0.06d0,UMAX,UBAR*GRID(PGRD(KP))%SIZE)

!           increment particle age after advection
            PAGE(KP)=PAGE(KP)+DT

! JCL:      default is .FALSE. (=0)
            SEEVEG=0

!-----------------------------------------------------------------------
! JCL:(11/03/03) store profile of wind errors
      IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        CONVFAC=METO%GDIS/60.0   !conversion factor for [grid/min]=>[m/s]
         CONVFACX=METO%GDISX/60.0
         CONVFACY=METO%GDISY/60.0

!        impose VERTICAL correlation
         UERR(1)=GASDEV(RSEED,SIGUVERR)
         VERR(1)=GASDEV(RSEED,SIGUVERR)
         UERR2(1)=GASDEV(RSEED,SIGUVERR)
         VERR2(1)=GASDEV(RSEED,SIGUVERR)
         DO KZZ=2,NLVL
           DELTAZ=(ZMDL-ZSFC)*(ZSG(KZZ-1)-ZSG(KZZ))
           RAUTO=DEXP(-1.0*DELTAZ/ZCORLEN)
           UU=GASDEV(RSEED,SIGUVERR)
           UERR(KZZ)=RAUTO*UERR(KZZ-1)+DSQRT(1.0-RAUTO*RAUTO)*UU
           UU2=GASDEV(RSEED,SIGUVERR)
           UERR2(KZZ)=RAUTO*UERR2(KZZ-1)+DSQRT(1.0-RAUTO*RAUTO)*UU2
           VV=GASDEV(RSEED,SIGUVERR)
           VERR(KZZ)=RAUTO*VERR(KZZ-1)+DSQRT(1.0-RAUTO*RAUTO)*VV
           VV2=GASDEV(RSEED,SIGUVERR)
           VERR2(KZZ)=RAUTO*VERR2(KZZ-1)+DSQRT(1.0-RAUTO*RAUTO)*VV2
         END DO

!        impose TEMPORAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
         RAUTO=DEXP(-1.0*DABS(DT)/TLUVERR)
         DO KZZ=1,NLVL
           UUERRPREV=UVERR(KP)%U(KZZ)*SIGUVERR
           VVERRPREV=UVERR(KP)%V(KZZ)*SIGUVERR
           !in [grid/min]
           METO%UUPREV(KZZ)=UUERRPREV/CONVFACX+METO%UUPREV(KZZ)
           UUERR_T(KZZ)=RAUTO*UUERRPREV+DSQRT(1.0-RAUTO*RAUTO)*UERR(KZZ)
           !in [grid/min]
           METO%VVPREV(KZZ)=VVERRPREV/CONVFACY+METO%VVPREV(KZZ)
           VVERR_T(KZZ)=RAUTO*VVERRPREV+DSQRT(1.0-RAUTO*RAUTO)*VERR(KZZ)
                  !DO KZZ=1,NLVL
         END DO

!        impose SPATIAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        CFAC=METO%GDIS/1000.0   ![m/grid]=>[km/grid]
                                   ![m/grid]=>[km/grid]
         CFACX=METO%GDISX/1000.0
                                   ![m/grid]=>[km/grid]
         CFACY=METO%GDISY/1000.0
         !distance travelled [km]
         DIST=DSQRT((CFACX*(XPOS(KP)-XPOSPREV))**2+                     &
     &              (CFACY*(YPOS(KP)-YPOSPREV))**2)
         RAUTOH=DEXP(-1.0*DIST/HORCORLEN)
         DO KZZ=1,NLVL
           UUERRNEXT=RAUTOH*UUERR_T(KZZ)+                               &
     &                DSQRT(1.0-RAUTOH*RAUTOH)*UERR2(KZZ)
           !in [grid/min]
           METO%UUNEXT(KZZ)=UUERRNEXT/CONVFACX+METO%UUNEXT(KZZ)

           VVERRNEXT=RAUTOH*VVERR_T(KZZ)+                               &
     &                DSQRT(1.0-RAUTOH*RAUTOH)*VERR2(KZZ)
           !in [grid/min]
           METO%VVNEXT(KZZ)=VVERRNEXT/CONVFACY+METO%VVNEXT(KZZ)
!          store for next timestep
           UVERR(KP)%U(KZZ)=UUERRNEXT/SIGUVERR
           UVERR(KP)%V(KZZ)=VVERRNEXT/SIGUVERR
!           WRITE(*,*)"SPATIAL :",KZZ,UUERR_T(KZZ),
!     &           UUERRNEXT,METO%UUNEXT(KZZ)*CONVFAC,UERR2(KZZ),RAUTOH
                  !DO KZZ=1,NLVL
         END DO

              !IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
      END IF
!-----------------------------------------------------------------------
! CHG:(28/04/06) get zi errors: here relative error in percent
      IF(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3)THEN

!        impose TEMPORAL correlation
         RAUTO=DEXP(-1.0*DABS(DT)/TLZIERR)
         ZIERR(KP)=RAUTO*ZIERR(KP)+DSQRT(1.0-RAUTO*RAUTO)*GASDEV(RSEED,SIGZIERR)

!        impose SPATIAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
                                   ![m/grid]=>[km/grid]
         CFACX=METO%GDISX/1000.0
                                   ![m/grid]=>[km/grid]
         CFACY=METO%GDISY/1000.0
         !distance travelled [km]
         DIST=DSQRT((CFACX*(XPOS(KP)-XPOSPREV))**2+                     &
     &              (CFACY*(YPOS(KP)-YPOSPREV))**2)
         RAUTOH=DEXP(-1.0*DIST/HORCORZILEN)
         ZIERR(KP)=RAUTOH*ZIERR(KP)+DSQRT(1.0-RAUTOH*RAUTOH)*GASDEV(RSEED,SIGZIERR)
              !IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
      END IF
!-----------------------------------------------------------------------

!            WRITE(*,*)'before PARDSP, XPOS, YPOS:',XPOSPREV,YPOSPREV
!CCCCCCCCCCCFollowing comments are for PARDSP
! JCL:      add SEEVEG as argument--gets set to 1 if particle sees veg
!           also add TLFRAC--fraction of TL that fixes internal timestep
!           also add RSEED (random seed) as argument
!           also add SIGMAW,SIGMAU (std dev of vert&hor velocity) as argument
!           SAMPTT is amt of time [min.] that particle 'sees' the ground
!           METO%TL is vertical profile of Lagrangian timescales [s]
!           METO%SIGW is vertical profile of std dev of vert velocity
!           METO%ZMLNEXT is mixed-layer ht [m]
!           METO%DENSNEXT is vertical profile air density [kg/m3]
! JCL:(4/13/2000)
!           METO%TLPREV,METO%SIGWPREV,METO%ZMLPREV,METO%DENSPREV are values from
!                 BEFORE mean advection timestep
! CHG:      SAMPTT2 is variable for output from pardsp, now Tlagrange vert.
!               (different from 'SIGMAW'--the OUTPUT from PARDSP)
! JCL:(4/27/00)Since want vertical profiles to be CONSTANT, not have them as argument
! JCL:(5/2/00) also remove ZSG as argument
! JCL:(5/18/00)add as 1st argument--# of min. since start
! JCL:(5/18/00)WWPREV is array storing NORMALIZED turbulent velocity for each particle
!              WWPREV gets updated in PARDSP
! JCL:(5/18/00)Arguments when have TIME-VARYING vertical profiles of TL & SIGMAW
! CHG:      NTURB true switched off turbulence
! JCL:(6/29/00)'VEGHT' is ht [m] below which a particle would be counted as 'seeing' grd vegetation
!                   VEGHT determines whether SAMPTT is > 0 or not
! JCL:(5/9/01)added horizontal position before mean advection (XPOSPREV & YPOSPREV)
!                  + vertical profiles of horizontal winds (UUPREV,UUNEXT,VVPREV,VVNEXT)
!             -these are used to implement interaction between vertical turbulence and wind shear
! JCL:(4/3/02)added vertical profile of mass violation before & after mean advection
! JCL:(4/3/02)added weighting of particle due to mass violation experienced by particle in PARDSP (DMASSWT)
! CHG:(12/05/01) add 'CONVDUR' (calc. in ADVPNT, reset to 0 after conv. redistribution was done in PARDSP)
! CHG:(12/04/01)add 'ZLOC' limit of convection altitude (m)
! CHG(09/11/03) Pass on RAMSFLG
! CHG:(03/17/2004) pass on rel. humidity and temperature profile to get dry air density column integrated
! CHG:(03/17/2004) also add output variable for specific humidity SPHU and FOOT
! JCL:(07/14/2004) split up GDIS=>GDISX&GDISY for global grids (non-conformal)
            SMIN=JET+INT(DT)-(CONC(1)%START%MACC)
            IF (NTURB == 0) THEN
               CALL PARDSP(SMIN,                                                 &
                  METO%HMIX, METO%GDISX, METO%GDISY, DT, ZMDL, ZSFC, NLVL,       &
                  METO%VMIX, ZSG, XPOS(KP), YPOS(KP), ZPOS(KP), SIGH(KP),        &
                  SIGV(KP), SIGX(KP), HDWP(KP), METO%ZNDX, ISOT, SEEVEG,         &
                  METO%TLPREV, METO%SIGWPREV, METO%ZMLPREV, METO%DENSPREV,       &
                  METO%TL, METO%SIGW, METO%ZMLNEXT, METO%DENS, TLFRAC,           &
                  XPOSPREV, YPOSPREV,                                            &
                  METO%UUPREV, METO%UUNEXT, METO%VVPREV, METO%VVNEXT,            &
                  METO%DMASSPREV, METO%DMASSNEXT, DMASSWT(KP), BACK,             &
                  RSEED, SIGMAW, SIGMAU, SAMPTT, SAMPTT2, WWPREV(KP), KP, VEGHT, &
                  CONVDUR(KP), METO%ZLOCNEXT, RAMSFLG, ECMFLG,                   &
                  METO%TEMP, METO%TEMPPREV, METO%RHFR, METO%RHFRPREV,            &
                  SPHU, FOOT)
            ELSE
               SIGMAW  = 0d0                                ! define output values
               SAMPTT  = 0d0
               SAMPTT2 = 0d0
               FOOT    = 0d0
            END IF
!CCCCCCCCCCC
! JCL:(11/03/03) apply transport error by directly changing wind profile
! JCL:(09/01/03) apply displacement resulting from transport error
!            IF(WINDERRTF.EQ.1)THEN
!               XPOS(KP)=XPOS(KP)+DXERR
!               YPOS(KP)=YPOS(KP)+DYERR
!            END IF

! JCL:(8/28/03)calculate total mass violation experienced by particle (fraction of mass in gridcell) for MEAN trajs
! CHG(10/10/03) different for RAMS
            IF (.NOT. RAMSFLG .AND.  NTURB == 1) THEN
               KZ=MAX0(1, MIN0(NLVL,INT(METO%ZNDX)))
               KT=MIN0(KZ+1,NLVL)
               IF(ZPOS(KP).LT.ZSG(1))THEN
                  DMASS=(METO%DMASSPREV(KT)+METO%DMASSNEXT(KT))/2.0
               ELSE
                  DMASS=(METO%DMASSPREV(1)+METO%DMASSNEXT(1))/2.0
               END IF
               DMASSWT(KP)=DMASSWT(KP)*(1.0+DMASS*DABS(DT))
            END IF
            IF(RAMSFLG.AND.NTURB.EQ.1)THEN
              KT=MAX0(1, MIN0(NLVL,INT(METO%ZNDX)+1))
              DMASS=(METO%DMASSPREV(KT)+METO%DMASSNEXT(KT))/2.0
              DMASSWT(KP)=DMASSWT(KP)*(1.0+DMASS*DABS(DT))
            END IF

! JCL:(6/29/00)increment cumulative SAMPTT
            SAMPTTCUM(KP)=SAMPTTCUM(KP)+SAMPTT
! CHG:(28/04/2006) apply ZI error to FOOT
            IF(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3)THEN
              RELZIERR=(100+ZIERR(KP))/100
              RELZIERR=MAX(RELZIERR,0d0)
              FOOT=FOOT*RELZIERR
!         increment cumulative FOOT with ZI error
              FOOTCUM(KP)=FOOTCUM(KP)+FOOT
!              WRITE(45,*)KP,JET,RELZIERR,FOOTCUM(KP),FOOT
            ELSE
! CHG:(03/18/2004)increment cumulative FOOT
              FOOTCUM(KP)=FOOTCUM(KP)+FOOT
            END IF

!           computes vertical and horizontal puff mixing
            CALL PUFDSP(METO%HMIX,DT,ZMDL,ZSFC,NLVL,METO%VMIX,ZSG,      &
     &         MASS(:,KP),ZPOS(KP),SIGH(KP),SIGV(KP),HDWP(KP),          &
     &         METO%ZNDX,ISOT)


!           compute solar flux from solar angle and humidity
!           required for resistance gaseous dry deposition
            IF(RDEP) CALL SUNFLX(NLVL,SEA,METO%RHFR,SWF,TR)

!           optional dry, wet, and decay routines
            IF(CDEP)                                                    &
     &         CALL DEPELM(NLVL,METO%ZNDX,DT,ZMDL,ZSFC,ZSG,             &
     &         MASS(:,KP),DEPT,ZPOS(KP),SIGV(KP),PTYP(KP),              &
     &         METO%LAND,METO%AERO,SFCL,METO%USTR,METO%PSI,             &
     &         SWF,HDWP(KP),METO%RAIN,METO%DENS,METO%TEMP,METO%RHFR)

! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
            IF(GRID(PGRD(KP))%LATLON)THEN
               CALL GBL2LL(PGRD(KP),XPOS(KP),YPOS(KP),                  &
     &                              METO%PLAT,METO%PLON)
            ELSE
!              convert end-point grid position to true coordinates
               CALL CXY2LL(GRID(PGRD(KP))%GBASE,XPOS(KP),YPOS(KP),      &
     &            METO%PLAT,METO%PLON, GRID(PGRD(KP))%proj)
            END IF

! JCL:      add BACK as argument to have more complicated conditions
!           that would still alter conc array when JET is decreasing
!           sum values to array
            CALL CONSUM(NUMGRD,METO%PLAT,METO%PLON,DT,JET,ZMDL,ZSFC,    &
     &         CGSIZE,MASS(:,KP),DEPT,ZPOS(KP),SIGH(KP),SIGV(KP),       &
     &         HDWP(KP),PTYP(KP),CSUM,BACK)


!           sum mass for diagnostic analysis
            IF(INITD.EQ.1.OR.INITD.EQ.2)THEN
!              vertical puffs summed into all levels within puff
               SGT=DMAX1(ZPOS(KP)-1.54*SIGV(KP),DBLE(0.0))
               ZZ=ZMDL*(1.0-DMIN1(DBLE(1.0),SGT))
               ZX=(25.0+SQRT(625.0-120.0*(5.0-ZZ)))/60.0
               KT=MIN(MAX(1,NINT(ZX)),NLVL)

               SGB=DMIN1(ZPOS(KP)+1.54*SIGV(KP),DBLE(1.0))
               ZZ=ZMDL*(1.0-DMIN1(DBLE(1.0),SGB))
               ZX=(25.0+SQRT(625.0-120.0*(5.0-ZZ)))/60.0
               KB=MIN(MAX(1,INT(ZX)),NLVL)
            ELSE
!              particles summed into nearest meteo index
               KB=MAX0(1, MIN0(NZM, NINT(METO%ZNDX) ))
               KT=KB
            END IF

            FRAC=KT-KB+1.0
            DO KZ=KB,KT
               PMASS=MASS(1,KP)
               ZMASS(KZ)=ZMASS(KZ)+MASS(1,KP)/FRAC
               MM=MAXDIM
               DO WHILE(MM.GT.1)
                  PMASS=PMASS+MASS(MM,KP)
                  ZMASS(KZ)=ZMASS(KZ)+MASS(MM,KP)/FRAC
                  MM=MM-1
               END DO
               TMASS=TMASS+PMASS/FRAC
            END DO

!           mark zero mass for removal
            IF(PMASS.EQ.0.0)PGRD(KP)=0

! JCL:      skip the steps reserved for particles which have moved off grid
            GO TO 400

! JCL:      Following are steps reserved for particles that have moved off grid
!CCCCCCCCCCCCCCCCCCCCCCC
!           terminated particles skip to here
  200       CONTINUE

! JCL:      also increment particle age even if it has moved off grid
!           increment particle age after advection
            PAGE(KP)=PAGE(KP)+DT

! JCL:(4/27/01)increment counter of number of particles that have moved off grid
            COUNTNPAROUT=COUNTNPAROUT+1

! JCL:      not write out particle results when particle has moved off grid
!           so branch off to end of particle loop
            GO TO 500

!CCCCCCCCCCCCCCCCCCCCCCC
! JCL:
  400       CONTINUE

! JCL:(11/1/02)use Draxler formulation, with 'terrain compression factor'
! JCL:      follwing lines are to write out each particle's position to file
!           convert sigma to AGL
! CHG(09/10/03) correct transformation between sigma and agl
!           ZOUT=(1.0-ZPOS(KP))*(ZMDL-ZSFC)
!           ZOUT=(1.0-ZPOS(KP))*ZMDL*(ZMDL/(ZMDL-ZSFC))
            ZOUT=(1.0-ZPOS(KP))*(ZMDL-ZSFC)

! JCL:(3/1/01)convert 'WWOUT' from [sigma/min]=>[m/s]
!           multiply by (-1.0*DT/ABS(DT)) b/c want to keep sign of direction
            WWOUT=(-1.0*DT/ABS(DT))*WWOUT*(ZMDL-ZSFC)/60.0

! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
            IF(GRID(PGRD(KP))%LATLON)THEN
              CALL GBL2LL(PGRD(KP),XPOS(KP),YPOS(KP),YOUT,XOUT)
            ELSE
! JCL:        call CXY2LL to convert positions to lat/lon before dump
              CALL CXY2LL(GRID(PGRD(KP))%GBASE,XPOS(KP),YPOS(KP),       &
     &                     YOUT,XOUT, GRID(PGRD(KP))%proj)
            END IF

! JCL:(2/28/01)add call to SUNANG b/c need solar angle to output incident solar radiation
!               feed current position instead of starting position, as original call does
! CHG:(11/20/01) use downward short.wave radiation from met data,
!           if available
! CHG:(9/24/02) this is flux interpolated from nearest 2 analysis times
! NEW: interpolation of cloud effect rather then interpolation of radiation
         IF(METO%DSWF.GE.0)THEN
               SWF=METO%DSWF
         ELSE
            CALL SUNANG                                                 &
     &         (SPOT(1)%IBYR,JET,YOUT,XOUT,EA,SEA)

! JCL:(2/28/01)add call to SUNFLX to derive solar flux
!           compute solar flux from solar angle and humidity
!           required for resistance gaseous dry deposition

            CALL SUNFLX(NLVL,SEA,METO%RHFR,SWF,TR)
         END IF

! JCL:      write out results to PARTICLE.DAT only in specified intervals
!             and OUTDT=0 means that results would be output EACH timestep
!           'WRITEOUT' is a logical flag that gets set before particle loop
            IF(WRITEOUT)THEN

! CHG: (20/11/01) write also rain
! CHG: (08/11/02) CHANGED FOR H2O budget study
! JCL:      write to 'particle.dat'  1)mins since start
!           2)particle index  3)particle LAT  4)particle LON
!           5)particle altitude  6) terrain height
!           7)standard deviation of vertical velocity [m/s]
!           8)air temperature at lowest model layer [K]
!           9)amount of time that particle 'sees' the ground [min]
!           10)downward solar short wave radiation [W/m2]
!           11)vertical mean wind [m/s] -->CHANGED to temp [k]
!           12)ht of mixed-layer [m]
!           13)total rain fall rate [m/min]
!           14)convective rain fall rate [m/min]
!           15)sensible heat flux [w/m2] --> CHANGED to loc. dens
!           16)latent heat flux [w/m2] (NGM: kg/m2/s)
!           17)total cloud cover (%)
!           18)particle weighting due to mass violation
!           19)Lagrangian timescale [s] -->CHANGED to RH fraction
! CHG:         Only write all the stuff for particle trajectories
! CHG: (20/11/01) write also rain


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG(11/27/02): get average RH between surface and zmix
! NEED average H2O MIXING RATIO between surface and zmix
! NEED average H2O MIXING RATIO between surface and zmix
! NEED average H2O MIXING RATIO between surface and zmix
! NEED average H2O MIXING RATIO between surface and zmix!
! mixed-layer ht in SIGMA coordinates
!            ZMLSIGMA=1.0-((METO%ZMLNEXT+METO%ZMLPREV)/2.0/ZMDL)*
!     &               ((ZMDL-ZSFC)/ZMDL)
! determine vertical index of mixed-layer ht
!            KK=1
!            RHMIX=0
!            ZBOT=0
!            DO WHILE(DNINT(ZSG(KK)*1E4).GT.DNINT(ZMLSIGMA*1E4))
!               KZM=KK+1
!               KK=KK+1
!               ZTOP=(1.0-ZSG(KK))*ZMDL*(ZMDL/(ZMDL-ZSFC))
!               RHMIX=RHMIX+METO%RHFR(KK)*(ZTOP-ZBOT)
!         WRITE(45,*)'RHp',METO%RHFR(KK),ZTOP-ZBOT,RHMIX
!               ZBOT=ZTOP
!            END DO
!            RHMIX=RHMIX/ZTOP
!         WRITE(45,*)'RHfinal',RHMIX
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TEST ONLY

!-------------------------------------------------------------------------
! JCL(060929): determine vertical indices to properly extract the values of variables at proper height
             KZZ=1
             DO WHILE(ZSG(KZZ+1).GT.ZPOS(KP))
               !CHG i.e. if particle below HGT(2), KB=1 and KT=2; also if particle below HGT(1)
               KZZ=KZZ+1
             END DO
             KTT=KZZ+1
             KBB=KTT-1
!            for RAMS use also different vertical index to extract scalars
             IF(RAMSFLG)THEN
               KZZ=1
               DO WHILE(ZSG(KZZ).GT.ZPOS(KP))
                 KZZ=KZZ+1
               END DO
               KTT=KZZ
               KBB=KTT-1
             END IF

! JCL(060929): vertical interpolation
             IF(KBB.EQ.0.AND.RAMSFLG)KBB=1
             ZF=(ZSG(KBB)-ZPOS(KP))/(ZSG(KBB)-ZSG(KTT))
!            if below 1st level, near the ground
             IF(ZPOS(KP)>ZSG(1))ZF=0.0
!            since don't have pressure profile @ t (only have at t+dt), just use the values extracted after mean advection
             PRESLOCAL=METO%PRES(KBB)+ZF*(METO%PRES(KTT)-METO%PRES(KBB))
             TEMPLOCAL=METO%TEMP(KBB)+ZF*(METO%TEMP(KTT)-METO%TEMP(KBB))
!            NOTE: this is more appropriate for density at final transported particle location than
!            METO%DENSLOCAL (which uses vertical position before PARDSP)
             DENSLOCAL=METO%DENS(KBB)+ZF*(METO%DENS(KTT)-METO%DENS(KBB))
             RHFRLOCAL=METO%RHFR(KBB)+ZF*(METO%RHFR(KTT)-METO%RHFR(KBB))

             SURFTEMP = METO%TEMP(1)
             IF(DREC(KG)%TFLG .and. awrfflg) SURFTEMP = METO%T02M

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! JCL(02/28/2004):  Call output routine written by Marcos
! CHG:(03/17/2004) also add output variable for specific humidity SPHU and FOOT
             CALL OUTPUT(50,IVMAX,VARSIWANT(1:IVMAX),NTURB,             &
     &                   JET+INT(DT)-(CONC(1)%START%MACC),KP,ICNDX(KP), &
     &                   XOUT,YOUT,ZOUT,ZSFC,SIGMAW,SURFTEMP,           &
     &                   SAMPTTCUM(KP),FOOTCUM(KP),SWF,WWOUT,           &
     &                   (METO%ZMLNEXT+METO%ZMLPREV)/2.0,METO%RAIN,     &
     &                   METO%CRAI,ZFXCUM(KP),METO%SHTF,METO%LHTF,METO%TCLD, &
     &                   DMASSWT(KP),SAMPTT2,DENSLOCAL,                 &
     &                   RHFRLOCAL,SPHU,METO%SOLW,METO%LCLD,            &
     &                   METO%ZLOCNEXT,PRESLOCAL,TEMPLOCAL)
! JCL:              reinialize the cumulative SAMPTT
                    SAMPTTCUM(KP)=0.0
                    FOOTCUM(KP)=0.0
                    ZFXCUM(KP)=0.0
!-------------------------------------------------------------------------


!CCCCCCCCCCCCCCCCCCCCCCCCCCC JCL(02/28/2004) not needed for dynamic output CCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!              IF(NTURB.EQ.0)THEN
!C JCL(03/31/03):     CHG's output
!C                    WRITE(50,1111) JET+INT(DT)-
!C     &              (CONC(1)%START%MACC),
!C     &              KP,YOUT,XOUT,ZOUT,ZSFC,
!C     &              SIGMAW,METO%TEMP(1),SAMPTTCUM(KP),
!C     &              SWF,METO%TEMP(KB),(METO%ZMLNEXT+METO%ZMLPREV)/2.0,
!C     &              METO%RAIN,METO%CRAI,METO%DENSLOCAL,WHTF,METO%TCLD,
!C     &              DMASSWT(KP),METO%RHFR(KB),METO%SOLW
!C JCL(03/31/03):     JCL's output
!                   WRITE(50,1111) JET+INT(DT)-
!    &              (CONC(1)%START%MACC),
!    &              KP,
!C    &              ICNDX(KP),
!    &              YOUT,XOUT,ZOUT,ZSFC,
!    &              SIGMAW,METO%TEMP(1),SAMPTTCUM(KP),
!    &              SWF,WWOUT,(METO%ZMLNEXT+METO%ZMLPREV)/2.0,
!    &              METO%RAIN,METO%CRAI,SHTF,WHTF,METO%TCLD,
!    &              DMASSWT(KP),SAMPTT2
!C    &              ,RAUP
!C CHG:(9/12/02)not needed all the time...
!C     &              METO%LCLD,METO%ZLOCNEXT
!1111          FORMAT(I9,I7,F11.4,F11.4,F11.2,F9.2,F10.5,
!    &               F10.2,F10.4,F10.2,F12.6,F12.2,F10.7,
!    &               F14.3,F10.2,F10.2,F6.1,F13.7,F17.3,F13.7)
!C     &               ,F13.7)

!C JCL:(4/3/02)      reinialize the cumulative mass violation
!C                    DMASSCUM(KP)=0.0
!C JCL:              reinialize the cumulative SAMPTT
!                   SAMPTTCUM(KP)=0.0
!              END IF

!C CHG:         Don't write all the stuff for mean trajectories (not 7-12)
!              IF(NTURB.EQ.1)THEN
!                   WRITE(50,1112) JET+INT(DT)-
!    &              (CONC(1)%START%MACC),
!    &              KP,YOUT,XOUT,ZOUT,ZSFC,
!    &              METO%TEMP(1),SWF,WWOUT,
!    &              (METO%ZMLNEXT+METO%ZMLPREV)/2.0,METO%RAIN,
!    &              METO%CRAI,SHTF,WHTF,METO%TCLD,METO%DENSLOCAL,
!C     &              METO%SOLW,
!    &              DMASSWT(KP)
!C     &              DMASS


!C JCL:         write out results in nice formatted way
!1112          FORMAT(I9,I7,F15.8,F15.8,F15.6,F9.2,
!    &               F10.2,F10.2,F14.6,F12.2,F10.7,F11.7,F10.2,
!    &               F10.2,F6.1,F10.6,F13.7,F17.3)

! JCL:(4/3/02)      reinialize the cumulative mass violation
!                    DMASSCUM(KP)=0.0
! JCL:              reinialize the cumulative SAMPTT
!                    SAMPTTCUM(KP)=0.0
!               END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCC JCL(02/28/2004) not needed for dynamic output CCCCCCCCCCCCCCCCCCCCCCCCCCCC

! JCL:      end of IF(WRITEOUT) statement
            END IF

! JCL:      branch from particles which have moved off grid
  500       CONTINUE
!        particle/puff loop
         END DO
!CCCCCCCCC

! JCL:   reset flag and output counter
         IF(WRITEOUT)THEN
            WRITEOUT=.FALSE.
            COUNTOUT=0.0
         END IF

! CHG:(12/05/01) reset temporary convtmp variable
            CONVTMP=0

!        decay of ground-level deposition amounts
         IF(CDEP) CALL DEPRAD(NUMGRD,NUMTYP,DT,CSUM)

!        increment elapsed time
         JET=JET+INT(DT)

!        check each time step for output to disk and zero-out array
! JCL: add BACK as argument to subroutine to avoid some tests that
!       would cause program to terminate when program is run backwards
         CALL CONDSK(NUMGRD,NUMTYP,JET,IFHR,CSUM,BACK)
!        zero-out array if required
         CALL CONZRO(NUMGRD,NUMTYP,DT,JET,IFHR,CSUM,BACK)

!        diagnostics to message file
         WRITE(30,'(1X,I4,I12,I6,E12.2)')                               &
     &      KH,JET,KPM,TMASS

! JCL:(4/27/01) branch out of timestep loop and stop program if fraction of particles
!           leaving model area exceeds specified fraction
!         WRITE(45,*)KS,COUNTNPAROUT,IDINT(OUTFRAC*KPM),
!     &                    MAX(IDINT(OUTFRAC*KPM),1)
         IF(COUNTNPAROUT > MAX(IDINT(OUTFRAC*KPM),1))GO TO 300

!     terminate number of time steps per hour loop
      END DO
!CCCCCCCCCCCCCCCCCCCC

!     special diagnostics dump
      IF(MOD(KH, 6).EQ.0.AND.TMASS.GT.0.0)THEN
         DO KZ=NLVL,1,-1
            PRCNT=100.0*ZMASS(KZ)/TMASS
            IF(PRCNT.GT.0.0)THEN
               HEIGHT=(1.0-ZSG(KZ))*ZMDL
               WRITE(30,'(40X,I3,F8.1,F6.2)')KZ,HEIGHT,PRCNT
            END IF
         END DO
      END IF

! JCL: comment out the following PUFF splitting & management statements, since
!      working mainly with PARTICLES
!=>puff splitting routines called each hour

!     vertical puff splitting
!     CALL PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,
!    :   HDWP,PAGE,PTYP,PGRD,NSORT,NUMSPL)

!     horizontal puff splitting
!     CALL PUFSPH(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,
!    :   HDWP,PAGE,PTYP,PGRD,NSORT,NUMSPL)

!=>puff management routines called each hour

!     eliminate unused puffs (pgrd=0)
!     CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,
!    :   PGRD,NSORT,KHMAX,0.0)

!     sort puffs by position before merge (fractions: h,v,t)
!     CALL PUFSRT(FRHS,FRVS,FRTS,0.0,KPM,ZMDL,XPOS,YPOS,ZPOS,
!    :   MASS,SIGH,SIGV,HDWP,PTYP,PGRD,NSORT)

!     merge puffs together if they are in the same position
!     CALL PUFMRG(FRHS,FRVS,FRTS,0.0,KPM,ZMDL,MASS,
!    :   XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,PGRD,NSORT)

!     eliminate unused puffs (pgrd=0)
!     CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,
!    :   PGRD,NSORT,KHMAX,0.0)

!=>less frequent but less restictive merge of low mass particles

!     IF(KPM.GT.MAXPAR/4.AND.MOD(KH,KRND).EQ.0)THEN
!        find the frme mass level for puff merging  - fmass
!        find the frmr mass level for puff deletion - rmass
!        CALL PUFMIN(KPM,MASS,NSORT,TMASS,FRME,FMASS,FRMR,RMASS)

!        sort puffs by position before merge (fractions: h,v,t)
!        CALL PUFSRT(FRHE,FRVE,FRTE,FMASS,KPM,ZMDL,XPOS,YPOS,ZPOS,
!    :      MASS,SIGH,SIGV,HDWP,PTYP,PGRD,NSORT)

!        merge puffs (mass<fmass) if they are in the same position
!        CALL PUFMRG(FRHE,FRVE,FRTE,FMASS,KPM,ZMDL,MASS,
!    :      XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,PGRD,NSORT)

!        eliminate unused puffs (pgrd=0) and those with mass <rmass
!        CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,
!    :      PGRD,NSORT,KHMAX,RMASS)
!     END IF

!     optional conversion of puffs to particles
      IF(PRESET)CALL PUFPAR(KPM,SIGH,SIGV,SIGX,HDWP,PGRD,NUMPAR)

!     if particles go to zero reset data file pointers
      IF(KPM.EQ.0.AND.KGRID.NE.0)THEN
         KGRID=0
         WRITE(30,*)'WARNING hymodelc: file pointer reset'
      END IF

!     percent completion
      PCNT=100.0*KH/NHRS
      WRITE(*,'(1X,A,F5.1)')'Percent complete: ',PCNT

!     terminate hours to run loop
      END DO

! JCL:(4/27/01) branch here if fraction of particles which have moved off grid exceeds OUTFRAC
  300 CONTINUE

! JCL:close 'JCLMESSAGE'
      CLOSE(45)

! JCL:close 'PARTICLE.DAT'
      CLOSE(50)

! JCL:(11/03/03)
      IF(WINDERRTF.EQ.1)DEALLOCATE(UVERR)

!=>optional simulation dump file for model reinitialization

      IF(NDUMP.EQ.1.OR.NDUMP.EQ.3)THEN

         WRITE(32)KPM,NUMTYP
         DO J=1,KPM
! CHG&JCL (03/10/2004) avoid using PGRD(KG) as index when eqal to 0 (i.e. off grid)
           IF(PGRD(J).NE.0)THEN

! JCL:(07/12/2004) implemented code from HYSPLIT Vers 45
!            missing not permitted in this context
             K=MAX0(KG,KGRID,PGRD(J))
! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
             IF(GRID(K)%LATLON)THEN
              CALL GBL2LL(K,XPOS(J),YPOS(J),TLAT,TLON)
             ELSE
! CHG (06/28/2004) changed from GRID(PGRD(KG)) to GRID(PGRD(J)) (was a bug fixed in Vers 4.5)
              CALL CXY2LL(GRID(PGRD(J))%GBASE,XPOS(J),YPOS(J),TLAT,TLON, GRID(PGRD(J))%proj)
                    !IF(GRID(K)%LATLON)THEN
             END IF

           ELSE
             XPOS(J)=0.0
             YPOS(J)=0.0
           END IF

! JCL:      this is NOT height in AGL! No surface height included in conversion
!           convert height to sigma ASSUMING TERRAIN HEIGHT = 0!!!
! jcl       simply the reverse transformation is done when reading in; as long as consistent, doesn't matter
            ZPOS(J)=(1.0-ZPOS(J))*ZMDL

            WRITE(32)(MASS(I,J),I=1,NUMTYP)
            WRITE(32)TLAT,TLON,ZPOS(J),SIGH(J),SIGV(J),SIGX(J)
            WRITE(32)PAGE(J),HDWP(J),PTYP(J),PGRD(J),NSORT(J)
! CHG(11/05/03) also write DMASSWT(KP)
            WRITE(32)DMASSWT(J)

         END DO
         CLOSE(32)

         WRITE(30,*)'NOTICE: created PARDUMP file for ',KPM,' particles'
      END IF

!=>termination message
      WRITE(*,*)'Complete Hysplit'

! JCL:(05/12/2004) turn dynamic allocation of starting locations back on
      DEALLOCATE(XARG)
      DEALLOCATE(YARG)

!=>required for NCEP operational implementation
!     CALL W3TAGE('HYMODELC')

END PROGRAM HYMODELC
