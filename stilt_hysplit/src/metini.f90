!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METINI           METeorological INItialization
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL INITIALIZATION OPENS THE METEO FILES THE FIRST TIME
!   ON EACH DEFINED GRID SYSTEM; DEFINES DEFAULT GRID STRUCTURE WHICH
!   CANNOT CHANGE WITH TIME FOR A DEFINED GRID NUMBER. MULTIPLE FILE
!   TIME ARE FOR THE SAME GRID NUMBER AND REQUIRE THE SAME STRUCTURE
!   DEFINITIONS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  15 Apr 1999 (RRD) - added surface height definition
!
! USAGE:  CALL METINI(HEADER,NGRD,OLAT,IBYR,IBMO,IBDA,IBHR,BACK,KVEL)
!   INPUT ARGUMENT LIST:
!     HEADER    - char  extended header of index record
!     NGRD      - int   number of meteo grids defined (<=mgrd)
!     OLAT      - real  calculation source latitude
!     IBYR...   - int   calculation starting time
!     BACK      - log   integration direction indicator
!     KVEL  - int vertical motion method indicator
!   OUTPUT ARGUMENT LIST:
!     COMMON GBLGRD
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metini.f90,v 1.13 2008-05-30 16:03:16 tnehrkor Exp $
!
!$$$

SUBROUTINE METINI (HEADER,NGRD,OLAT,IBYR,IBMO,IBDA,IBHR,BACK,KVEL)

   USE module_defgrid
   USE map_utils
   IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'

   CHARACTER HEADER*(*)
   LOGICAL BACK

!      COMMON /GBLGRD/ GRID, DREC, FILE

!  Need local real*4 variables for interface with WPS modules:
   REAL              :: dxm, true1, true2, sync_lat, synclon180, reflon180, sync_xp, sync_yp
   REAL, ALLOCATABLE :: xlat_array(:,:)
   INTEGER           :: iproj, i_xlat, j_xlat
   REAL              :: xi, xj, xlon, xlat


!---------------------------------------------------------------------------------------------------
   DO KG=1,NGRD

!==>file parameters

!     record counter offset usually set to zero (see metold)
!     however old style NH/SH data sets require the offset to
!     position to the SH data if the starting lat indicates SH

      DREC(KG)%OFFSET = 0
      IF (OLAT < 0.0) DREC(KG)%OFFSET = INT(OLAT)

!     unit numbers increment by two for each grid
      FILE(KG,1)%KUNIT = 70+(KG-1)*2
!     initializes time period for index=1
      CALL METSET (HEADER,KG,1)

!     copy structure to 2nd time period, because multiple time
!     periods using the same KG are required to be on the same grid
!     also, initialize array index pointer, indicies are reversed because
!     they are incremented (reversed) on call to METPOS
      FILE(KG,1)%PERIOD = 2
      FILE(KG,2)        = FILE(KG,1)
      FILE(KG,2)%PERIOD = 1

!==>surface parameters

!     set flags for the type of surface data available
      DREC(KG)%UFLG = .FALSE.
      DREC(KG)%TFLG = .FALSE.
!     DREC(KG)%SFLG=.FALSE.
!     DREC(KG)%EFLG = .FALSE.
!     DREC(KG)%STAR=.FALSE.

! JCL:(5/7/01) implement Draxler's 'flux test correction' made in July 28, 2000
!     surface exchange coefficient, heat, momentum flux
      DREC(KG)%EFLX = .FALSE.
      DREC(KG)%HFLX = .FALSE.
      DREC(KG)%UFLX = .FALSE.
!     friction velocity and temperature
      DREC(KG)%USTR = .FALSE.
      DREC(KG)%TSTR = .FALSE.

!     surface terrain height and pressure
      DREC(KG)%SHGT = .FALSE.
      DREC(KG)%PRSS = .FALSE.

!     set default precipitation accumulation time (default = none)
      DREC(KG)%ACYCLE = 0

! CHG:(20/11/01) set default conv. prec. flag
      DREC(KG)%CFLG = .FALSE.
! CHG:(20/11/01) set default total cloud flag
      DREC(KG)%TCLF = .FALSE.
! CHG:(12/04/01) set default low cloud flag
      DREC(KG)%LCLF = .FALSE.
! CHG:(20/11/01) set default downw. rad. flag
      DREC(KG)%RADF = .FALSE.
! CHG:(22/01/03) set default soil moisture flag
      DREC(KG)%SLMF = .FALSE.

!        scan index record settings for certain variable ids
!        if one component there assume the other is as well

      NVAR = DREC(KG)%NUM_VARB(1)
      DO I=1,NVAR
!        surface height (shgt) or pressure (prss) is available
         IF (DREC(KG)%VARB_ID(I,1) == 'SHGT') DREC(KG)%SHGT = .TRUE.
         IF (DREC(KG)%VARB_ID(I,1) == 'PRSS') DREC(KG)%PRSS = .TRUE.

!        two or ten meter temperatures available
         IF (DREC(KG)%VARB_ID(I,1) == 'TMPS') DREC(KG)%TFLG = .TRUE.
         IF (DREC(KG)%VARB_ID(I,1) == 'T02M') DREC(KG)%TFLG = .TRUE.

!        ten meter winds available
         IF (DREC(KG)%VARB_ID(I,1) == 'U10M') DREC(KG)%UFLG = .TRUE.

!        momentum exchange coefficients (assume heat available)
!        negative checksum (eta fix) indicates variable not valid
!        IF(DREC(KG)%VARB_ID(I,1).EQ.'UMOF'.AND.
!     :     DREC(KG)%CHK_SUM(I,1) >= 0) DREC(KG)%SFLG = .TRUE.

! JCL:(5/7/01) implement Draxler's 'flux test correction' made in July 28, 2000
!        momentum flux (checksum positive for valid)
         IF (DREC(KG)%VARB_ID(I,1) == 'UMOF' .AND.                     &
             DREC(KG)%CHK_SUM(I,1) >= 0) DREC(KG)%UFLX = .TRUE.

!        heat flux (checksum positive for valid)
         IF (DREC(KG)%VARB_ID(I,1) == 'SHTF' .AND.                     &
             DREC(KG)%CHK_SUM(I,1) >= 0) DREC(KG)%HFLX = .TRUE.
         IF (DREC(KG)%VARB_ID(I,1) == 'HFLX' .AND.                     &
             DREC(KG)%CHK_SUM(I,1) >= 0) DREC(KG)%HFLX = .TRUE.

!        special scalar momentum from NGM fields
         IF (DREC(KG)%VARB_ID(I,1) == 'EXCO') THEN
! JCL:(5/7/01) implement Draxler's 'flux test correction' made in July 28, 2000
!           DREC(KG)%EFLG = .TRUE.
!           DREC(KG)%SFLG=.TRUE.
            DREC(KG)%EFLX = .TRUE.
            DREC(KG)%UFLX = .TRUE.
         END IF

!           check for normalized momentum and temperature profiles
!           assume if USTR or TSTR available, then both fields are there
! JCL:(5/7/01) implement Draxler's 'flux test correction' made in July 28, 2000
!        IF(DREC(KG)%VARB_ID(I,1) == 'USTR') DREC(KG)%STAR = .TRUE.
         IF (DREC(KG)%VARB_ID(I,1) == 'USTR') DREC(KG)%USTR = .TRUE.
         IF (DREC(KG)%VARB_ID(I,1) == 'TSTR') DREC(KG)%TSTR = .TRUE.

!           cycle time that the precip bucket is emptied (min)
         IF (DREC(KG)%VARB_ID(I,1) == 'TPPA' .OR.                        &
             DREC(KG)%VARB_ID(I,1) == 'TPPS') THEN
!            precip accumulation or summed over entire data set
            DREC(KG)%ACYCLE = 99999
         ELSEIF(DREC(KG)%VARB_ID(I,1) == 'TPPD') THEN
!           accumulation daily over 24 hours
            DREC(KG)%ACYCLE = 1440
         ELSEIF(DREC(KG)%VARB_ID(I,1) == 'TPPT') THEN
!           accumulation over 12 hours
            DREC(KG)%ACYCLE = 720
         ELSEIF(DREC(KG)%VARB_ID(I,1) == 'TPP6') THEN
!           accumulation over 6 hours
            DREC(KG)%ACYCLE = 360
         ELSEIF(DREC(KG)%VARB_ID(I,1) == 'TPP3') THEN
!              accumlation over 3 hours
            DREC(KG)%ACYCLE = 180
         END IF

         ! convective precip available?
         IF (ANY(DREC(KG)%VARB_ID(I,1) == (/'CPPT','CPPD','CPP3','CPP6','CPRC'/))) &
            DREC(KG)%CFLG = .TRUE.

! CHG:(11/20/01) total cloud cover available
         IF (DREC(KG)%VARB_ID(I,1) == 'TCLD') DREC(KG)%TCLF = .TRUE.

! CHG:(12/04/01) low cloud cover available
         IF (DREC(KG)%VARB_ID(I,1) == 'LCLD') DREC(KG)%LCLF = .TRUE.

! CHG:(11/20/01) shortwave radiative flux available
         IF (DREC(KG)%VARB_ID(I,1) == 'DSWF') DREC(KG)%RADF = .TRUE.

! CHG:(11/20/01) soil moisture available
         IF (DREC(KG)%VARB_ID(I,1) == 'SOLW') DREC(KG)%SLMF = .TRUE.

      END DO


!==>upper level parameters (start at level 2)

      NLVL = GRID(KG)%NZ
      DO J=2,NLVL

!           default flags for other meteorological data
         DREC(KG)%QFLG      = .FALSE.
         DREC(KG)%RFLG(J-1) = .FALSE.
         DREC(KG)%WFLG(J-1) = .FALSE.

!           shift index down to corresspond with upper level variables
         NVAR = DREC(KG)%NUM_VARB(J)
         DO I=1,NVAR
!           specific or relative humidity
            IF (DREC(KG)%VARB_ID(I,J) == 'SPHU') DREC(KG)%QFLG = .TRUE.
!           set vertical motion flag if field found
            IF (DREC(KG)%VARB_ID(I,J) == 'WWND' .or. DREC(KG)%VARB_ID(I,J) == 'DZDT')  &
                DREC(KG)%WFLG(J-1) = .TRUE.
!           specific humidity available rather than relative humidity
            IF (DREC(KG)%VARB_ID(I,J) == 'SPHU' .OR.                   &
                DREC(KG)%VARB_ID(I,J) == 'RELH')                      &
                DREC(KG)%RFLG(J-1) = .TRUE.
         END DO
      END DO

!        check selected vertical motion method versus available field
!        no selection and no data default to isobaric (kvel = 1)
      IF ((.NOT.DREC(KG)%WFLG(2)) .AND. KVEL == 0) KVEL = 1

!        surface pressure or terrain required for calculations
      IF (.NOT.(DREC(KG)%SHGT .OR. DREC(KG)%PRSS)) THEN
         WRITE (*,*) 'ERROR Metini: input meteorology requires either'
         WRITE (*,*) '   Surface pressure: ',DREC(KG)%PRSS
         WRITE (*,*) '   Terrain height  : ',DREC(KG)%SHGT
         STOP
      END IF

!==>initialize grid conversion subroutines

      WRITE (*,*) "MODEL_ID: ",GRID(KG)%MODEL_ID
      WRITE (*,*) "TANG_LAT: ",GRID(KG)%TANG_LAT
      WRITE (*,*) "REF_LAT:  ",GRID(KG)%REF_LAT
      WRITE (*,*) "REF_LON:  ",GRID(KG)%REF_LON

      if (GRID(KG)%pole_lat .le. vmissle .and. GRID(KG)%pole_lon .le. vmissle &
             .and. GRID(KG)%ref_lat .le. vmissle) then
         GRID(KG)%LATLON = .FALSE.
         GRID(KG)%GLOBAL = .FALSE.
         ! using WPS geolocation routines
         dxm = GRID(KG)%size*1000.
         true1 = GRID(KG)%tang_lat
         true2 = GRID(KG)%dummy
         iproj = nint(GRID(KG)%orient)
         synclon180 = GRID(KG)%sync_lon
         if (synclon180 > 180) synclon180 = synclon180-360
         reflon180 = GRID(KG)%ref_lon
         if (reflon180 > 180) reflon180 = reflon180-360
         sync_lat = GRID(KG)%sync_lat
         sync_xp = GRID(KG)%sync_xp
         sync_yp = GRID(KG)%sync_yp
         call map_set(proj_code = iproj,proj = GRID(KG)%proj,&
                lat1 = sync_lat,lon1 = synclon180,knowni=sync_xp,knownj=sync_yp, &
                dx = dxm, stdlon = reflon180, truelat1=true1, truelat2=true2)
         GRID(KG)%gbase(:) = vmiss
         write(*,*) "metini: Using WPS geolocation routines with proj_code: ",iproj
         ! Initialize GRID(KG)%mapfactor
         allocate(xlat_array(grid(kg)%nx,grid(kg)%ny))
         allocate(grid(kg)%mapfactor(grid(kg)%nx,grid(kg)%ny))
         do j_xlat=1,grid(kg)%ny
            xj = j_xlat
            do i_xlat=1,grid(kg)%nx
               xi = i_xlat
               call ij_to_latlon(grid(kg)%proj, xi, xj, xlat, xlon)
               xlat_array(i_xlat,j_xlat) = xlat
            end do
         end do
         call get_map_factor(grid(kg)%proj,xlat_array,grid(kg)%mapfactor, &
                1, 1, grid(kg)%nx, grid(kg)%ny)
         deallocate(xlat_array)

      else

! JCL:(07/06/2004) test for global latlon grid system
         CALL GBLSET (KG)

         WRITE (*,*) 'after GBLSET:',GRID(KG)%LATLON,GRID(KG)%GLOBAL

!!!!!!!
!      GRID(KG)%REF_LAT=-3.5
!      GRID(KG)%REF_LON=-59.0
!        temporary patch for RAMS - assume that the RAMS oblique
!        stereographic grid is really Lambert Conformal until additional
!        oblique library routines are complete
!         IF(GRID(KG)%MODEL_ID.EQ.'RAMS')
!     :      GRID(KG)%TANG_LAT = GRID(KG)%POLE_LAT
!(03/11/03)JCL: use DMAPf routines to deal with oblique stereographic grid
         IF (GRID(KG)%MODEL_ID == 'RAMS') THEN
            WRITE (*,*) "yes RAMS"
!           initialize grid conversion variable (into gbase)
!            CALL SOBSTR(GRID(KG)%GBASE,
!     &                  GRID(KG)%TANG_LAT,GRID(KG)%REF_LON)
            CALL SOBSTR (GRID(KG)%GBASE,GRID(KG)%REF_LAT,GRID(KG)%REF_LON)
         ELSE
            WRITE (*,*) "not RAMS"
! JCL:(07/06/2004) not on global latlon grid system (flag set in GBLSET)
            IF (.NOT.GRID(KG)%LATLON) THEN

!           initialize grid conversion variable (into gbase)
! CHG(03/11/03) changed to version 2.0 of cmapf
               CALL STLMBR (GRID(KG)%GBASE,GRID(KG)%TANG_LAT, GRID(KG)%REF_LON)

!              CALL STLMBR(GRID(KG)%GBASE,
!     :                    GRID(KG)%REF_LAT, GRID(KG)%REF_LON)
!           CALL STCMAP(GRID(KG)%GBASE,
!     :         GRID(KG)%TANG_LAT, GRID(KG)%REF_LON)
            END IF
         END IF
!!!!!!!
!      GRID(KG)%SYNC_LON=-73.36
!      GRID(KG)%SYNC_LAT=-11.18

! JCL:(07/06/2004) not on global latlon grid system (flag set in GBLSET)
         IF (.NOT.GRID(KG)%LATLON) THEN
!           use single point grid definition
            CALL STCM1P (GRID(KG))
         END IF
      ENDIF
      WRITE (*,*) "SYNC_XP,SYNC_YP:",GRID(KG)%SYNC_XP,GRID(KG)%SYNC_YP
      WRITE (*,*) "SYNC_LON,SYNC_LAT:",GRID(KG)%SYNC_LON,GRID(KG)%SYNC_LAT
      WRITE (*,*) "SIZE,NX,NY:",GRID(KG)%SIZE,GRID(KG)%NX,GRID(KG)%NY

   END DO

!==>set proper starting time if mo=0 set the initial model start time
!   assume (grid 1 always defined) and day and hour relative

   IF (IBMO == 0) THEN
      IF (BACK) THEN
         IBYR = FILE(1,1)%LAST%YR
         IBMO = FILE(1,1)%LAST%MO
         IBDA = FILE(1,1)%LAST%DA-IBDA
         IBHR = FILE(1,1)%LAST%HR-IBHR
      ELSE
         IBYR = FILE(1,1)%FIRST%YR
         IBMO = FILE(1,1)%FIRST%MO
         IBDA = FILE(1,1)%FIRST%DA+IBDA
         IBHR = FILE(1,1)%FIRST%HR+IBHR
      END IF

!     adjust relative date for potential month crossing error
      IBMN = 0
      CALL TM2MIN (IBYR,IBMO,IBDA,IBHR,IBMN,MACC)
      CALL TM2DAY (MACC,IBYR,IBMO,IBDA,IBHR,IBMN)
   END IF


END SUBROUTINE METINI
