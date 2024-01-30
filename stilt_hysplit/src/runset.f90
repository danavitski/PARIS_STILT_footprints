!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  RUNSET           RUN SETup for basic model parameters
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   RUN SETUP SETS THE MOST BASIC MODEL SIMULATION PARAMETERS, BEFORE
!   ANY OTHER INFORMATION IS KNOWN. IF THERE EXISTS A FILE NAMED CONTROL
!   THEN ALL INPUT PARAMETERS WILL BE READ FROM THAT FILE.  OTHERWISE
!   INPUT IS EXPECTED ON STANDARD INPUT (UNIT 5).  IN THAT CASE AN OUTPUT
!   FILE CALLED STARTUP WILL BE CREATED THAT WILL CONTAIN ALL DATA ENTRIES
!   IT MAY BE USED IN SUBSEQUENT SIMULATIONS AS THE CONTROL FILE.
!   ENTRIES FROM STANDARD INPUT HAVE DEFAULT VALUES WHICH MAY BE SELECTED
!   BY JUST ENTERING "/".  THIS IS NOT VALID FOR CHARACTER STRINGS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 21 Jul 1998 (RRD)
!                 09 Mar 1999 (RRD) - additonal subgrid parameters added
!                 08 Apr 1999 (RRD) - location specific emissions
!                 22 Apr 1999 (RRD) - renamed surface layer depth variable
!                                   - optional source/area at release location
!                 26 Aug 1999 (RRD) - common block to define vertical grid
!
! USAGE:  CALL RUNSET(NLOC,NHRS,NGRD,NLVL,KSFC,ZSG,ZMDL,SFCL,
!              BACK,IUNIT,KVEL,KLEN,FNAME)
!   INPUT ARGUMENT LIST:
!     NONE
!   OUTPUT ARGUMENT LIST:
!     NLOC      - int   number of starting locations
!     NHRS      - int   simulation duration in hours
!     NGRD      - int   number of data grids for this run
!     NLVL      - int   number of internal levels
!     KSFC      - int   index number for top of sfc layer
!     ZSG       - real  internal terrain following sigma values
!     ZMDL      - real  maximum height for scaling coordinate system
!     SFCL      - real  height of the surface layer top (m)
!     BACK      - log   integration direction flag
!     IUNIT     - int   unit for file of input parameters
!     KVEL      - int   vertical velocity remapping
!     KLEN      - int   length of unique file name string
!     FNAME     - char  string giving unique file name to input data
!     AA,BB,CC  - real  coefficients of the quadratic equation for the vertical model levels
!     COMMON GBLSPT
!     COMMON GBLGRD
!   INPUT FILES:
!     UNIT 5 or 30 depending if CONTROL file is found
!   OUTPUT FILES:
!     UNIT 31 when input defined on unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: runset.f90,v 1.8 2007-05-03 13:09:13 skoerner Exp $
!
!$$$
      SUBROUTINE RUNSET (NLOC, NHRS, NGRD, NLVL, KSFC, ZSG, ZMDL, SFCL,         &
                         BACK, IUNIT, KVEL, KLEN, FNAME, AA, BB, CC)

      USE module_defspot
      USE module_defmeto
      USE module_defgrid

      IMPLICIT REAL*8 (A-H,O-Z)


!     internal model grid heights and corresponding sigma levels
      REAL(KIND(1d0)), INTENT(INOUT) :: AA, BB, CC

      REAL(KIND(1d0)) :: HGT(NZM), ZSG(NZM), ZDATA=0d0
      LOGICAL BACK, CNTL
      CHARACTER FNAME*(*), LABEL*80


!---------------------------------------------------------------------------------------------------
!=>generic model defaults

      NLOC=1
      SPOT(1)%IBYR=0
      SPOT(1)%IBMO=0
      SPOT(1)%IBDA=0
      SPOT(1)%IBHR=0
      NHRS=48
      KVEL=0
      ZMDL=10000.0
      NGRD=1
      FILE(1,1)%DIR='/main/sub/data/'
      FILE(1,1)%METEO='file_name'

!=>create special file name for input

      IF(KLEN.LE.0)THEN
         LABEL='CONTROL'
      ELSE
         LABEL='CONTROL.'//FNAME(1:KLEN)
      END IF

!=>check for input data control file in local directory

      INQUIRE(FILE=LABEL,EXIST=CNTL)
      IF(CNTL)THEN
         IUNIT=30
         OPEN(IUNIT,FILE=LABEL)
      ELSE
         IUNIT=5
         OPEN(31,FILE='STARTUP')
      END IF

!=>starting time and location information

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Enter starting time (year, month, day, hour)'
         WRITE(*,*)SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR
      END IF
      READ(IUNIT,*)                                                     &
     &   SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR
      IF(IUNIT.EQ.5)                                                    &
     &   WRITE(31,*)SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Enter number of starting locations'
         WRITE(*,*)NLOC
      END IF
      READ(IUNIT,*)NLOC
      IF(IUNIT.EQ.5)WRITE(31,*)NLOC
      IF(NLOC.GT.MLOC)THEN
         WRITE(*,*)'# of starting locations exceeds dimension:',MLOC
         STOP
      END IF

      DO N=1,NLOC
!        default starting locations
         IF(N.EQ.1)THEN
            SPOT(N)%OLAT=40.0
            SPOT(N)%OLON=-90.0
            SPOT(N)%OLVL=50.0
            SPOT(N)%AREA=0.0
            SPOT(N)%QTRM=0.0
         ELSE
            SPOT(N)%OLAT=SPOT(N-1)%OLAT
            SPOT(N)%OLON=SPOT(N-1)%OLON
            SPOT(N)%OLVL=SPOT(N-1)%OLVL
            SPOT(N)%AREA=SPOT(N-1)%AREA
            SPOT(N)%QTRM=SPOT(N-1)%QTRM
         END IF

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Enter starting location (lat, lon, m-agl)'
            WRITE(*,*)SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
         END IF

!=>read standard 3 parameters or variable parameters

         READ(IUNIT,*)SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
!        READ(IUNIT,'(A)')LABEL

!        extract each field from label
!        KK=1
!        KV=0
!        DO WHILE (KK.LT.80.AND.KV.LT.5)
!           KL=INDEX(LABEL(KK:),' ')
!           IF(KL.GT.0)THEN
!              KV=KV+1
!              IF(KV.EQ.1)READ(LABEL(KK:),*)SPOT(N)%OLAT
!              IF(KV.EQ.2)READ(LABEL(KK:),*)SPOT(N)%OLON
!              IF(KV.EQ.3)READ(LABEL(KK:),*)SPOT(N)%OLVL
!              IF(KV.EQ.4)READ(LABEL(KK:),*)SPOT(N)%QTRM
!              IF(KV.EQ.5)READ(LABEL(KK:),*)SPOT(N)%AREA
!              KK=KK+KL-1
!           END IF
!           KK=KK+1
!        END DO

         IF(IUNIT.EQ.5)                                                 &
     &      WRITE(31,*)SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL

!        for now only permit one starting time for all locations
         SPOT(N)%IBYR=SPOT(1)%IBYR
         SPOT(N)%IBMO=SPOT(1)%IBMO
         SPOT(N)%IBDA=SPOT(1)%IBDA
         SPOT(N)%IBHR=SPOT(1)%IBHR
      END DO

!=>simulation run time, vertical coordinate

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Enter total run time (hours)'
         WRITE(*,*)NHRS
      END IF
      READ(IUNIT,*)NHRS
      IF(IUNIT.EQ.5)WRITE(31,*)NHRS

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Vertical (0:data 1:isob 2:isen 3:dens 4:sigma)'
         WRITE(*,*)KVEL
      END IF
      READ(IUNIT,*)KVEL
      IF(IUNIT.EQ.5)WRITE(31,*)KVEL

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Top of model domain (internal coordinates m-msl)'
         WRITE(*,*)ZDATA
      END IF
      READ(IUNIT,*)ZDATA
      IF(IUNIT.EQ.5)WRITE(31,*)ZDATA

!=>meteorological grid information

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Number of input data grids'
         WRITE(*,*)NGRD
      END IF
      READ(IUNIT,*)NGRD
      IF(IUNIT.EQ.5)WRITE(31,*)NGRD

!     test limits
      IF(NGRD.GT.MGRD)THEN
         WRITE(*,*)'ERROR runset: Number of grids exceed limit'
         STOP
      END IF

      DO KG=1,NGRD
         IF(KG.GT.1)THEN
            FILE(KG,1)%DIR=FILE(KG-1,1)%DIR
            FILE(KG,1)%METEO=FILE(KG-1,1)%METEO
         END IF

         IF(IUNIT.EQ.5)THEN
! JCL:(6/1/00)'vf90' compiler doesn't like the '\.' characters
            WRITE(*,'(A,I2,A)')' Enter grid #',KG,' directory'
!           WRITE(*,'(A,I2,A)')' Enter grid #',KG,' directory (\...\)'
            WRITE(*,'(1X,A)')FILE(KG,1)%DIR
         END IF
         READ(IUNIT,'(A)')FILE(KG,1)%DIR
         IF(IUNIT.EQ.5)WRITE(31,'(A)')FILE(KG,1)%DIR

         IF(IUNIT.EQ.5)THEN
            WRITE(*,'(A,I2,A)')' Enter grid #',KG,' file name (?????)'
            WRITE(*,'(1X,A)')FILE(KG,1)%METEO
         END IF
         READ(IUNIT,'(A)')FILE(KG,1)%METEO
         IF(IUNIT.EQ.5)WRITE(31,'(A)')FILE(KG,1)%METEO
      END DO

      PRINT *, 'NHRS in RUNSET:  ', NHRS
!     set integration direction
      BACK=.FALSE.
      IF(NHRS.LT.0)THEN
         BACK=.TRUE.
         NHRS=-NHRS
      END IF

      PRINT *, 'BACK in RUNSET:  ', BACK

!     initialize meteorological subgrid (1,1) location to missing
!     will be updated to correct location after first advection step
      METO%LX1=-1
      METO%LY1=-1

!     subgrid center position
      METO%LXC=-1
      METO%LYC=-1

!     subgrid range
      METO%LXR=NXM
      METO%LYR=NYM


!---------------------------------------------------------------------------------------------------
!  set vertical grid coefficients
!  mechanism: if set in SETUP.CFG - no further change,
!  otherwise set depending on model ID of the meteo input file
!  ARL files in old format get the defaults
      IF (AA < -9999d0) THEN

         IF (BB >= -9999d0 .OR. CC >= -9999d0) &
            STOP 'In SETUP.CFG, either all of AA, BB, CC need to be set, or none. Stop.'

         SELECT CASE (MODEL_ID(TRIM(FILE(1,1)%DIR)//TRIM(FILE(1,1)%METEO)))
         CASE ('ALAG')                                      ! Meteo France ALADIN
            AA =  10d0
            BB =   5d0
            CC =   0d0
            PRINT *, 'Non-default coefficients for the vertical coordinate system'
         CASE DEFAULT
            AA =  30d0
            BB = -25d0
            CC =   5d0
         END SELECT

      END IF

!=>internal model grid terrain following levels m-agl
!  are defined according to the following quadratic equation
!  z =  aa k^2 + bb k + cc  where default aa=30 bb=-25 cc=5

      DO K=1,NZM
         HGT(K)=AA*K*K+BB*K+CC
      END DO

!=>height of the top of the surface layer, used in a variety of
!  scaling applications in several different subroutines

      KSFC=2
      SFCL=HGT(KSFC)

!     check upper limit
      IF(ZDATA.GT.HGT(NZM))THEN
         WRITE(*,*)'WARNING runset: model top exceeds compiled limit'
         WRITE(*,*)'   recompile after increasing NZM in module_defgrid'
         WRITE(*,*)'   or lower model top to below - ',HGT(NZM)
         ZDATA=HGT(NZM)
      END IF

!=>ZMDL also is the terrain scaling parameter, the height at which the
!  sigma surfaces go flat relative to terrain.

      ZMDL=MAX(ZDATA,25000d0)

!=>restrict internal grid to max input array levels as set by NZM
!  and ZDATA because variables share the same array space. ZDATA
!  parameter used to set the limit to how much vertical data
!  will be processed through the NLVL parameter.

      DO K=1,NZM
         IF(HGT(K).LT.ZDATA)THEN
            NLVL=K
!           compute levels as sigma coordinate for terrain=0
!jcl!!!!    correction applied in each subroutine as needed
            ZSG(K)=1.0-HGT(K)/ZMDL
         END IF
      END DO
! CHG(09/16/03) not for RAMS, do it all later in hymodelc
!      DO N=1,NLOC
!        starting height limit
!         IF(SPOT(N)%OLVL.GT.HGT(NLVL))THEN
!            WRITE(*,*)'ERROR runset: Start height above mdl domain'
!            WRITE(*,*)'   Source - ',N,'    Height - ',SPOT(N)%OLVL
!            WRITE(*,*)'   Model domain - ',INT(HGT(NLVL))
!            WRITE(*,*)'   Model top ht - ',ZMDL
!            STOP
!         END IF
!      END DO

      

      CONTAINS


! return the model-ID of ARL file
! doesn't work with old ARL format

      FUNCTION MODEL_ID(FNAME)

      IMPLICIT NONE

      CHARACTER(4)             :: MODEL_ID
      CHARACTER(*), INTENT(IN) :: FNAME

      INTEGER       :: lun
      LOGICAL       :: ftest
      CHARACTER(54) :: str

      INQUIRE(FILE=FNAME,EXIST=ftest)
      IF(.NOT.ftest) THEN
         WRITE(*,*)'Unable to find file: ', FNAME
         WRITE(*,*)'Check input CONTROL file for correct values'
         STOP
      END IF

      !  get free lun
      lun = 10
      DO
         INQUIRE(UNIT=lun, OPENED=ftest)
         IF (.NOT. ftest) EXIT
         lun = lun + 1
      END DO
      OPEN (lun, FILE=FNAME, STATUS='OLD',                              &
            RECL=54,ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
      READ (lun,REC=1) str
      CLOSE (lun)
      READ (str(51:54),'(A4)') MODEL_ID

      END FUNCTION MODEL_ID


      END SUBROUTINE RUNSET
