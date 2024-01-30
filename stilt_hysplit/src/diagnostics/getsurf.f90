      PROGRAM GETSURF
!======================================================================
!     Prints the meteorological profile at a selected grid point from
!     any ARL formatted file.  The grid point nearest to the location
!     is printed.  No interpolation is performed.  Primary purpose of
!     program is for diagnostic evaluation of meteo data files.
!----------------------------------------------------------------------
!     $Id: getsurf.f90,v 1.9 2007-05-03 13:09:13 skoerner Exp $
!=======================================================================

      USE module_defsize
      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL BACK, SXTNBIT, RAMSFLG

!     packed data arrays
      CHARACTER                 :: HEADER*10000, GETLL*1
      CHARACTER(1), ALLOCATABLE :: CDATA(:)

!     unpacked data profile at specific point
      REAL*8 VDATA(MVAR,MLVL)
!     true wind component profile
      REAL*8 UTW(MLVL), VTW(MLVL)


!-----------------------------------------------------------------------
!==>define input meteorological data file

! CHG(3/21/02): don't write too much to std output
      WRITE(*,*)'Enter meteorological directory ...'
      READ(5,'(A)')FILE(1,1)%DIR
      WRITE(*,*)'Enter meteorological file name ...'
      READ(5,'(A)')FILE(1,1)%METEO
      WRITE(*,*)'want lat-lon? ("y" or "n")'
      READ(5,'(A)')GETLL

!     read and initialize data file
      KVEL = 0
      IBMO = 0
      IBDA = 0
      IBHR = 0
      CALL METINI(HEADER,1,0d0,IBYR,IBMO,IBDA,IBHR,.FALSE.,KVEL)
      RAMSFLG = GRID(1)%MODEL_ID == 'RAMS'
      SXTNBIT = RAMSFLG .OR. GRID(1)%MODEL_ID(1:3) == 'ECX'
!     open output file
!      OPEN(30,FILE='profile.txt')
!      WRITE(30,'(2A)')' Meteorological Profile: ',FILE(1,1)%METEO
! CHG(28/3/02): file for whole grid output
      OPEN(35,FILE='info.txt')
      OPEN(40,FILE='p.txt')
!      OPEN(41,FILE='u.txt')
!      OPEN(42,FILE='v.txt')
!      OPEN(45,FILE='w.txt')
      IF (GETLL == "y".OR.GETLL == "Y")OPEN(46,FILE='lat.txt')
      IF (GETLL == "y".OR.GETLL == "Y")OPEN(47,FILE='lon.txt')
      OPEN(48,FILE='dswr.txt')
      OPEN(49,FILE='sfct.txt')
      OPEN(50,FILE='solw.txt')
      OPEN(51,FILE='CPRC.txt')
      OPEN(52,FILE='PPRC.txt')

!==>set default location and starting time of profile

      NXGP=GRID(1)%NX
      NYGP=GRID(1)%NY
      WRITE(35,'(A)')'  NX  NY  NZ'
      WRITE(35,'(3I4)')NXGP,NYGP,GRID(1)%NZ
      IF (SXTNBIT) THEN
         NXYP = NXGP*NYGP*2
      ELSE
         NXYP = NXGP*NYGP
      ENDIF
      ALLOCATE (CDATA(NXYP))

!     set grid center as initial position
      XCP=(NXGP-1.0)/2.0+1.0
      YCP=(NYGP-1.0)/2.0+1.0

      IF (GRID(1)%LATLON)THEN
        CALL GBL2LL(1,XCP,YCP,CLAT,CLON)
      ELSE
        CALL CXY2LL(GRID(1)%GBASE,XCP,YCP,CLAT,CLON)
      END IF

! CHG(3/21/02): don't write too much to std output
      WRITE(*,'(A,5I3)')' File start time :',FILE(1,1)%FIRST%YR,        &
     &   FILE(1,1)%FIRST%MO, FILE(1,1)%FIRST%DA,                        &
     &   FILE(1,1)%FIRST%HR, FILE(1,1)%FIRST%MN
      WRITE(*,'(A,5I3)')' File ending time:',FILE(1,1)%LAST%YR,         &
     &   FILE(1,1)%LAST%MO, FILE(1,1)%LAST%DA,                          &
     &   FILE(1,1)%LAST%HR, FILE(1,1)%LAST%MN

!     set file beginning as initial time
      IYR=FILE(1,1)%FIRST%YR
      IMO=FILE(1,1)%FIRST%MO
      IDA=FILE(1,1)%FIRST%DA
      IHR=FILE(1,1)%FIRST%HR
      IMN=FILE(1,1)%FIRST%MN

!==>select time and location of profile

  100 WRITE(*,*)' '
! CHG(3/21/02): don't write too much to std output
         WRITE(*,'(A,5I3)')' Enter Day & Hour:',IDA,IHR
         READ(*,*,END=200)IDA,IHR
         IF (IDA == 0)GOTO 200

! CHG(3/21/02): don't write too much to std output
!         WRITE(*,'(A,2F10.4)')' Enter Lat & Lon: ',CLAT,CLON
!         READ(*,*,END=200)CLAT,CLON
!         CALL CLL2XY(GRID(1)%GBASE,CLAT,CLON,XCP,YCP)

! CHG(3/28/02): try to read forecasts
!      WRITE(*,*)'XCP:',XCP,' YCP:',YCP
!      WRITE(*,*)'REC_LEN:',FILE(1)%REC_LEN
!      WRITE(*,*)'PERIOD:',FILE(1)%PERIOD
!      WRITE(*,*)'REC_PER:',DREC(1)%REC_PER
!      WRITE(*,*)'DELTA:',DREC(1)%DELTA
!        compute record number of selected time period's index record
!        and last data record for that time period
         CALL TM2MIN(IYR,IMO,IDA,IHR,IMN,JET)
!         JREC1=1+DREC(1)%REC_PER*(JET-FILE(1,1)%FIRST%MACC)
!     :         /DREC(1)%DELTA
! CHG & JCL: 3/21/2002  extra parentheses to divide first
! !!!!!!     always returns profile from nearest earlier analysis
         JREC1=1+DREC(1)%REC_PER*((JET-FILE(1,1)%FIRST%MACC)            &
     &         /DREC(1)%DELTA)

         JREC2=JREC1+DREC(1)%REC_PER-1

! CHG(3/28/02): try to read forecasts
!      WRITE(*,*)'JREC1:',JREC1,' JREC2',JREC2,
!     &           ' REC_PER:',DREC(1)%REC_PER

         IF (JREC1.LT.1.OR.JREC2.GT.FILE(1,1)%ENDREC)THEN
            WRITE(*,*)'Selected time outside of file boundary'
            GOTO 100
         END IF

         IF (NINT(XCP).LT.1.OR.NINT(XCP).GT.NXGP.OR.                     &
     &      NINT(YCP).LT.1.OR.NINT(YCP).GT.NYGP)THEN
            WRITE(*,*)'Selected location outside of file boundary'
            GOTO 100
         END IF

!        dump location information to output file
!         WRITE(30,'(A)')' '
!         WRITE(30,'(A,5I3)')' Profile Time: ',IYR,IMO,IDA,IHR,IMN
!         WRITE(30,'(A,2F10.3)')' Profile Location: ',CLAT,CLON

!==>loop through all the records for that time period

!        load data array for that time period
         CALL DLOAD(SXTNBIT,RAMSFLG,GETLL,CDATA,VDATA,NXGP,NYGP,NXYP,     &
     &              JREC1,JREC2,XCP,YCP,SFCP,SFCT,UTW,VTW)
         GETLL = 'n'                                        ! write lon/lat only once

!        print sounding after all data have been processed
         CALL SOUND(VDATA,SFCP,SFCT,UTW,VTW)

         PRINT *
         PRINT *, 'data dumped to files, type 0,0 to exit'
      GOTO 100

  200 STOP
      END

!--------------------------------------------------------------------


      SUBROUTINE DLOAD(SXTNBIT,RAMSFLG,GETLL,CDATA,VDATA,NXGP,NYGP,NXYP,  &
     &                 JREC1,JREC2,XCP,YCP,SFCP,SFCT,UTW,VTW)

      USE module_defsize
      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

! JCL:'VDATA' stores the processed data

      CHARACTER LABEL*50, CDATA(NXYP)*1, VARB*4, GETLL*1
      REAL*8 RDATA(NXGP,NYGP), VDATA(MVAR,MLVL), UTW(MLVL), VTW(MLVL)
      LOGICAL, INTENT(IN) :: SXTNBIT, RAMSFLG


! CHG: get lat and lon array
      REAL*8 RLAT(NXGP,NYGP), RLON(NXGP,NYGP)


!---------------------------------------------------------------------------------------------------
!     initialize surface pressure and terrain
      SFCP=1013.0
      SFCT=0.0

!     initialize true wind profile
      LP=0
      DO L=1,MLVL
         UTW(L)=0.0
         VTW(L)=0.0
      END DO

!     initialize the data profile
      DO L=1,MLVL
      DO K=1,MVAR
         VDATA(K,L)=0.0
      END DO
      END DO

      IF (GETLL == "y".OR.GETLL == "Y") THEN
         WRITE(*,*)"LAT and LON"
         DO I=1,NXGP
           XX=0.0+I
           DO J=1,NYGP
             YY=0.0+J
             IF (GRID(1)%LATLON)THEN
               CALL GBL2LL(1,XX,YY,RLAT(I,J),RLON(I,J))
             ELSE
               CALL CXY2LL(GRID(1)%GBASE,XX,YY,RLAT(I,J),RLON(I,J))
             END IF
           END DO
         END DO
         WRITE(46,*) RLAT(1:NXGP,1:NYGP)
         WRITE(47,*) RLON(1:NXGP,1:NYGP)
      END IF

!==>loop through each data record in the selected time period

!      DO IREC=JREC1,JREC2
! CHG&JCL: 3/21/2002 record number wrong;rotates wind at 50 mbar now
      DO IREC=JREC1,JREC2+2

!         READ(FILE(1,1)%KUNIT,REC=IREC)LABEL,(CDATA(KK),KK=1,NXYP)


         READ(FILE(1,1)%KUNIT,REC=IREC)LABEL,(CDATA(KK),KK=1,NXYP)
         READ(LABEL,'(7I2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,IL,IG,          &
     &                                   VARB,NEXP,PREC,VAR1

! CHG(3/28/02): try to read forecasts
!      WRITE(*,*)LABEL
!      WRITE(*,*)NXGP,NYGP,NXYP,1,1,NXGP,NYGP,
!     :                  PREC,NEXP,VAR1
!      WRITE(*,*)'NXGP:',NXGP
!      WRITE(*,*)' IL:',IL,' IY:',IY,' IM:',IM,' ID:',ID,' IH:',IH,
!     &          ' IC:',IC,' IL:',IL,' IG:',IG,
!     &          ' VARB:',VARB,' NEXP:',NEXP,' PREC:',PREC,' VAR1:',VAR1


!        only process data records with valid observations (non-missing)
         IF (VARB.NE.'INDX' .AND. VARB.NE.'NULL' .AND. IC.NE.-1)THEN

!           unpack the data
            CALL PAKINP(SXTNBIT,RDATA,CDATA,NXGP,NYGP,NXYP,1,1,         &
     &                  NXGP,NYGP,PREC,NEXP,VAR1,-1)

! CHG(3/28/02): try to read forecasts
!      WRITE(*,*)'LABEL after PAKINP:',LABEL

!           convert level number to array index because input data
!           level index starts at 0 for the surface
            LL=IL+1

! JCL:      Here's where rotation of winds happens!!
!           rotate winds when vertical index changes
            IF (LL.NE.LP  .AND.  .NOT. GRID(1)%LATLON) THEN
               IF (LP.NE.0)THEN
                  CALL CG2CXY(GRID(1)%GBASE,                            &
     &                 XCP,YCP,UTW(LP),VTW(LP),UT,VT)
                  UTW(LP)=UT
                  VTW(LP)=VT
               END IF
               LP=LL
            END IF

!           find the variable array element number - match the input
!           variable with its position as indicated in index record
            NVAR=DREC(1)%NUM_VARB(LL)
            KVAR=0
            DO KK=1,NVAR
               IF (VARB == DREC(1)%VARB_ID(KK,LL))KVAR=KK
!               IF(LL == 1 .AND. IREC == JREC1)WRITE(*,*)VARB
!               IF(LL == 1 .AND. IREC == JREC1)WRITE(*,*)
!     :           DREC(1)%VARB_ID(KK,LL)
            END DO

!           load the meteo data into the array according to index record
!           data structure at the selected grid location
            VDATA(KVAR,LL)=RDATA(NINT(XCP),NINT(YCP))

! CHG(3/28/02): try to read forecasts: RDATA has whole horizontal field!
! CHG(4/03/02): try identifying variables by label
!     IF(DREC(1)%VARB_ID(KK,LL) == VARB(NN))NT(KK)=NN
! Surface pressure
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'PRSS')WRITE(40,*)     &
     &   RDATA(1:NXGP,1:NYGP)
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'PRSS')WRITE(*,*)      &
     &   DREC(1)%VARB_ID(KVAR,LL)
! winds
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'UWND')WRITE(41,*)
!     :   RDATA(1:NXGP,1:NYGP)
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'VWND')WRITE(42,*)
!     :   RDATA(1:NXGP,1:NYGP)
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'WWND')WRITE(45,*)
!     :   RDATA(1:NXGP,1:NYGP)
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'UWND')WRITE(*,*)
!     :   DREC(1)%VARB_ID(KVAR,LL)
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'VWND')WRITE(*,*)
!     :   DREC(1)%VARB_ID(KVAR,LL)
!      IF(LL.GE.1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'WWND')WRITE(*,*)
!     :   DREC(1)%VARB_ID(KVAR,LL)
! Surface temperature
      IF ((.NOT.RAMSFLG) .AND. LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'T02M') WRITE(49,*) &
         RDATA(1:NXGP,1:NYGP)
      IF ((.NOT.RAMSFLG) .AND. LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'T02M') WRITE(*,*)  &
         DREC(1)%VARB_ID(KVAR,LL)
      IF ((RAMSFLG .OR. GRID(1)%MODEL_ID(1:2) == 'EC' .OR. GRID(1)%MODEL_ID == 'GFSN'  &
                   .OR. GRID(1)%MODEL_ID == 'ALAG')                                 &
           .AND.  LL == 1 .AND.  DREC(1)%VARB_ID(KVAR,LL) == 'T02M') WRITE(49,*) RDATA(1:NXGP,1:NYGP)
      IF ((RAMSFLG .OR. GRID(1)%MODEL_ID(1:2) == 'EC' .OR. GRID(1)%MODEL_ID == 'GFSN'  &
                   .OR. GRID(1)%MODEL_ID == 'ALAG')                                 &
           .AND. LL == 1 .AND.  DREC(1)%VARB_ID(KVAR,LL) == 'T02M') WRITE(*,*) DREC(1)%VARB_ID(KVAR,LL)
! radiation
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'DSWF')WRITE(48,*)     &
     &   RDATA(1:NXGP,1:NYGP)
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'DSWF')WRITE(*,*)      &
     &   DREC(1)%VARB_ID(KVAR,LL)
! Soil moisture
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'SOLW')WRITE(50,*)     &
     &   RDATA(1:NXGP,1:NYGP)
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'SOLW')WRITE(*,*)      &
     &   DREC(1)%VARB_ID(KVAR,LL)
! convective precipitation from ECMWF forecast
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'CPRC') THEN
         WRITE(51,*) RDATA(1:NXGP,1:NYGP)
         WRITE(*,*) DREC(1)%VARB_ID(KVAR,LL)
      ENDIF
! convective precipitation computed by the "convect"-module of GRIB2ARL
      IF (LL == 1 .AND. DREC(1)%VARB_ID(KVAR,LL) == 'PPRC') THEN
         WRITE(52,*) RDATA(1:NXGP,1:NYGP)
         WRITE(*,*) DREC(1)%VARB_ID(KVAR,LL)
      ENDIF

!==>some variables need to be saved for other calculations or require
!   additional units conversions not accounted by simple multiplication
!   factors defined in the sound subroutine

!           convert units of temperature to oC
            IF (VARB == 'TEMP'.OR.VARB == 'T02M'.OR.VARB == 'TMPS')      &
     &         VDATA(KVAR,LL)=VDATA(KVAR,LL)-273.16

!           save the surface pressure and terrain height for scaling
!           of the vertical coordinate system (height = sigma * scaling)
            IF (VARB == 'PRSS')SFCP=VDATA(KVAR,LL)
            IF (VARB == 'SHGT')SFCT=VDATA(KVAR,LL)

!           load winds for subsequent rotation to true
            IF (VARB == 'UWND'.OR.VARB == 'U10M')UTW(LL)=VDATA(KVAR,LL)
            IF (VARB == 'VWND'.OR.VARB == 'V10M')VTW(LL)=VDATA(KVAR,LL)

         END IF
      END DO

      RETURN
      END
