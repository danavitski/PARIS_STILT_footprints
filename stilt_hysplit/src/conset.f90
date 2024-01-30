!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONSET           CONcentration grid SET data entry
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION GRID SET IS THE DATA ENTRY FOR SAMPLING GRID LOCATION
!   LEVELS, AND TIME INFORMATION FOR START, STOP, AND DURATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Feb 1998 (RRD)
!                 07 May 1998 (RRD) - y2k mod
!                 20 Oct 1999 (RRD) - concentration grid dateline correction
!                 08 Nov 1999 (RRD) - concentration grid span correction
!
! USAGE:  CALL CONSET(ZMDL,NUMGRD,IUNIT,CGSIZE,OLAT,OLON,
!              IBYR,IBMO,IBDA,IBHR)
!   INPUT ARGUMENT LIST:
!     ZMDL      - real  model domain top in meters AGL
!     NUMGRD    - int   number of concentration grids
!     IUNIT     - int   unit number for input information
!     CGSIZE    - real  returns minimum grid spacing (km) for integration
!     OLAT,OLON - real  starting point
!     IBYR,IBMO,IBDA,IBHR  - int  starting time
!   OUTPUT ARGUMENT LIST:
!     COMMON GBLCON
!   INPUT FILES:
!     UNIT 5 or UNIT 30 - input CONTROL data file
!   OUTPUT FILES:
!     UNIT 31 - when input unit=5 then write STARTUP file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: conset.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

! JCL: add BACK as argument to make sure that the sampling
!     STOP time is earlier than the sampling START time when
!     BACK is T
      SUBROUTINE CONSET(ZMDL,NUMGRD,IUNIT,CGSIZE,OLAT,OLON,             &
     &   IBYR,IBMO,IBDA,IBHR,BACK)

!      SUBROUTINE CONSET(ZMDL,NUMGRD,IUNIT,CGSIZE,OLAT,OLON,
!     :   IBYR,IBMO,IBDA,IBHR)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

      LOGICAL BACK

!      COMMON /GBLCON/ CONC, DIRT

!          km/deg-lat        deg/radianT

      DATA YMPD/111.198323d0/, DGPR/57.295828d0/

!==>generic defaults

      NUMGRD=1
      CNLAT=OLAT
      CNLON=OLON
      CONC(1)%DELT_LAT=1.0
      CONC(1)%DELT_LON=1.0
      CONC(1)%LEVELS=1
      CONC(1)%HEIGHT(1)=50
      IF(MAXZP.GE.2)CONC(1)%HEIGHT(2)=100
      IF(MAXZP.GE.3)CONC(1)%HEIGHT(3)=200

      CONC(1)%START%YR   = IBYR
      CONC(1)%START%MO   = IBMO
      CONC(1)%START%DA   = IBDA
      CONC(1)%START%HR   = IBHR
      CONC(1)%START%MN   = 0
      CONC(1)%START%IC   = 0
      CONC(1)%START%MACC = 0

! JCL: not understand why original program had MOD statement
!           that gave weird years--e.g., 0 when IBYR=99
      CONC(1)%STOP%YR=IBYR
!      CONC(1).STOP.YR=MOD(IBYR+1,100)
      CONC(1)%STOP%MO=12
      CONC(1)%STOP%DA=31
      CONC(1)%STOP%HR=24
      CONC(1)%STOP%MN=60

      CONC(1)%SNAP=0
      CONC(1)%DELTA%HR=24
      CONC(1)%DELTA%MN=0

      CONC(1)%DIR='/main/sub/output/'
      CONC(1)%FILE='file_name'

!==>primary loop for the number of concentration grids

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Number of simultaneous concentration grids'
         WRITE(*,*)NUMGRD
      END IF
      READ(IUNIT,*)NUMGRD
      IF(IUNIT.EQ.5)WRITE(31,*)NUMGRD

!     test limits
      IF(NUMGRD.GT.MAXGRD)THEN
         WRITE(*,*)'ERROR conset: Number of grids exceed limit'
         STOP
      END IF

      DO KK=1,NUMGRD

         IF(KK.GT.1)CONC(KK)=CONC(KK-1)

!==>horizontal grid definition

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Enter the follow for grid #:',KK
            WRITE(*,*)'Center Latitude, Longitude'
            WRITE(*,*)CNLAT,CNLON
            READ(IUNIT,*)CNLAT,CNLON
         ELSE
            READ(IUNIT,*)CNLAT,CNLON
            IF(CNLAT.EQ.0.0.AND.CNLON.EQ.0.0)THEN
               CNLAT=OLAT
               CNLON=OLON
            END IF
         END IF
         IF(IUNIT.EQ.5)WRITE(31,*)CNLAT,CNLON

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Grid spacing (deg) Latitude, Longitude'
            WRITE(*,*)CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON
         END IF
         READ(IUNIT,*)CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON
         IF(IUNIT.EQ.5)WRITE(31,*)CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON

!        maximum concentration grid span
         SPLATM=MAXYP*CONC(KK)%DELT_LAT
         SPLONM=MAXXP*CONC(KK)%DELT_LON
         IF(IUNIT.EQ.5)THEN
            SPLAT=SPLATM
            SPLON=SPLONM
            WRITE(*,*)'Grid span (deg) Latitude, Longitide'
            WRITE(*,*)SPLAT,SPLON
            READ(IUNIT,*)SPLAT,SPLON
         ELSE
            READ(IUNIT,*)SPLAT,SPLON
         END IF
         IF(SPLAT.EQ.0.0) SPLAT=SPLATM
         IF(SPLON.EQ.0.0) SPLON=SPLONM
         SPLAT=DMIN1(SPLAT,SPLATM)
         SPLON=DMIN1(SPLON,SPLONM)
         IF(IUNIT.EQ.5)WRITE(31,*)SPLAT,SPLON

!        determine minimum grid spacing for all grids
!        km/degree-longitude
         XMPD=YMPD*DCOS(CNLAT/DGPR)
         CONC(KK)%SIZE=DMIN1(XMPD*CONC(KK)%DELT_LON,                    &
     &                     YMPD*CONC(KK)%DELT_LAT)
!        return minimum spacing over all grids
         CGSIZE=DMIN1(CONC(KK)%SIZE,CGSIZE)

!        determine number of grid points (not to exceed dimension)
         CONC(KK)%NUMB_LAT=MIN0(1+NINT(SPLAT/CONC(KK)%DELT_LAT),MAXYP)
         CONC(KK)%NUMB_LON=MIN0(1+NINT(SPLON/CONC(KK)%DELT_LON),MAXXP)

!        compute lower left corner
         CONC(KK)%X1Y1_LAT=CNLAT-SPLAT/2.0
         CONC(KK)%X1Y1_LON=CNLON-SPLON/2.0
!        dateline correction (RRD - 20/10/99)
         IF(CONC(KK)%X1Y1_LON.LT.-180.0)                                &
     &      CONC(KK)%X1Y1_LON=CONC(KK)%X1Y1_LON+360.0

!==>output file names and directory

         IF(IUNIT.EQ.5)THEN
! JCL:(6/1/00)'vf90' compiler doesn't like the '\.' characters
            WRITE(*,'(A,I2,A)')' Enter grid #',KK,' directory '
!            WRITE(*,'(A,I2,A)')' Enter grid #',KK,' directory (\...\)'
            WRITE(*,'(A)')CONC(KK)%DIR
         END IF
         READ(IUNIT,'(A)')CONC(KK)%DIR
         IF(IUNIT.EQ.5)WRITE(31,'(A)')CONC(KK)%DIR

         IF(IUNIT.EQ.5)THEN
            WRITE(*,'(A,I2,A)')' Enter grid #',KK,' file name (?????)'
            WRITE(*,'(A)')CONC(KK)%FILE
         END IF
         READ(IUNIT,'(A)')CONC(KK)%FILE
         IF(IUNIT.EQ.5)WRITE(31,'(A)')CONC(KK)%FILE

!==>vertical grid levels or spacing

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Number of vertical concentration levels'
            WRITE(*,*)CONC(KK)%LEVELS
         END IF
         READ(IUNIT,*)CONC(KK)%LEVELS
         IF(IUNIT.EQ.5)WRITE(31,*)CONC(KK)%LEVELS

!        test limits
         IF((CONC(KK)%LEVELS).GT.MAXZP)THEN
            WRITE(*,*)'ERROR conset: Number of levels exceed limit'
            STOP
         END IF

         NLVL=CONC(KK)%LEVELS
         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Height of each level m AGL'
            WRITE(*,*)(CONC(KK)%HEIGHT(KL),KL=1,NLVL)
         END IF
         READ(IUNIT,*)(CONC(KK)%HEIGHT(KL),KL=1,NLVL)
         IF(IUNIT.EQ.5)                                                 &
     &      WRITE(31,'(20I6)')(CONC(KK)%HEIGHT(KL),KL=1,NLVL)

         DO KL=1,NLVL
!           note level height "0" assumed to be for deposition
            IF(CONC(KK)%HEIGHT(KL).EQ.0.AND.KL.NE.1)THEN
               WRITE(*,*)'ERROR conset: height=0 should be level=1'
               STOP
            END IF
            IF(CONC(KK)%HEIGHT(KL).GT.ZMDL)THEN
               WRITE(*,*)'ERROR conset: level above model top'
               STOP
            END IF
         END DO

!==>set sampling time intervals

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Sampling start time: year month day hour minute'
            WRITE(*,*)CONC(KK)%START%YR, CONC(KK)%START%MO,             &
     &         CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN
            READ(*,*) CONC(KK)%START%YR, CONC(KK)%START%MO,             &
     &         CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN
         ELSE
            READ(IUNIT,*)IYR,IMO,IDA,IHR,IMN
            IF(IMO.EQ.0)THEN
               CONC(KK)%START%DA=IBDA+IDA
               CONC(KK)%START%HR=IBHR+IHR
            ELSE
               CONC(KK)%START%YR=IYR
               CONC(KK)%START%MO=IMO
               CONC(KK)%START%DA=IDA
               CONC(KK)%START%HR=IHR
               CONC(KK)%START%MN=IMN
            END IF
         END IF

         IF(IUNIT.EQ.5)WRITE(31,*)                                      &
     &      CONC(KK)%START%YR, CONC(KK)%START%MO,                       &
     &      CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Sampling stop time: year month day hour minute'
            WRITE(*,*)CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,               &
     &         CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN
            READ(*,*) CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,               &
     &         CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN
         ELSE
            READ(IUNIT,*)IYR,IMO,IDA,IHR,IMN

            IF(IYR+IMO+IDA+IHR+IMN.NE.0)THEN
               CONC(KK)%STOP%YR=IYR
               CONC(KK)%STOP%MO=IMO
               CONC(KK)%STOP%DA=IDA
               CONC(KK)%STOP%HR=IHR
               CONC(KK)%STOP%MN=IMN
! JCL: if all zeros in CONTROL file, then assign default values that
!      are before the START%values when BACK is T (Jan 1st of PREVIOUS YR)
            ELSEIF(BACK)THEN
               CONC(KK)%STOP%YR=IBYR-1
               CONC(KK)%STOP%MO=1
               CONC(KK)%STOP%DA=1
               CONC(KK)%STOP%HR=1
               CONC(KK)%STOP%MN=1
            END IF
         END IF

         IF(IUNIT.EQ.5)WRITE(31,*)                                      &
     &      CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,                         &
     &      CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Sampling interval: type hour minute'
            WRITE(*,*)CONC(KK)%SNAP, CONC(KK)%DELTA%HR,                 &
     &         CONC(KK)%DELTA%MN
         END IF
         READ(IUNIT,*)CONC(KK)%SNAP, CONC(KK)%DELTA%HR,                 &
     &      CONC(KK)%DELTA%MN

         IF(CONC(KK)%DELTA%MN.NE.0)                                     &
     &      WRITE(*,*)'WARNING: conset - Non-zero minutes'
         IF(IUNIT.EQ.5)WRITE(31,*)                                      &
     &      CONC(KK)%SNAP, CONC(KK)%DELTA%HR, CONC(KK)%DELTA%MN

!==>convert dates to accumulated time

         CALL TM2MIN(CONC(KK)%START%YR, CONC(KK)%START%MO,              &
     &      CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN,    &
     &      CONC(KK)%START%MACC)

         CALL TM2MIN(CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,                &
     &      CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN,       &
     &      CONC(KK)%STOP%MACC)

!        convert back to date to fix date+delta errors
         CALL TM2DAY(CONC(KK)%START%MACC,                               &
     &      CONC(KK)%START%YR, CONC(KK)%START%MO,                       &
     &      CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN)

         CALL TM2DAY(CONC(KK)%STOP%MACC,                                &
     &      CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,                         &
     &      CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN)

!==>set remaining variables

!        sampling intervals in years, months, days not used
         CONC(KK)%DELTA%YR=0
         CONC(KK)%DELTA%MO=0
         CONC(KK)%DELTA%DA=0

!        convert sampling interval to minutes
         CONC(KK)%DELTA%MACC=60*CONC(KK)%DELTA%HR+CONC(KK)%DELTA%MN

!        save sampling start time in conc output file marker variable
         CONC(KK)%NOW=CONC(KK)%START

      END DO

! JCL:
      WRITE(45,*)'CONC(1)%START%YR:  ',CONC(1)%START%YR
      WRITE(45,*)'CONC(1)%START%MO:  ',CONC(1)%START%MO
      WRITE(45,*)'CONC(1)%START%DA:  ',CONC(1)%START%DA
      WRITE(45,*)'CONC(1)%START%HR:  ',CONC(1)%START%HR
      WRITE(45,*)'CONC(1)%START%MN:  ',CONC(1)%START%MN
      WRITE(45,*)'CONC(1)%STOP%YR:  ',CONC(1)%STOP%YR
      WRITE(45,*)'CONC(1)%STOP%MO:  ',CONC(1)%STOP%MO
      WRITE(45,*)'CONC(1)%STOP%DA:  ',CONC(1)%STOP%DA
      WRITE(45,*)'CONC(1)%STOP%HR:  ',CONC(1)%STOP%HR
      WRITE(45,*)'CONC(1)%STOP%MN:  ',CONC(1)%STOP%MN

      WRITE(45,*)'CONC(1)%START%MACC:  ',CONC(1)%START%MACC
      WRITE(45,*)'CONC(1)%STOP%MACC:  ',CONC(1)%STOP%MACC


      RETURN
      END
