!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SFCINP           SurFaCe data INPut reads lat/lon files
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SURFACE DATA INPUT READS A LAT/LON BASED SFC CHARACTERISTICS FILE
!   CURRENT COMPILED VERSION ASSUMES A FILE WITH 1-DEG RESOLUTION WITH
!   THE FIRST RECORD STARTING AT THE NORTHWEST CORNER (CNTR 179.5W,89.5N)
!   RETURNS THE ROUGHNESS LENGTH AND LAND-USE AT SPECIFIED POINT.
!
!   Land-Use Values are defined as follows:
!     1-urban           2-agricultural    3-dry range 4-deciduous
!     5-coniferous      6-mixed forest    7-water           8-desert
!     9-wetlands  10-mixed range    11-rocky
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Apr 1997 (RRD)
!                 14 Jun 1999 (RRD) - generalized data file input
!                 21 Oct 1999 (RRD) - LREC added to save statement
!
! USAGE:  CALL SFCINP(CLAT,CLON,ZNOT,LUSE)
!   INPUT ARGUMENT LIST:
!     CLAT,CLON   - real      Latitude/Longitude of required point
!   OUTPUT ARGUMENT LIST:
!     ZNOT  - real      Aerodynamic rougness length (m)
!     LUSE  - int Land-use categories (1-11)
!   INPUT FILES:
!     UNIT 60 - landuse categories (LANDUSE.ASC)
!     UNIT 62 - roughness length file (ROUGLEN.ASC)
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: sfcinp.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

      SUBROUTINE SFCINP(CLAT,CLON,ZNOT,LUSE)

      IMPLICIT REAL*8 (A-H,O-Z)

!==>maximum file configuration set to 0.25 degrees (360 x 4 data points)
!   record length number of data points x 4 bytes

      PARAMETER (MLON=1440)
      PARAMETER (MREC=MLON*4+1)

!     dummy variable to contain data of maximum length MREC
      CHARACTER LABEL*5761, FDIR*80
      REAL*8 ROUG(MLON)
      INTEGER LAND(MLON)

!     flag for files
      LOGICAL TFILE, QFILE, RFILE

      SAVE QFILE, TFILE, RFILE, KREC, ROUG, LAND, ALATB, ALONL,         &
     &     DLAT, DLON, NLAT, NLON, LUSE0, ZNOT0, LREC

      DATA RFILE/.FALSE./, QFILE/.FALSE./, TFILE/.FALSE./, KREC/0/

!==>open the files only once

      IF(.NOT.QFILE)THEN

!        the configuration file

         INQUIRE(FILE='../bdyfiles/ASCDATA.CFG',EXIST=QFILE)
         IF(QFILE)THEN
            OPEN(60,FILE='../bdyfiles/ASCDATA.CFG')
         ELSE
            INQUIRE(FILE='ASCDATA.CFG',EXIST=QFILE)
            IF(QFILE)OPEN(60,FILE='ASCDATA.CFG')
         END IF

         IF(QFILE)THEN
            WRITE(30,*)'NOTICE sfcinp: reading ASCDATA.CFG'
            READ(60,*)ALATB,ALONL
            READ(60,*)DLAT,DLON
            READ(60,*)NLAT,NLON
            READ(60,*)LUSE0
            READ(60,*)ZNOT0
            READ(60,*)FDIR
            CLOSE(60)
         ELSE
            WRITE(30,*)'NOTICE sfcinp: no ASCDATA.CFG, using default'
!           lower left corner of the file
            ALATB=-90.0
            ALONL=-180.0
!           incremental spacing
            DLAT=1.0
            DLON=1.0
!           number of points in each direction
            NLAT=180
            NLON=360
!           content default values
            LUSE0=2
            ZNOT0=0.2
!           data file directory
            FDIR='../bdyfiles/'
         END IF

!        record structure

         KLEN=INDEX(FDIR,' ')-1
         LREC=1+NLON*4
         IF(NLON.GT.MLON.OR.LREC.GT.MREC)THEN
            WRITE(30,*)'ERROR sfcinp: input data record length exceeded'
            WRITE(30,*)'Requested   : ',NLON,LREC
            WRITE(30,*)'Compiled    : ',MLON,MREC
            STOP
         END IF

!        land use file

         INQUIRE(FILE=FDIR(1:KLEN)//'LANDUSE.ASC',EXIST=QFILE)
         IF(QFILE)THEN
            OPEN(60,FILE=FDIR(1:KLEN)//'LANDUSE.ASC',ACCESS='DIRECT',   &
     &           RECL=LREC,FORM='UNFORMATTED')
            TFILE=.TRUE.
         ELSE
            WRITE(30,*)'NOTICE sfcinp: LANDUSE.ASC file not found'
            WRITE(30,*)'On default directory:',FDIR(1:KLEN)
            WRITE(30,*)'Using default value category 2 agricultural'
         END IF

!        roughness length file

         INQUIRE(FILE=FDIR(1:KLEN)//'ROUGLEN.ASC',EXIST=QFILE)
         IF(QFILE)THEN
            OPEN(62,FILE=FDIR(1:KLEN)//'ROUGLEN.ASC',ACCESS='DIRECT',   &
     &         RECL=LREC,FORM='UNFORMATTED')
            RFILE=.TRUE.
         ELSE
            WRITE(30,*)'NOTICE sfcinp: ROUGLEN.ASC file not found'
            WRITE(30,*)'On default directory:',FDIR(1:KLEN)
            WRITE(30,*)'Using default value 0.20 m'
         END IF

!        set logical file test so that files are not opened again
         QFILE=.TRUE.

      END IF

!==>load data from file as required

      IF(RFILE.OR.TFILE)THEN

!        determine required record number based upon latitude
         JREC=NLAT-INT((CLAT-ALATB)/DLAT)
         JREC=MIN0(MAX0(1,JREC),NLAT)

!        compute element number from longitude
         KK=INT((CLON-ALONL)/DLON)+1
         KK=MIN0(MAX0(1,KK),NLON)

         IF(JREC.NE.KREC)THEN
!           load new data into buffers if record different from previous
            KREC=JREC

            IF(TFILE)THEN
               READ(60,REC=JREC)LABEL(1:LREC)
               READ(LABEL,'(4(360I4))') (LAND(K),K=1,NLON)
               LUSE=LAND(KK)
            END IF

            IF(RFILE)THEN
               READ(62,REC=JREC)LABEL(1:LREC)
               READ(LABEL,'(4(360F4.0))') (ROUG(K),K=1,NLON)
!              convert from cm to m
               ZNOT=ROUG(KK)*0.01
            END IF

         ELSE
!           existing data in buffer OK
            IF(TFILE)LUSE=LAND(KK)
            IF(RFILE)ZNOT=ROUG(KK)*0.01
         END IF

      ELSE
!        default values
         LUSE=LUSE0
         ZNOT=ZNOT0

      END IF

      RETURN
      END
