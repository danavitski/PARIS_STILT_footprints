!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSDAT           EMiSsion DATa read of special file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INPUT - READS A LAT/LON BASED EMISSION INVENTORY FILE.
!   THE EMISSION ARRAY BASED UPON A TWO-POINT SOURCE LIMIT DEFINED IN
!   MODEL INPUT SECTION. THE ROUTINE IS CALLED ONLY ONCE.
!   THE FILE IS ONE RECORD PER EMISSION POINT WITH EACH EMISSION POINT
!   DEFINED AT A LAT/LON POINT. EMISSIONS ARE SUMMED INTO A GRID.
!   EMISSION FILE PROVIDED BY AERONOMY LABORATORY.  NATURAL EMISSIONS
!   OF ISOPRENE ARE COMPUTED FOR ALL FOREST LAND-USE CATEGORIES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Feb 1998 (RRD)
!                 22 Dec 1998 (RRD) - generalized subroutine
!
! USAGE:  CALL EMSDAT(NLOC,NUMTYP,QFILE)
!   INPUT ARGUMENT LIST:
!     NLOC      - int   total number of source locations
!     NUMTYP    - int   number of pollutant types
!   OUTPUT ARGUMENT LIST:
!     QFILE - log flag to indicate area source emission file
!   INPUT FILES:
!     UNIT 50 defines emission input data file
!   OUTPUT FILES:
!     UNIT 30 for diagnostic messages to MESSAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: emsdat.f90,v 1.5 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE EMSDAT(NLOC,NUMTYP,QFILE)

      use module_defspot
      use module_defgrid

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology file and grid
!      INCLUDE 'DEFGRID.INC'
!     multiple source information
!      INCLUDE 'DEFSPOT.INC'

!     number of emission grid points
      PARAMETER (MQLAT=150, MQLON=150, MQVAL=4, MQHRS=24)

      REAL*8 QAREA(MQLON,MQLAT,MQVAL,MQHRS)
      CHARACTER POLID(MQVAL)*4, FNAME*80
      REAL*8 QVAL(MQHRS)

!     diagnostic variable
      REAL*8 QSUM(MQVAL)

!     emission file test, record test
      LOGICAL QFILE, SKIP

!      COMMON /GBLSPT/ SPOT
!      COMMON /GBLGRD/ GRID, DREC, FILE
      COMMON /QARRAY/ KLAT, KLON, DLAT, DLON, QAREA, POLID
      SAVE /QARRAY/

!==>input the emissions from the inventory file
!   flag set to point source emission if file is not found

      FNAME='emission.txt'
      INQUIRE(FILE=FNAME,EXIST=QFILE)
      IF(QFILE)THEN
         OPEN(50,FILE=FNAME)
      ELSE
         RETURN
      END IF

!==>check inputs and dimensions for consistency

      IF(NUMTYP.GT.MQVAL.OR.NUMTYP.GT.MAXDIM)THEN
         WRITE(30,*)'ERROR emsdat: # pollutants exceeds dimension'
         WRITE(30,*)'      MQVAL = ',mqval,'    MAXDIM = ',maxdim
         STOP
      END IF

      IF(NLOC.NE.2)THEN
         WRITE(30,*)'ERROR emsdat: 2 source locations required'
         WRITE(30,*)'      to define grid array lat/lon limits'
         STOP
      END IF

!==>zero out the emissions array

      DO L=1,MQHRS
      DO K=1,MQVAL
         DO J=1,MQLAT
         DO I=1,MQLON
            QAREA(I,J,K,L)=0.0
         END DO
         END DO
      END DO
      END DO

      DO K=1,MQVAL
!        diagnostic variable
         QSUM(K)=0.0
         POLID(K)='    '
      END DO

!==>set emission file default values from index record (#1)
!   NQVAL - the number of pollutants in file
!   UNITS - defines conversion from file units to kg/hour
!   DLAT  - selected resolution for input file summation, which should
!   DLON    be set to some a value equal to or coarser than file
!   POLID - definitions must match pollutants defined in control file

      READ(50,'(I4,3F10.4,10A4)')                                       &
     &     NQVAL,UNITS,DLAT,DLON,(POLID(K),K=1,NQVAL)

!==>internal grid limits according to selected resolution

!     average grid emission cell area (km^2)
      ALAT=0.5*(SPOT(1)%OLAT+SPOT(2)%OLAT)
      ALON=0.5*(SPOT(1)%OLON+SPOT(2)%OLON)
      AREA=(DLAT*111.0)*(DLON*111.0*COS(ALAT/57.3))

!     maximum number of grid points requested
      KLAT=(SPOT(2)%OLAT-SPOT(1)%OLAT)/DLAT+1
      KLON=(SPOT(2)%OLON-SPOT(1)%OLON)/DLON+1

      IF(KLAT.GT.MQLAT.OR.KLON.GT.MQLON)THEN
         WRITE(30,*)'ERROR emsdat: requested emission area too large'
         WRITE(30,*)'Lat/Lon pnts: ',KLAT,KLON
         WRITE(30,*)'Dimensions  : ',MQLAT,MQLON
         STOP
      END IF

!==>loop thru entire file

  100 READ(50,'(2I4,2F10.4)',END=900)II,JJ,XLON,XLAT

!     check if record is within selected limits
      SKIP=.FALSE.
      IF(XLAT.LT.SPOT(1)%OLAT.OR.XLAT+DLAT.GT.SPOT(2)%OLAT)SKIP=.TRUE.
      IF(XLON.LT.SPOT(1)%OLON.OR.XLON+DLON.GT.SPOT(2)%OLON)SKIP=.TRUE.

!     compute index on internal emission grid
!     when SW corner falls within cell then accumulate in that cell
      II=(XLON-SPOT(1)%OLON)/DLON+1
      JJ=(XLAT-SPOT(1)%OLAT)/DLAT+1

!     check limits again
      IF(II.LT.1.OR.II.GT.MQLON.OR.JJ.LT.1.OR.JJ.GT.MQLAT)SKIP=.TRUE.

      DO NN=1,NQVAL
         READ(50,'(12E10.3)')(QVAL(IH),IH=1,MQHRS)

         IF(.NOT.SKIP)THEN
!           fill array after converting to standard units (/hour)
            DO IH=1,MQHRS

!              convert emissions to kg/hr
               QAREA(II,JJ,NN,IH)=QVAL(IH)*UNITS

!              sum average hourly rate for diagnostics
               QSUM(NN)=QSUM(NN)+QAREA(II,JJ,NN,IH)/MQHRS

            END DO
         END IF

      END DO
      GO TO 100

!==>diagnostic dump of generated emissions

  900 KPTS=0
      DO J=1,MQLAT
      DO I=1,MQLON
         SKIP=.FALSE.
         DO L=1,MQHRS
         DO K=1,MQVAL
!           set flag if any hour or pollutant is non-zero
            IF(QAREA(I,J,K,L).GT.0.0)SKIP=.TRUE.
         END DO
         END DO
!        count the number of non-zero grid points
         IF(SKIP)KPTS=KPTS+1
      END DO
      END DO

!==>terminate and print message according to status

      IF(KPTS.EQ.0)THEN
         WRITE(30,*)'ERROR emsdat: lat/lon emission selection'
         WRITE(30,*)'There are no emissions in the region'
         STOP
      ELSE
         WRITE(30,*)'NOTICE emsdat: gridded emissions initialized'
         WRITE(30,'(5A10)')'Kg/hr',(POLID(N),N=1,NQVAL)
         WRITE(30,'(10X,4E10.2)') (QSUM(N),N=1,NQVAL)
         WRITE(30,*)'Number of emission lat/lon points: ',KLAT,KLON
         WRITE(30,*)'Total number of non-zero points  : ',KPTS
         WRITE(30,*)'Emission grid cell area (km^2)   : ',AREA
      END IF

      RETURN
      END
