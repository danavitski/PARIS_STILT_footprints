!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METOLD           METeorological OLD data file analysis
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL OLD DATA FILE DETERMINES THE GRID OF OLD STYLE
!   FILE AND SETS THE STRUCTURE.  AN OLD STYLE FILE IS ONE THAT
!   DOES NOT HAVE AN INDEX RECORD AT THE BEGINNING OF EACH TIME
!   PERIOD.  THESE FILES MUST BE COMPLETELY IDENTIFIED IN THE
!   DATA STATEMENTS OF THIS SUBROUTINE.  PROPER FUNCTIONING OF
!   THIS ROUTINE REQUIRES THAT A UNIQUE GRID NUMBER BE WRITTEN IN
!   THE 50 BYTE ASCII HEADER PORTION OF EACH DATA RECORD.
!
! PROGRAM HISTORY LOG:
!   last revision: 23 Dec 1998 (RRD)
!                  21 Oct 1999 (RRD) - move structure variables out of read
!
! USAGE:  CALL METOLD(KG,KT,KLEN)
!   INPUT ARGUMENT LIST:
!     KG        - int   sequential grid identfication number
!     KT        - int   data time period identification (1 or 2)
!     KLEN      - int   length of directory string
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
! $Id: metold.f90,v 1.5 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE METOLD(KG,KT,KLEN)

      use module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'

!     real and integer grid information for three grids
      REAL*8   RGRID(12,3),IGRID(4,3)
!     vertical information
      REAL*8   VGRID(15,3)
!     grid number index, number variables
      INTEGER   KDX(0:9), NVARB(15,3)

      CHARACTER KVAR*4, LABEL*50, MTYPE(3)*4, FNAME*80

!      COMMON /GBLGRD/ GRID, DREC, FILE

!     associate grid number with array index (only 0-9)
!     grid num 0 1 2 3 4 5 6 7 8 9
      DATA KDX/0,0,2,3,0,0,1,0,0,0/

!     meteo model data names
      DATA MTYPE/'NGM ','MRFN','MRFS'/

!     grid location parameters
      DATA RGRID/                                                       &
     &   90.,0., 60.,-105.,182.9, 0., 90., 13.25,42.75, 90.,0.,0.,      &
     &   90.,0., 60., -80.,381.0, 0., 90., 33.00,33.00, 90.,0.,0.,      &
     &  -90.,0.,-60., -80.,381.0, 0.,-90., 33.00,33.00,-90.,0.,0./

!     vertical structure, x and y dimensions, number levels
      DATA IGRID/                                                       &
     &   1,  33, 28, 11,                                                &
     &   2,  65, 65, 13,                                                &
     &   2,  65, 65, 13/

!     number of variables at each level
      DATA NVARB/                                                       &
     &   11, 10*5, 4*0,                                                 &
     &   03,  6*6, 4*5, 2*4, 2*0,                                       &
     &   03,  6*6, 4*5, 2*4, 2*0/

!     coordinates of vertical levels
      DATA VGRID/                                                       &
     &   1.,.9823,.9432,.8967,.8437,.7848,.721,.6531,.582,              &
     &      .5086,.4341,4*0.0,                                          &
     &   0.,1000.,850.,700.,500.,400.,300.,250.,100.,150.,              &
     &      100.,70.,50.,2*0.0,                                         &
     &   0.,1000.,850.,700.,500.,400.,300.,250.,100.,150.,              &
     &      100.,70.,50.,2*0.0/

!     simplify structure variables for IO statements
      LABEL=FILE(KG,KT)%DIR
      KUNIT=FILE(KG,KT)%KUNIT
      FNAME=LABEL(1:KLEN)//FILE(KG,KT)%METEO

      OPEN(KUNIT,FILE=FNAME,RECL=50,ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')

!     decode the fixed portion of the index record
      READ(KUNIT,REC=1)LABEL
      READ(LABEL,'(5I2,2X,I2,A4)')                                      &
     &   FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,                    &
     &   FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR,                    &
     &   FILE(KG,KT)%FIRST%IC, MG, KVAR
      FILE(KG,KT)%FIRST%MN=0
      CLOSE (KUNIT)

      IF(KDX(MG).EQ.0)THEN
         WRITE(*,*)'ERROR metold: meteo grid not pre-defined'
         STOP
      END IF

      IF(KVAR.EQ.'NULL'.OR.KVAR.EQ.'    ')THEN
         WRITE(*,*)'ERROR metold: meteo data file first record'
         WRITE(*,*)'contains NULL or BLANK field for variable ID'
         STOP
      END IF

      MG1=MG
!     over-ride grid and set for SH when (offset = olat < 0)
      IF((DREC(KG)%OFFSET).LT.0.0)MG1=3

!     set the spatial grid information
      GRID(KG)%MODEL_ID=MTYPE(KDX(MG1))
      GRID(KG)%POLE_LAT=RGRID(1,KDX(MG1))
      GRID(KG)%POLE_LON=RGRID(2,KDX(MG1))
      GRID(KG)%REF_LAT= RGRID(3,KDX(MG1))
      GRID(KG)%REF_LON= RGRID(4,KDX(MG1))
      GRID(KG)%SIZE=    RGRID(5,KDX(MG1))
      GRID(KG)%ORIENT=  RGRID(6,KDX(MG1))
      GRID(KG)%TANG_LAT=RGRID(7,KDX(MG1))
      GRID(KG)%SYNC_XP= RGRID(8,KDX(MG1))
      GRID(KG)%SYNC_YP= RGRID(9,KDX(MG1))
      GRID(KG)%SYNC_LAT=RGRID(10,KDX(MG1))
      GRID(KG)%SYNC_LON=RGRID(11,KDX(MG1))

      GRID(KG)%NUMBER=MG1
      DREC(KG)%Z_FLAG=IGRID(1,KDX(MG1))
      GRID(KG)%NX=IGRID(2,KDX(MG1))
      GRID(KG)%NY=IGRID(3,KDX(MG1))
      GRID(KG)%NZ=IGRID(4,KDX(MG1))

!     check dimensions
      IF((GRID(KG)%NZ).GT.MLVL)THEN
         WRITE(*,*)'ERROR metset: exceeding structure dimension MLVL'
         STOP
      END IF

!     fill remainder of structure
      NLVL=GRID(KG)%NZ
      DO L=1,NLVL
         DREC(KG)%HEIGHT(L)=VGRID(L,KDX(MG))
         DREC(KG)%NUM_VARB(L)=NVARB(L,KDX(MG))
      END DO

!     reopen with proper length
      FILE(KG,KT)%REC_LEN = GRID(KG)%NX*GRID(KG)%NY+50
      LREC = FILE(KG,KT)%REC_LEN
      OPEN(KUNIT,FILE=FNAME,RECL=LREC,ACCESS='DIRECT',                  &
     &     FORM='UNFORMATTED',ACTION='READ')

!     loop through and decode the remainder
      NREC=0
      NLVL=GRID(KG)%NZ
      DO L=1,NLVL
         NVAR=DREC(KG)%NUM_VARB(L)
         DO K=1,NVAR
            NREC=NREC+1
            READ(KUNIT,REC=NREC)LABEL
            READ(LABEL,'(14X,A4)')KVAR
            DREC(KG)%VARB_ID(K,L)=KVAR
            DREC(KG)%CHK_SUM(K,L)=0
         END DO
      END DO
      DREC(KG)%REC_PER=NREC

!     skip to the next time period index record to find the time
!     interval between data periods (minutes)
      NREC=NREC+1
      READ(KUNIT,REC=NREC)LABEL
      READ(LABEL,'(5I2,2X,I2,A4)')                                      &
     &   FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                      &
     &   FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                      &
     &   FILE(KG,KT)%LAST%IC, MG2, KVAR
      FILE(KG,KT)%LAST%MN=0

!     check additional hemisphere (special SH data available)
      IF(MG2.NE.MG)THEN
!        skip to new time past SH data
         NREC=NREC+DREC(KG)%REC_PER
         READ(KUNIT,REC=NREC)LABEL
         READ(LABEL,'(5I2)')                                            &
     &      FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                   &
     &      FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                   &
     &      FILE(KG,KT)%LAST%IC
!        enable data skipping
         DREC(KG)%REC_PER=DREC(KG)%REC_PER*2
!        shm start request requires positioning offset
         IF((DREC(KG)%OFFSET).LT.0)DREC(KG)%OFFSET=DREC(KG)%REC_PER/2
      ELSE
!        no multiple grids per time period
         DREC(KG)%OFFSET=0
      END IF

!     compute the minute difference first for hires files
      DREC(KG)%DELTA=FILE(KG,KT)%LAST%MN-FILE(KG,KT)%FIRST%MN
!     add the hour difference in minutes
      DREC(KG)%DELTA=DREC(KG)%DELTA+                                    &
     &   (FILE(KG,KT)%LAST%HR-FILE(KG,KT)%FIRST%HR)*60

!     continue on to find the last record
      NREC=NREC+DREC(KG)%REC_PER
      DO WHILE (NREC.GE.1)
         READ(KUNIT,REC=NREC,ERR=500)LABEL
         READ(LABEL,'(5I2,4X,A4)',ERR=500)                              &
     &      FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                   &
     &      FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                   &
     &      FILE(KG,KT)%LAST%IC, KVAR
         FILE(KG,KT)%LAST%MN=0
         NREC=NREC+DREC(KG)%REC_PER

      END DO
  500 FILE(KG,KT)%ENDREC=NREC-1

!     set the year's accumulated clock time in minutes
!     for the first time period on the file
      CALL TM2MIN(FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,           &
     &            FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR,           &
     &            FILE(KG,KT)%FIRST%MN, FILE(KG,KT)%FIRST%MACC)

!     for the last time period on the file
      CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,             &
     &            FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,             &
     &            FILE(KG,KT)%LAST%MN, FILE(KG,KT)%LAST%MACC)

      RETURN
      END
