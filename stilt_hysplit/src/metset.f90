!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METSET           METeorological data structure SET
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL STRUCTURE SET READS AND UNPACKS THE METEOROLOGICAL
!   INDEX RECORD AND INITIALIZES THE FILE STRUCTURE DEFINITIONS.
!
! PROGRAM HISTORY LOG:
! LAST%REVISION: 23 Dec 1998 (RRD)
!                  20 Jan 1999 (RRD) - added status= on open statement
!                  09 Feb 1999 (RRD) - test for max dir string length
!                  29 Oct 1999 (RRD) - moved structure variables out of read
!                                    - fixed problem with single time files
!
! USAGE:  CALL METSET(HEADER,KG,KT)
!   INPUT ARGUMENT LIST:
!     HEADER      - char      extended index record header (length unknown)
!     KG    - int sequential grid identfication number
!     KT    - int data time period identification (1 or 2)
!   OUTPUT ARGUMENT LIST:
!     COMMON GBLGRD
!   INPUT FILES:
!     UNIT 10,12,14,etc depending on number of input files defined
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metset.f90,v 1.12 2007-05-03 13:09:13 skoerner Exp $
!
!$$$

      SUBROUTINE METSET(HEADER,KG,KT)

      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'

      CHARACTER KVAR*4, LABEL*50, HEADER*10999, FNAME*80
      LOGICAL FTEST
      integer :: ierr


!-------------------------------------------------------------------------------
!     test for meteo file existence
      KLEN=INDEX(FILE(KG,KT)%DIR,' ')-1
      IF(KLEN.LE.0)THEN
         WRITE(*,*)'Directory: ',FILE(KG,KT)%DIR
         WRITE(*,*)'String length exceeds 40 character maximum'
         STOP
      END IF
      LABEL=FILE(KG,KT)%DIR
      FNAME=LABEL(1:KLEN)//FILE(KG,KT)%METEO
      INQUIRE(FILE=FNAME,EXIST=FTEST)
      IF(.NOT.FTEST)THEN
         WRITE(*,*)'Unable to find file: ',FILE(KG,KT)%METEO
         WRITE(*,*)'On local directory : ',FILE(KG,KT)%DIR(1:KLEN)
         WRITE(*,*)'Check input CONTROL file for correct values'
         STOP
      END IF

!     set IO unit number
      KUNIT=FILE(KG,KT)%KUNIT

!     open file to decode the standard label (50) plus the
!     fixed portion (108) of the extended header
      OPEN(KUNIT, FILE=FNAME,STATUS='OLD',                              &
     &   RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')

!     decode the standard portion of the index record
      READ(KUNIT,REC=1)LABEL,HEADER(1:108)
      READ(LABEL,'(5I2,2X,I2,A4)')                                      &
     &   FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,                    &
     &   FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR,                    &
     &   FILE(KG,KT)%FIRST%IC, IGRID, KVAR

      WRITE(*,*)'LABEL in METSET:',LABEL
      WRITE(*,*)'HEADER in METSET:',HEADER(1:108)

      IF(KVAR.NE.'INDX')THEN
         WRITE(*,*)'WARNING metset: Old format meteo data grid ',KG
         CLOSE(KUNIT)
         CALL METOLD(KG,KT,KLEN)
         DREC(KG)%TYPE=1
         RETURN
      END IF
      DREC(KG)%TYPE=2


!     decode extended portion of the header
      READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)')                 &
     &   GRID(KG)%MODEL_ID, ICX,                FILE(KG,KT)%FIRST%MN,   &
     &   GRID(KG)%POLE_LAT, GRID(KG)%POLE_LON,  GRID(KG)%REF_LAT,       &
     &   GRID(KG)%REF_LON,  GRID(KG)%SIZE,      GRID(KG)%ORIENT,        &
     &   GRID(KG)%TANG_LAT, GRID(KG)%SYNC_XP,   GRID(KG)%SYNC_YP,       &
     &   GRID(KG)%SYNC_LAT, GRID(KG)%SYNC_LON,  GRID(KG)%DUMMY,         &
     &   GRID(KG)%NX,       GRID(KG)%NY,        GRID(KG)%NZ,            &
     &   DREC(KG)%Z_FLAG,   LENH

  ! automatic detect values for lenh > 9999;  these values are stored as negative numbers
  if (lenh < 0) then
         WRITE(*,*) 'lenh:',lenh
    lenh=-lenh + 10000
  endif

!     check number of levels against structure dimension
      IF((GRID(KG)%NZ).GT.MLVL)THEN
         WRITE(*,*) 'ERROR metset: exceeding structure dimension MLVL'
         WRITE(*,*) 'Recompile after increasing MLVL in module_defsize. Stop.'
         STOP
      END IF

!     grid id variable (needed for old data sets)
      GRID(KG)%NUMBER=IGRID

!     close file and reopen with proper length
      CLOSE (KUNIT)

      !  change to 16 bit for RAMS and ECMWF (2 bytes rqd)
      IF (GRID(KG)%MODEL_ID  == 'RAMS'  .OR.  GRID(KG)%MODEL_ID(1:3)  == 'ECX' .or. &
           & GRID(KG)%MODEL_ID  == 'DWRF') THEN
         LDAT=GRID(KG)%NX*GRID(KG)%NY*2
      ELSE
         LDAT=GRID(KG)%NX*GRID(KG)%NY
      END IF
      WRITE(*,*)'LDAT:  ',LDAT


      FILE(KG,KT)%REC_LEN = LDAT+50
      LREC = FILE(KG,KT)%REC_LEN
      OPEN(KUNIT, FILE=FNAME,STATUS='OLD',                              &
     &   RECL=LREC,ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')

!     determine number of index records
      NNDX=LENH/LDAT+1

!     read extended character string over multiple index records
      NHL1=1
      KREC=1
      DO N=1,NNDX
         NHL2=NHL1+LDAT-1
         IF(N.EQ.NNDX)NHL2=NHL1+(LENH-(NNDX-1)*LDAT)-1
         READ(KUNIT,REC=KREC)LABEL,HEADER(NHL1:NHL2)
         WRITE(*,*)'metset nhl',NHL1,NHL2,N,NNDX,LENH
         KREC=KREC+1
         NHL1=NHL2+1
      END DO

!     loop through and decode the remainder of the index string
      KOL=109
      NREC=NNDX
      NLVL=GRID(KG)%NZ
      DO L=1,NLVL
         READ(HEADER(KOL:KOL+7),'(F6.2,I2)') DREC(KG)%HEIGHT(L), DREC(KG)%NUM_VARB(L)
         KOL=KOL+8
         NVAR=DREC(KG)%NUM_VARB(L)
         DO K=1,NVAR
! CHG (10/01/03) change from 8 to 16 bit
            IF(GRID(KG)%MODEL_ID /= 'RAMS' .AND. GRID(KG)%MODEL_ID(1:3) /= 'ECX' .AND. &
                 & GRID(KG)%MODEL_ID  /= 'DWRF')THEN
               READ(HEADER(KOL:KOL+7),'(A4,I3)',iostat=ierr) DREC(KG)%VARB_ID(K,L), DREC(KG)%CHK_SUM(K,L)
               if (ierr .ne. 0) then
                  DREC(KG)%CHK_SUM(K,L) = 0
                  READ(HEADER(KOL:KOL+7),'(A4,I3)',iostat=ierr) DREC(KG)%VARB_ID(K,L)
                  if (ierr .ne. 0) stop 'variable ID read error in metset'
               endif
               KOL=KOL+8
            ELSE
               READ(HEADER(KOL:KOL+9),'(A4,I5)',iostat=ierr) DREC(KG)%VARB_ID(K,L), DREC(KG)%CHK_SUM(K,L)
               if (ierr .ne. 0) then
                  DREC(KG)%CHK_SUM(K,L) = 0
                  READ(HEADER(KOL:KOL+9),'(A4,I5)',iostat=ierr) DREC(KG)%VARB_ID(K,L)
                  if (ierr .ne. 0) stop 'variable ID read error in metset'
               endif
               KOL=KOL+10
            END IF
            NREC=NREC+1
         END DO
      END DO
      DREC(KG)%REC_PER=NREC
      DREC(KG)%OFFSET=0

!     set the year's accumulated clock time in minutes
!     for the first time period on the file
      CALL TM2MIN(FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,           &
     &            FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR,           &
     &            FILE(KG,KT)%FIRST%MN, FILE(KG,KT)%FIRST%MACC)

!     skip to the next time period index record to find the time
!     interval between data periods (minutes)
      NREC=NREC+1

      READ(KUNIT,REC=NREC,ERR=900)LABEL,HEADER(1:108)
      READ(LABEL,'(5I2,4X,A4)')                                         &
     &   FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                      &
     &   FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                      &
     &   FILE(KG,KT)%LAST%IC, KVAR
      IF(KVAR.NE.'INDX')GO TO 950

!     decode the minutes from the extended header
      READ(HEADER(1:108),'(7X,I2)') FILE(KG,KT)%LAST%MN

!     set the year's accumulated clock time in minutes
!     for the next time period on the file
      CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,             &
     &            FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,             &
     &            FILE(KG,KT)%LAST%MN, FILE(KG,KT)%LAST%MACC)


!     compute the minute difference between time periods
      DREC(KG)%DELTA=FILE(KG,KT)%LAST%MACC-FILE(KG,KT)%FIRST%MACC

!     continue on to find the last record
      NREC=NREC+DREC(KG)%REC_PER
      DO WHILE (NREC.GE.1)

         READ(KUNIT,REC=NREC,ERR=500) LABEL, HEADER(1:108)
         READ(LABEL,'(5I2,4X,A4)',ERR=500)                              &
     &      FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                   &
     &      FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                   &
     &      FILE(KG,KT)%LAST%IC, KVAR
         IF(KVAR.NE.'INDX')GO TO 500
         READ(HEADER(1:108),'(7X,I2)')FILE(KG,KT)%LAST%MN
         NREC=NREC+DREC(KG)%REC_PER

      END DO
  500 FILE(KG,KT)%ENDREC=NREC-1

!     set the year's accumulated clock time in minutes
!     for the last time period on the file
      CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,             &
     &            FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,             &
     &            FILE(KG,KT)%LAST%MN, FILE(KG,KT)%LAST%MACC)

      RETURN

!     abnormal terminations
  900 WRITE(*,*)'WARNING metset: Only one time period of meteo data'
      FILE(KG,KT)%LAST=FILE(KG,KT)%FIRST
      FILE(KG,KT)%ENDREC=NREC-1
!     set to dummy value to avoid division by zero
      DREC(KG)%DELTA=1
      RETURN

  950 WRITE(*,*)'ERROR metset: 2nd time period INDX record missing'
      RETURN

      END
