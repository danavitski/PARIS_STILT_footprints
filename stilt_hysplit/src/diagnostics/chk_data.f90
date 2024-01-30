PROGRAM CHK_DATA

!-------------------------------------------------------------------------------
! Simple program to dump the first few elements of the data array for each
! record of an ARL packed meteorological file. Used for diagnostic testing.
! Created: 23 Nov 1999 (RRD)
!          14 Dec 2000 (RRD) - fortran90 upgrade
!          18 Oct 2001 (RRD) - expanded grid domain
!-------------------------------------------------------------------------------

  REAL,          ALLOCATABLE :: RDATA(:,:)   
  CHARACTER(1),  ALLOCATABLE :: CPACK(:)

  CHARACTER(4)               :: KVAR, MODEL 
  CHARACTER(50)              :: LABEL          
  CHARACTER(80)              :: FDIR, FILE
  CHARACTER(8000)            :: HEADER
  LOGICAL                    :: FTEST, SXTNBIT

!-------------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE UNPACK (SXTNBIT,CPACK,RDATA,NX,NY,NEXP,VAR1)
     LOGICAL     ,INTENT(IN)  :: SXTNBIT
     CHARACTER(1),INTENT(IN)  :: CPACK(:)
     REAL,        INTENT(OUT) :: RDATA(:,:)
     INTEGER,     INTENT(IN)  :: NX,NY,NEXP
     REAL,        INTENT(IN)  :: VAR1
     END SUBROUTINE
  END INTERFACE
!-------------------------------------------------------------------------------

! directory and file name
  WRITE(*,*)'Enter directory name:'
  READ(*,'(a)')FDIR
  WRITE(*,*)'Enter file name:'
  READ(*,'(a)')FILE

! test for meteo file existence
  KLEN=INDEX(FDIR,' ')-1
  INQUIRE(FILE=FDIR(1:KLEN)//FILE,EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(*,*)'Unable to find file: ',FILE
     WRITE(*,*)'On local directory : ',FDIR(1:KLEN)
     STOP
  END IF

! open file to decode the standard label (50) plus the
! fixed portion (108) of the extended header
  OPEN(10,FILE=FDIR(1:KLEN)//FILE,RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

! decode the standard portion of the index record
  READ(10,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,4X,A4)')IYR,IMO,IDA,IHR,IFC,KVAR
  WRITE(*,'(A,4I5)')'Opened file       : ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'WARNING Old format meteo data grid'
     WRITE(*,*)LABEL
     WRITE(*,*)HEADER(1:108)
     STOP
  END IF

! decode extended portion of the header
  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)',ERR=900)                    &
         MODEL,    ICX,       MN,                                              &
         POLE_LAT, POLE_LON,  REF_LAT,                                         &
         REF_LON,  SIZE,      ORIENT,                                          &
         TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
         SYNC_LAT, SYNC_LON,  DUMMY,                                           &
         NX,       NY,        NZ,                                              &
         K_FLAG,   LENH

! close file and reopen with proper length
  CLOSE (10)
  NXY = NX*NY
!RAMS need more bits...
  SXTNBIT = MODEL == 'RAMS' .OR. MODEL(1:3) == 'ECX'
  IF (SXTNBIT) THEN
     NXY=NXY*2
     WRITE (*,*)'*** 16 bit file ***'
  END IF
  LEN = NXY+50
  OPEN(10,FILE=FDIR(1:KLEN)//FILE,RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')

! print file diagnostic
  WRITE(*,'(A,4I5)')'Grid size and lrec: ',NX,NY,NXY,LEN
  WRITE(*,'(A,I5)') 'Header record size: ',LENH

! allocate array space
  ALLOCATE (RDATA(NX,NY), STAT=KRET)   
  ALLOCATE (CPACK(NXY),   STAT=KRET)

! read entire file and print headers
  KREC=1
100 READ(10,REC=KREC,ERR=800)LABEL,(CPACK(K),K=1,NXY)
    READ(LABEL,'(6I2,2X,A4,I4,2E14.7)',ERR=900) IY,IM,ID,IH,IF,KL,  &
                                             KVAR,NEXP,PREC,VAR1

    WRITE(*,'(A)')LABEL
    IF(KVAR.NE.'INDX') CALL UNPACK (SXTNBIT,CPACK,RDATA,NX,NY,NEXP,VAR1)

    READ(*,*,END=800)
    KREC=KREC+1
  GO TO 100

800 STOP

900 WRITE(*,*)'ERROR: decoding header'
    WRITE(*,*)LABEL
    WRITE(*,*)HEADER(1:108)

END PROGRAM chk_data

!-------------------------------------------------------------------------------

SUBROUTINE UNPACK (SXTNBIT,CPACK,RDATA,NX,NY,NEXP,VAR1)

  IMPLICIT NONE

  LOGICAL     ,INTENT(IN)  :: SXTNBIT
  CHARACTER(1),INTENT(IN)  :: CPACK(:)  
  REAL,        INTENT(OUT) :: RDATA(:,:)   
  INTEGER,     INTENT(IN)  :: NX,NY,NEXP
  REAL,        INTENT(IN)  :: VAR1

  INTEGER :: INDX, I, J
  REAL    :: SCALE, VOLD

! only required when dealing with SUN F90 compiler
! replace ICHAR below with internally defined JCHAR function
! CHARACTER MYCHR*1
! JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

!-------------------------------------------------------------------------------
  IF (SXTNBIT) STOP 'UNPACK: 16 bit files not implemented. Stop.'
  SCALE=2.0**(7-NEXP)
  VOLD=VAR1
  INDX=0
  DO J=1,NY
     DO I=1,NX
        INDX=INDX+1
        RDATA(I,J)=(ICHAR(CPACK(INDX))-127.)/SCALE+VOLD
        VOLD=RDATA(I,J)
        IF(I.LE.2.AND.J.LE.2)  &
           WRITE(*,'(3I5,g12.4)')J,I,ICHAR(CPACK(INDX)),RDATA(I,J)
        IF(I.GE.(NX-1).AND.J.GE.(NY-1))  &
           WRITE(*,'(3I5,g12.4)')J,I,ICHAR(CPACK(INDX)),RDATA(I,J)
     END DO
     VOLD=RDATA(1,J)
  END DO

END SUBROUTINE unpack
