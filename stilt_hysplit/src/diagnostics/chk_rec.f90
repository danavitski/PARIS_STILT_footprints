PROGRAM CHK_REC

!-------------------------------------------------------------------------------
! simple program to dump the 50 byte ascii header on each record
!-------------------------------------------------------------------------------
! LAST REVISED: 19 Jan 1999 (RRD)
!               15 Dec 2000 (RRD) - fortran90 upgrade
!               18 Oct 2001 (RRD) - extended grid domain
!-------------------------------------------------------------------------------

!  CHARACTER KVAR*4, LABEL*50, HEADER*3072, FDIR*80, FILE*80, MODEL*4
  CHARACTER KVAR*4, LABEL*50, HEADER*8000, FDIR*80, FILE*80, MODEL*4

  LOGICAL :: FTEST, SXTNBIT

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
  OPEN(10,FILE=FDIR(1:KLEN)//FILE,STATUS='OLD',RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

  WRITE(*,*)'directory/file : ',FDIR(1:KLEN)//FILE

! decode the standard portion of the index record
  READ(10,REC=1)LABEL
  READ(10,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,4X,A4)')IYR,IMO,IDA,IHR,IFC,KVAR
  WRITE(*,*)'Opened file: ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'WARNING Old format meteo data grid '
     WRITE(*,*)LABEL
     WRITE(*,*)HEADER(1:108)
     STOP
  END IF

! decode extended portion of the header
  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)',ERR=900)                &
     MODEL,    ICX,       MN,                                              &
     POLE_LAT, POLE_LON,  REF_LAT,                                         &
     REF_LON,  SIZE,      ORIENT,                                          &
     TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
     SYNC_LAT, SYNC_LON,  DUMMY,                                           &
     NX,       NY,        NZ,                                              &
     K_FLAG,   LENH

! close file and reopen with proper length
  CLOSE (10)

  SXTNBIT = MODEL == 'RAMS' .OR. MODEL(1:3) == 'ECX'
  LDAT = NX*NY
  IF (SXTNBIT) THEN
     LDAT=LDAT*2
     WRITE (*,*) '*** 16 bit file ***'
  END IF
  LEN = LDAT+50
  OPEN(10,FILE=FDIR(1:KLEN)//FILE,RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')

! print file diagnostic
  WRITE (*,*) 'Grid size and lrec: ', NX, NY, LDAT, LEN
  WRITE (*,*) 'Header record size: ', LENH

! read entire file and print headers
  KREC=1
  DO
     READ(10,REC=KREC,ERR=800)LABEL,HEADER(1:20)
     READ(LABEL,'(5I2,4X,A4)',ERR=900)IYR,IMO,IDA,IHR,IFC,KVAR
     IF(MOD(KREC,20).EQ.0)READ(*,*)
     IF(KVAR.EQ.'INDX') WRITE(*,'(A)')HEADER(1:80)
     WRITE(*,'(A)')LABEL
     KREC=KREC+1
  END DO

800 STOP

900 WRITE(*,*)'ERROR: decoding header'
    WRITE(*,*)LABEL
    WRITE(*,*)HEADER(1:108)

END PROGRAM chk_rec
