      SUBROUTINE SOUND (VDATA,SFCP,SFCT,UTW,VTW)

! JCL:'SOUND' stores the results in a variable called 'LABEL', then
!      writes out to 'profile.txt'

      USE module_defsize
      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)


      PARAMETER (NTAB=12)
! JCL:
      CHARACTER VARB(0:NTAB)*4, UNIT(0:NTAB)*4, LABEL*500
      INTEGER NT(MVAR)
      REAL*8 VDATA(MVAR,MLVL),FACT(0:NTAB),                             &
     &           UTW(MLVL), VTW(MLVL)

!==>table of predefined variables, corresponding label, and units
!   conversion factors to go from file units to printer units

      DATA VARB/'    ','PRSS','TPPA','TPPT','TPP6',                     &
                'UWND','VWND','WWND','SPHU','TEMP','RELH','CPRC','PPRC'/

      DATA UNIT/'    ',' hPa','  mm','  mm','  mm',                     &
                ' m/s',' m/s','mb/h','g/kg','  oC','   .','mm/s','mm/s'/

      DATA FACT/1.0   ,1.0   ,1000. ,1000. ,1000. ,                     &
                1.0   ,1.0   ,3600. ,1000. ,1.0   ,1.0   ,   1e3,   1.0/

!--------------------------------------------------------------------
!     data levels including surface
      NLVL=GRID(1)%NZ

      DO LL=1,NLVL
         NVAR=DREC(1)%NUM_VARB(LL)

!==>convert height units according to vertical coordinate system

!        default vertical motion units in mb/s
         DO NN=1,NTAB
            IF ('WWND' == VARB(NN))UNIT(NN)='mb/h'
         END DO

         IF (DREC(1)%Z_FLAG == 1)THEN
!           pressure sigma levels
            PLEVEL=SFCP*DREC(1)%HEIGHT(LL)

         ELSEIF(DREC(1)%Z_FLAG == 2)THEN
!           absolute pressure units
            PLEVEL=DREC(1)%HEIGHT(LL)
            IF (LL == 1)PLEVEL=SFCP

         ELSEIF(DREC(1)%Z_FLAG == 3)THEN
!           terrain following height units
            ZTOP=20000.0
            IF (DREC(1)%HEIGHT(NLVL).GT.20000.0)ZTOP=34800.0
            FACTOR=1.0-SFCT/ZTOP
            PLEVEL=FACTOR*DREC(1)%HEIGHT(LL)

!           terrain following Z system units in m/s
            DO NN=1,NTAB
               IF ('WWND' == VARB(NN))UNIT(NN)=' m/h'
            END DO

         ELSEIF(DREC(1)%Z_FLAG == 4)THEN
!           ecmwf hybrid coordinate system
            OFFSET=INT(DREC(1)%HEIGHT(LL))
            PSIGMA=DREC(1)%HEIGHT(LL)-OFFSET
            PLEVEL=SFCP*PSIGMA+OFFSET
            IF (LL == 1)PLEVEL=SFCP
         END IF

!        by default assume level = pressure unless PRES variable appears
!        (i.e. terrain data (type=3) will have local pressure variable
! JCL:   not use integer, to keep more significant digits
         LEVEL=NINT(PLEVEL)

!==>match variables defined in file's index record with those variables
!   that have been defined in this subroutine and create a variable number
!   for simple table lookup

         DO KK=1,NVAR
            NT(KK)=0
            DO NN=1,NTAB
               IF (DREC(1)%VARB_ID(KK,LL) == VARB(NN))NT(KK)=NN
            END DO
         END DO

!==>write out variable and label information before records 1 and 2
!   which correspond to the surface and first vertical level

         IF (LL == 1)THEN
!           write out surface variables and units
! CHG(3/21/02): don't write too much to std output
            WRITE(*,*)' '
!            WRITE(30,'(A)')' '
! JCL:
!            WRITE(30,'(10X,12(6X,A4))')
!           WRITE(30,'(8X,12(4X,A4))')
!     &           (DREC(1)%VARB_ID(KK,LL),KK=1,NVAR)
! JCL:
!            WRITE(30,'(10X,12(6X,A4))')
!           WRITE(30,'(8X,12(4X,A4))')
!     &           (UNIT(NT(KK)),KK=1,NVAR)

         ELSEIF(LL == 2)THEN
!           write out upper level variables and units
! CHG(3/21/02): don't write too much to std output
!            WRITE(*,*)' '
!            WRITE(30,'(A)')' '

! JCL:
!            WRITE(LABEL,'(10X,12(6X,A4))')                              &
!     &           (DREC(1)%VARB_ID(KK,LL),KK=1,NVAR)
!            WRITE(LABEL(81:),'(3A10)')'      TPOT','      UWND',       &
!     &                                '      VWND'
!            WRITE(30,'(A)')LABEL(1:110)

!            WRITE(LABEL,'(10X,12(6X,A4))')(UNIT(NT(KK)),KK=1,NVAR)
!            WRITE(LABEL(81:),'(3A10)')'        oK','      W->E',        &
!     &                                '      S->N'
!            WRITE(30,'(A)')LABEL(1:110)
         END IF

!==>convert each variable at that level to standard units as defined
!   from the table lookup. Variables not found are not converted and
!   have no specific units label.

         DO KK=1,NVAR
            VDATA(KK,LL)=VDATA(KK,LL)*FACT(NT(KK))
         END DO


! JCL:*****************
!==>write results using a variable precision format depending upon value
! JCL:have separate formats for the first level vs the upper levels
      IF (LL == 1)THEN
         WRITE(LABEL,'(I10)')LEVEL
         DO KK=1,NVAR

!           real and integer equivalents
            RNUM=VDATA(KK,LL)
            KNUM=NINT(RNUM)

!           print column for variable
            K1=KK*10+1
            K2=KK*10+10

!           select format according to value
            IF (ABS(RNUM).GE.1000.0)THEN
               WRITE(LABEL(K1:K2),'(F10.3)')RNUM
            ELSEIF(RNUM.GE.100.0)THEN
               WRITE(LABEL(K1:K2),'(F10.4)')RNUM
            ELSEIF(ABS(RNUM).GE.1.0)THEN
               WRITE(LABEL(K1:K2),'(F10.4)')RNUM
            ELSEIF(ABS(RNUM).GE.0.1)THEN
               WRITE(LABEL(K1:K2),'(F10.5)')RNUM
            ELSEIF(RNUM == 0.0)THEN
               WRITE(LABEL(K1:K2),'(I10)')KNUM
            ELSE
               NEXP=0
               DO WHILE (KNUM.LE.0 .AND. NEXP.LT.9)
                  NEXP=NEXP+1
                  RNUM=RNUM*10.0
                  KNUM=NINT(ABS(RNUM))
               END DO
               IF (NEXP == 9)KNUM=0
               KNUM=NINT(RNUM)
               WRITE(LABEL(K1:K2),'(I7,A2,I1)')KNUM,'E-',NEXP
! JCL:         WRITE(LABEL(K1:K2),'(I3,A2,I1)')KNUM,'E-',NEXP
            END IF
         END DO
!==>write results using a variable precision format depending upon value
! JCL:***
      ELSE
! JCL:   write out 'PLEVEL' instead of 'LEVEL', b/c has more sign digits
         WRITE(LABEL,'(F10.4)')PLEVEL
         DO KK=1,NVAR

!           real and integer equivalents
            RNUM=VDATA(KK,LL)
            KNUM=NINT(RNUM)

!           print column for variable
            K1=KK*10+1
            K2=KK*10+10

!           select format according to value
            IF (ABS(RNUM).GE.100.0)THEN
               WRITE(LABEL(K1:K2),'(I10)')KNUM
            ELSEIF(ABS(RNUM).GE.1.0)THEN
               WRITE(LABEL(K1:K2),'(F10.5)')RNUM
            ELSEIF(ABS(RNUM).GE.0.1)THEN
               WRITE(LABEL(K1:K2),'(F10.5)')RNUM
            ELSEIF(RNUM == 0.0)THEN
               WRITE(LABEL(K1:K2),'(I10)')KNUM
            ELSE
               NEXP=0
               DO WHILE (KNUM.LE.0 .AND. NEXP.LT.9)
                  NEXP=NEXP+1
                  RNUM=RNUM*10.0
                  KNUM=NINT(ABS(RNUM))
               END DO
               IF (NEXP == 9)KNUM=0
               KNUM=NINT(RNUM)
! JCL          WRITE(LABEL(K1:K2),'(I3,A2,I1)')KNUM,'E-',NEXP
               WRITE(LABEL(K1:K2),'(I7,A2,I1)')KNUM,'E-',NEXP
            END IF

!           save temperature and pressure for potential temp
            IF (DREC(1)%VARB_ID(KK,LL) == 'TEMP')                        &
     &         TEMP=VDATA(KK,LL)+273.16
            IF (DREC(1)%VARB_ID(KK,LL) == 'PRES')                        &
     &         PLEVEL=VDATA(KK,LL)
         END DO
! JCL:***
      END IF
! JCL:*****************

!        add rotated winds to profile above surface values
!         IF(LL.GE.2)THEN
!            TPOT=TEMP*(1000.0/PLEVEL)**0.286
! JCL:      change to output MORE significant figs
!            WRITE(LABEL(81:),'(3F10.5)')TPOT,UTW(LL),VTW(LL)
!            WRITE(LABEL(81:),'(3F6.1)')TPOT,UTW(LL),VTW(LL)
!         END IF

! JCL:
!        print data record
!         DO KK=1,K2,78
!            WRITE(*,'(A)')LABEL(KK:KK+77)
!            WRITE(30,'(A)')LABEL(KK:KK+77)
!         END DO
!        print data record
         DO KK=1,K2,140
! CHG(3/21/02): don't write too much to std output
!            WRITE(*,'(A)')LABEL(KK:KK+139)
!            WRITE(30,'(A)')LABEL(KK:KK+139)
         END DO

      END DO

      END SUBROUTINE SOUND
