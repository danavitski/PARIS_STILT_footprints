!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKSET           PACking SETup to initialize routines
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACKING SETUP INITIALIZE PACKING ROUTINES PAKREC AND PAKNDX
!   ONE ENTRY PER DATA FILE TO PACK.  MULTIPLE FILES CAN BE PACKED A
!   SAME TIME (UP TO MGRD - SEE STRUCTURE). A UNIQUE UNIT NUMBER MUST
!   BE DEFINED FOR EACH OUTPUT FILE.  THE FILE NAME CONTAINS THE METEO
!   STRUCTURE INFORMATION. AFTER THIS ROUTINE THE OUTPUT FILE SHOULD
!   OPENED TO THE SPECIFIED UNIT.
!
! PROGRAM HISTORY LOG:
!   Last Revision: 14 Feb 1997 (RRD)
!                  13 Jul 1999 (RRD)
!
! USAGE:  CALL PAKSET(LUNIT,FNAME,KREC,NXP,NYP,NZP)
!   INPUT ARGUMENT LIST:
!     LUNIT - int output unit number
!     FNAME - char      file name of METDATA.CFG
!   OUTPUT ARGUMENT LIST:
!     KREC  - int position of index record at time-1
!     NXP,NYP     - int horizontal grid dimensions
!     NZP   - int vertical grid dimension (includes surface)
!   INPUT FILES:
!     UNIT as LUNIT file named METDATA.CFG unless otherwise defined
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE PAKSET(LUNIT,FNAME,KREC,NXP,NYP,NZP)

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'DEFPACK.INC'

!     configuration file name
      CHARACTER FNAME*(*)

!     file existence test
      LOGICAL FTEST

!     save structure between routines
      COMMON / PAKCOM / GV, KG

!==>save grid number counter

      SAVE NG
      DATA NG /0/

      NG=NG+1
      IF(NG.GT.MGRD)THEN
         WRITE(*,*)'ERROR pakset: '
         WRITE(*,*)'Trying to intialize grid numbr: ', NG
         WRITE(*,*)'Exceeding max numbers of grids: ', MGRD
         STOP
      END IF
      KG=NG

!==>check for output format configuration file

      INQUIRE(FILE=FNAME,EXIST=FTEST)
      IF(.NOT.FTEST)THEN
         FNAME='METDATA.CFG'
         INQUIRE(FILE=FNAME,EXIST=FTEST)
         IF(.NOT.FTEST)THEN
            WRITE(*,*)'Unable to find data configuration file: ',FNAME
            WRITE(*,*)'Re-enter the file name ...'
            READ(*,'(A)') FNAME
         END IF
      END IF

      OPEN(LUNIT,FILE=FNAME)
!     character field meteo model identification
      READ(LUNIT,'(20X,A4)') GV(NG)%MODEL

!     integer grid number and vertical coordinate system
!     ksys - 1:sigma  2: pressure  3: z-terrain  4: hybrid
      READ(LUNIT,'(20X,I4)') GV(NG)%IG, GV(NG)%KSYS

!     read 12 grid mapping variables
      READ(LUNIT,'(20X,F10.0)')(GV(NG)%GRIDS(I),I=1,12)

!     grid size in x,y,z
      READ(LUNIT,'(20X,I4)') GV(NG)%NXG, GV(NG)%NYG, GV(NG)%NLVL

!     level heights, variables per level, variable id
      NLVL=GV(NG)%NLVL
      DO NL=1,NLVL
         READ(LUNIT,'(20X,F6.0,I3,99(1X,A4))')                          &
     &   GV(NG)%HEIGHT(NL),NVAR,(GV(NG)%VARB(NV,NL),NV=1,NVAR)
         GV(NG)%NVAR(NL)=NVAR
      END DO
      CLOSE (LUNIT)

!==>compute extended header length

      GV(NG)%LENH=108
      NLVL=GV(NG)%NLVL
      DO L=1,NLVL
         GV(NG)%LENH=GV(NG)%LENH+8
         NVAR=GV(NG)%NVAR(L)
         DO K=1,NVAR
            GV(NG)%LENH=GV(NG)%LENH+8
         END DO
      END DO

!==>check for multiple extended header records

      GV(NG)%LREC=(GV(NG)%NXG*GV(NG)%NYG)
!     check limits
      IF((GV(NG)%LENH).GT.MLEN.OR.(GV(NG)%LREC).LT.108)THEN
         WRITE(*,*)'ERROR pakset: Extended header exceeds limits...'
         WRITE(*,*)'Maximum header length : ',MLEN
         WRITE(*,*)'Required HEADER length: ',GV(NG)%LENH
         WRITE(*,*)'Available (nx*xy)     : ',GV(NG)%LREC
         STOP
      END IF

!     number of header records
      GV(NG)%NHREC=GV(NG)%LENH/GV(NG)%LREC+1
!     bytes in last header record
      GV(NG)%NHBYT=GV(NG)%LENH-(GV(NG)%NHREC-1)*GV(NG)%LREC

!==>compute record count offset from index record for each level

      NTOT=GV(NG)%NHREC
      GV(NG)%NREC(1)=NTOT
      NLVL=GV(NG)%NLVL
      DO K=2,NLVL
         NTOT=NTOT+GV(NG)%NVAR(K-1)
         GV(NG)%NREC(K)=NTOT
      END DO

!     records per time period
      NTOT=NTOT+GV(NG)%NVAR(NLVL)
      GV(NG)%NRPT=NTOT

!     check validity of index record
      KREC=MAX0(KREC,1)
      IF(MOD(KREC-1,NTOT).NE.0)THEN
         WRITE(*,*)'ERROR pakset: position record not even multiple'
         WRITE(*,*)'  Record numb in argument:',KREC
         WRITE(*,*)'  Records per time period:',NTOT
         STOP
      ELSE
         WRITE(*,*)'NOTICE pakset'
         WRITE(*,*)'  Number of index records = ',GV(NG)%NHREC
         WRITE(*,*)'  Number of records /time = ',NTOT
      END IF

!==>set remaining variables

!     internal record counter points to index record
      GV(NG)%MREC=KREC-NTOT

!     set flag for first time group for initialization
      GV(NG)%NEWT=.TRUE.

!     unit number
      GV(NG)%KUNIT=LUNIT

!==>return sizes for variable dimensions

      NXP=GV(NG)%NXG
      NYP=GV(NG)%NYG
      NZP=GV(NG)%NLVL

      RETURN
      END
