!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONDSK           CONcentration output to DiSK
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION TO DISK WRITES OUT THE CONCENTRATION MATRIX TO
!   DISK FILES AT PRESELECTED SAMPLING INTERVALS.  AFTER OUTPUT
!   ELEMENTS ARE SET TO ZERO FOR THE NEXT CYLCE ACCUMULATION.
!   APPROPRIATE TIME INDEX RECORDS ARE WRITTEN TO EACH FILE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL CONDSK(NUMGRD,NUMTYP,JET,IYR,IFHR,CSUM)
!   INPUT ARGUMENT LIST:
!     NUMGRD    - int   number of concentration grids
!     NUMTYP    - int   number of pollutants
!     JET       - int   current elapsed time
!     IYR       - int   current calculation year
!     IFHR      - int   current forecast hour
!     CSUM      - real  concentration summation matrix
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNITS 21,22,etc - as defined in input CONTROL file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: condsk.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

! JCL:  add an extra argument 'BACK' to skip over tests
!       below if BACK==T
      SUBROUTINE CONDSK(NUMGRD,NUMTYP,JET,IFHR,CSUM,BACK)
!      SUBROUTINE CONDSK(NUMGRD,NUMTYP,JET,IFHR,CSUM)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

! JCL:
      LOGICAL BACK

!      COMMON /GBLCON/ CONC, DIRT

!     master concentration array (x,y,z,grids,species)
      REAL*8   CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!     convert current time as end sampling period
      CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)

!     go through each grid
      DO KG=1,NUMGRD

! JCL: differentiate between when BACK is T or F
         IF (.NOT.BACK) THEN
!        current time must be within output window
            IF((JET.LT.CONC(KG)%START%MACC).OR.                         &
     &         (JET.GT.CONC(KG)%STOP%MACC))THEN

              GO TO 100
            END IF

!        current time must be at the output interval
            MTIME=JET-CONC(KG)%START%MACC
            IF(MTIME.LE.0)THEN

              GO TO 100
            END IF
! JCL:      checks to see if current time is multiple of sampling
!             interval; if not, then not write output to disk
            IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0)THEN

              GO TO 100
            END IF
! JCL: add the following lines to provide a different test if
!      BACK is TRUE
          ELSE
            IF((JET.GT.CONC(KG)%START%MACC).OR.                         &
     &         (JET.LT.CONC(KG)%STOP%MACC))THEN

              GO TO 100
            END IF
! JCL:  current time must be at the output interval
!       reverse subtraction direction
            MTIME=CONC(KG)%START%MACC-JET
            IF(MTIME.LE.0)THEN
              GO TO 100
            END IF

! JCL:      checks to see if current time is multiple of sampling
!             interval; if not, then not write output to disk
            IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0)THEN
              GO TO 100
            END IF
! JCL: end of IF/ELSE statement for BACK
         END IF


         KUNIT=CONC(KG)%UNIT
!        write index record for sampling start time
         WRITE(KUNIT,*)                                                 &
     &      CONC(KG)%NOW%YR, CONC(KG)%NOW%MO, CONC(KG)%NOW%DA,          &
     &      CONC(KG)%NOW%HR, CONC(KG)%NOW%MN, CONC(KG)%NOW%IC

!        write index record for sampling stop time
         WRITE(KUNIT,*)IYR,IMO,IDA,IHR,IMN,IFHR

!        determine loop indicies
         NXP=CONC(KG)%NUMB_LON
         NYP=CONC(KG)%NUMB_LAT
         NZP=CONC(KG)%LEVELS


         DO KT=1,NUMTYP
         DO KL=1,NZP
!           write data record by pollutant type and level
            WRITE(KUNIT,*)DIRT(KT)%IDENT, CONC(KG)%HEIGHT(KL),          &
     &         ((CSUM(II,JJ,KL,KT,KG),II=1,NXP),JJ=1,NYP)

         END DO
         END DO

         WRITE(30,*)'NOTICE condsk: output at ',IMO,IDA,IHR

  100 CONTINUE
      END DO

      RETURN
      END
