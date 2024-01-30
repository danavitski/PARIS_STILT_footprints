!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSSET           EMiSsion SET data entry
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION TERM SET IS THE DATA ENTRY FOR EMISSION RATES, POLLUTANT
!   INFORMATION, STARTING TIME AND RELEASE DURATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 01 Apr 1997 - RRD
!
! USAGE:  CALL EMSSET(NUMTYP,IUNIT,IBYR,IBMO,IBDA,IBHR)
!   INPUT ARGUMENT LIST:
!     NUMTYP      - int number of pollutant types
!     IUNIT - int unit number for input data
!     IBYR...   - int   starting time
!   OUTPUT ARGUMENT LIST:
!     COMMON GBLCON
!   INPUT FILES:
!     UNIT 5 or UNIT 30 if input from file CONTROL
!   OUTPUT FILES:
!     UNIT 31 if input from unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: emsset.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE EMSSET(NUMTYP,IUNIT,IBYR,IBMO,IBDA,IBHR)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

!      COMMON /GBLCON/ CONC, DIRT

!==>generic defaults

      NUMTYP=1
      DIRT(1)%IDENT='TEST'
      DIRT(1)%QRATE=1.0
      DIRT(1)%QHRS=1.0
      DIRT(1)%START%YR=IBYR
      DIRT(1)%START%MO=IBMO
      DIRT(1)%START%DA=IBDA
      DIRT(1)%START%HR=IBHR
      DIRT(1)%START%MN=0

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Number of different pollutants'
         WRITE(*,*)NUMTYP
      END IF
      READ(IUNIT,*)NUMTYP
      IF(IUNIT.EQ.5)WRITE(31,*)NUMTYP

!     test limits
      IF(NUMTYP.GT.MAXTYP)THEN
         WRITE(*,*)'ERROR emsset: Number of pollutants exceed limit'
         WRITE(*,*)'numtyp: ',NUMTYP
         WRITE(*,*)'maxtyp: ',MAXTYP
         STOP
      END IF

      DO KK=1,NUMTYP

!        multiple pollutants copy over default values
         IF(KK.GT.1)THEN
            DIRT(KK)%IDENT=DIRT(KK-1)%IDENT
            DIRT(KK)%QRATE=DIRT(KK-1)%QRATE
            DIRT(KK)%QHRS=DIRT(KK-1)%QHRS
            DIRT(KK)%START%YR=DIRT(KK-1)%START%YR
            DIRT(KK)%START%MO=DIRT(KK-1)%START%MO
            DIRT(KK)%START%DA=DIRT(KK-1)%START%DA
            DIRT(KK)%START%HR=DIRT(KK-1)%START%HR
            DIRT(KK)%START%MN=DIRT(KK-1)%START%MN
         END IF

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Pollutant 4-Character Identification'
            WRITE(*,'(A4)')DIRT(KK)%IDENT
         END IF
         READ(IUNIT,'(A4)')DIRT(KK)%IDENT
         IF(IUNIT.EQ.5)WRITE(31,'(A)')DIRT(KK)%IDENT

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Emission rate (per hour)'
            WRITE(*,*)DIRT(KK)%QRATE
         END IF
         READ(IUNIT,*)DIRT(KK)%QRATE
         IF(IUNIT.EQ.5)WRITE(31,*)DIRT(KK)%QRATE

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Hours of emission'
            WRITE(*,*)DIRT(KK)%QHRS
         END IF
         READ(IUNIT,*)DIRT(KK)%QHRS
         IF(IUNIT.EQ.5)WRITE(31,*)DIRT(KK)%QHRS

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Release start time: year month day hour minute'
            WRITE(*,*)DIRT(KK)%START%YR, DIRT(KK)%START%MO,             &
     &         DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN
            READ(*,*) DIRT(KK)%START%YR, DIRT(KK)%START%MO,             &
     &         DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN
         ELSE
            READ(IUNIT,*)IYR,IMO,IDA,IHR,IMN
            IF(IYR+IMO+IDA+IHR+IMN.NE.0)THEN
               DIRT(KK)%START%YR=IYR
               DIRT(KK)%START%MO=IMO
               DIRT(KK)%START%DA=IDA
               DIRT(KK)%START%HR=IHR
               DIRT(KK)%START%MN=IMN
            END IF
         END IF

         IF(IUNIT.EQ.5)WRITE(31,*)                                      &
     &      DIRT(KK)%START%YR, DIRT(KK)%START%MO,                       &
     &      DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN

         CALL TM2MIN(DIRT(KK)%START%YR, DIRT(KK)%START%MO,              &
     &      DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN,    &
     &      DIRT(KK)%START%MACC)

      END DO

      RETURN
      END
