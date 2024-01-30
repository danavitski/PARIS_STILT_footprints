!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPRAD           DEPosition RADioactive decay
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION RADIOACTIVE COMPUTES THE RADIOACTIVE DECAY OF
!   POLLUTANTS THAT HAVE ALREADY BEEN DEPOSITED.  THE DECAY IS
!   APPLIED TO THE VALUES THAT ARE BEING SUMMED IN THE DEPOSITION
!   ARRAY EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL DEPRAD(NUMGRD,NUMTYP,DT,CSUM)
!   INPUT ARGUMENT LIST:
!     NUMGRD      - int number of concentration grids
!     NUMTYP      - int number of pollutants
!     DT    - real      intgration time step (min)
!     CSUM  - real      concentration summation matrix
!   OUTPUT ARGUMENT LIST:
!     CSUM  - real      concentration summation matrix
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: deprad.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE DEPRAD(NUMGRD,NUMTYP,DT,CSUM)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

!     master concentration array (x,y,z,grids,species)
      REAL*8   CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!      COMMON /GBLCON/ CONC, DIRT

!     go through each grid
      DO KG=1,NUMGRD

!        determine loop indicies
         NXP=CONC(KG)%NUMB_LON
         NYP=CONC(KG)%NUMB_LAT
         NZP=CONC(KG)%LEVELS

         DO KT=1,NUMTYP
            IF(DIRT(KT)%DORAD)THEN
!              convert half-life to time constant
               RTC=-ALOG(0.5)/(DIRT(KT)%RHALF*86400.0)
!              convert time constant to removal fraction
               RFR=DEXP(-60.0*DT*RTC)

               DO KL=1,NZP
!                 zero height is for deposition
                  IF(CONC(KG)%HEIGHT(KL).EQ.0.0)THEN
                     DO JJ=1,NYP
                     DO II=1,NXP
                        CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)*RFR
                     END DO
                     END DO
                  END IF
               END DO

            END IF
         END DO

      END DO
      RETURN
      END
