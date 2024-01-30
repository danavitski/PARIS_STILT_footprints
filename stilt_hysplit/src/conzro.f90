!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONZRO           CONcentration array set to ZeRO
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   MATRIX ELEMENTS ARE SET TO ZERO FOR THE NEXT CYLCE ACCUMULATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL CONZRO(NUMGRD,NUMTYP,DT,JET,IFHR,CSUM)
!   INPUT ARGUMENT LIST:
!     NUMGRD    - int   number of concentration grids
!     NUMTYP    - int   number of pollutants
!     DT        - real  integration time step (min)
!     JET       - int   current elapsed time
!     IFHR      - int   current forecast hour
!     CSUM      - real  concentration summation matrix
!   OUTPUT ARGUMENT LIST:
!     CSUM      - real  concentration summation matrix
!     COMMON GBLCON
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: conzro.f90,v 1.4 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

! JCL:  add BACK as argument
      SUBROUTINE CONZRO(NUMGRD,NUMTYP,DT,JET,IFHR,CSUM,BACK)
!      SUBROUTINE CONZRO(NUMGRD,NUMTYP,DT,JET,IFHR,CSUM)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

!      COMMON /GBLCON/ CONC, DIRT

! JCL:
      LOGICAL BACK

!     master concentration array (x,y,z,grids,species)
      REAL*8   CSUM(MAXXP,MAXYP,MAXZP,MAXTYP,MAXGRD)

!     go through each grid
      DO KG=1,NUMGRD

!        convert current time as end sampling period
         IF(CONC(KG)%SNAP.EQ.0)THEN
            CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)
         ELSE
            CALL TM2DAY(JET+INT(DT),IYR,IMO,IDA,IHR,IMN)
         END IF

! JCL:   add BACK into logic
!        current time must be within output window
         IF((.NOT.BACK).AND.((JET.LT.CONC(KG)%START%MACC).OR.           &
     &      (JET.GT.CONC(KG)%STOP%MACC)))GO TO 100
! JCL:   when BACK is T
         IF((BACK).AND.((JET.GT.CONC(KG)%START%MACC).OR.                &
     &      (JET.LT.CONC(KG)%STOP%MACC)))GO TO 100

!        snapshot maps always zero out array, when output averaged
!        test for end of averaging interval
         IF(CONC(KG)%SNAP.EQ.0)THEN
! JCL:      have two different scenarios depending on BACK
!           current time must be at the output interval
            IF (.NOT.BACK)THEN
               MTIME=JET-CONC(KG)%START%MACC
! JCL:      if BACK is T
            ELSE
               MTIME=CONC(KG)%START%MACC-JET
            END IF
            IF(MTIME.LE.0)GO TO 100
            IF(MOD(MTIME,CONC(KG)%DELTA%MACC).NE.0)GO TO 100
         END IF

!        determine loop indicies
         NXP=CONC(KG)%NUMB_LON
         NYP=CONC(KG)%NUMB_LAT
         NZP=CONC(KG)%LEVELS

         DO KT=1,NUMTYP
         DO KL=1,NZP

!           zero out that level pollutant
            DO JJ=1,NYP
            DO II=1,NXP
               CSUM(II,JJ,KL,KT,KG)=0.0
            END DO
            END DO

         END DO
         END DO

!        save end time as start of next sampling period
         CONC(KG)%NOW%YR=IYR
         CONC(KG)%NOW%MO=IMO
         CONC(KG)%NOW%DA=IDA
         CONC(KG)%NOW%HR=IHR
         CONC(KG)%NOW%MN=IMN
         CONC(KG)%NOW%IC=IFHR

  100 CONTINUE
      END DO

      RETURN
      END
