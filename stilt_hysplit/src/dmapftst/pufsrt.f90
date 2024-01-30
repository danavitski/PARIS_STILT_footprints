!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSRT           binary PUFf SoRTting by position
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SORTTING IS THE BINARY SORT OF PUFFS BY ROUNDED POSITION
!   NEEDS TO BE CALLED IMMEDIATELY AFTER PUFDEL
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 - RRD
!
! USAGE:  CALL PUFSRT(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,
!              XPOS,YPOS,ZPOS,MASS,SIGH,SIGV,HDWP,PTYP,PGRD,NSORT)
!   INPUT ARGUMENT LIST:
!     FRACH - real horizontal position rounding fraction
!     FRACV - real vertical position rounding fraction
!     FRACT - real travel-time rounding fraction
!     FMASS - real mass sort cutoff value (0 for none)
!     KPM   - int  total number of puffs or particles
!     ZMDL  - real vertical model scaling height
!     XPOS,YPOS   - real puff center positions (grid units)
!     ZPOS  - real puff center height (sigma)
!     MASS  - real pollutant mass
!     SIGH,SIGV   - real horiz (meters) and vert sigma (sigma)
!     HDWP  - int  Horizontal distribution within pollutant
!     PTYP  - int  pollutant type index number
!     PGRD  - int  meteo grid index for puff
!   OUTPUT ARGUMENT LIST:
!     NSORT - int  sortted array by position
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE PUFSRT(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,               &
     &   XPOS,YPOS,ZPOS,MASS,SIGH,SIGV,HDWP,PTYP,PGRD,NSORT)

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
      INCLUDE 'DEFSIZE.INC'
!     meteorology grid and file
      INCLUDE 'DEFGRID.INC'

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR)
!     distribution type, pollutant index, meteo index
      INTEGER   HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass
      REAL*8   MASS(MAXDIM,MAXPAR)

!     sort index and temporary sort index
      INTEGER   NSORT(MAXPAR), VAR(MAXPAR)
!     sort flag and equality flags
      LOGICAL SORTF, EQUAL

      COMMON /GBLGRD/ GRID, DREC, FILE

!     need at least two elements to sort
      IF(KPM.LT.2)RETURN

!     length of primary sort group
      KP=1

!     there are 2^n passes through n data points
      DO WHILE (KP.LT.KPM)
!        start index of sortted output
         KN=1

!        pass through the data by group
         DO WHILE (KN.LE.KPM)
!           first group starting pointer
            IG1=KN
!           second group starting pointer
            IG2=KN+KP

!           number of elements in group being tested
            LEN1=MAX0(MIN0(KPM-IG1+1,KP),0)
            LEN2=MAX0(MIN0(KPM-IG2+1,KP),0)

!           pointer to elements in each group
            L1=0
            L2=0

!           loop through the elements within a group
            DO WHILE (L1.LT.LEN1.AND.L2.LT.LEN2)
!              convert pointer to array index number
               INDEX1=NSORT(IG1+L1)
               INDEX2=NSORT(IG2+L2)

!              horizontal meteo grid distance (m) for rounding
               HGD=GRID(PGRD(INDEX1))%SIZE*1000.0

!              round position according to sigma and determine sort
               CALL PUFRND(ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,            &
     &            XPOS(INDEX1),YPOS(INDEX1),ZPOS(INDEX1),SIGH(INDEX1),  &
     &            SIGV(INDEX1),HDWP(INDEX1),PTYP(INDEX1),               &
     &            MASS(1,INDEX1),TMASS1,                                &
     &            XPOS(INDEX2),YPOS(INDEX2),ZPOS(INDEX2),SIGH(INDEX2),  &
     &            SIGV(INDEX2),HDWP(INDEX2),PTYP(INDEX2),               &
     &            MASS(1,INDEX2),TMASS2,SORTF,EQUAL)

               IF(SORTF)THEN
                  VAR(KN)=NSORT(IG2+L2)
                  L2=L2+1
                  KN=KN+1
               ELSE
                  VAR(KN)=NSORT(IG1+L1)
                  L1=L1+1
                  KN=KN+1
               END IF
            END DO

!           when one group is finished ... copy the other
            DO WHILE (L1.LT.LEN1.AND.L2.EQ.LEN2)
               VAR(KN)=NSORT(IG1+L1)
               L1=L1+1
               KN=KN+1
            END DO

            DO WHILE (L2.LT.LEN2.AND.L1.EQ.LEN1)
               VAR(KN)=NSORT(IG2+L2)
               L2=L2+1
               KN=KN+1
            END DO
         END DO

!        save sortted index structure for next pass
         DO I=1,KPM
            NSORT(I)=VAR(I)
         END DO

!        double size of sort group
         KP=KP*2
      END DO
      RETURN
      END
