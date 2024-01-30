!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFMRG           PUFf MeRGes puffs together at same position
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF MERGE MERGES PUFFS TOGETHER WHEN AT SAME POSITION AND
!   OF THE SAME HORIZONATL DISTRIBUTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 - RRD
!
! USAGE:  CALL PUFMRG(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,
!              MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,PGRD,NSORT)
!   INPUT ARGUMENT LIST:
!     FRACH - real horizontal position rounding fraction
!     FRACV - real vertical position rounding fraction
!     FRACT - real travel-time rounding fraction
!     FMASS - real mass cutoff for sortting
!     KPM   - int  total number of puffs or particles
!     NSORT - int  sortted array by position
!   OUTPUT ARGUMENT LIST:
!     ZMDL  - real model top scaling depth
!     MASS  - real mass of pollutant (arbitrary units)
!     XPOS,YPOS   - real puff center positions (grid units)
!     ZPOS  - real puff center height (sigma)
!     SIGH,SIGV   - real horiz (meters) and vert sigma (sigma)
!     HDWP  - int  Horizontal distribution within pollutant
!     PAGE  - int  pollutant age since release (min)
!     PTYP  - int  pollutant type index number
!     PGRD  - int  meteorological grid of puff position
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

      SUBROUTINE PUFMRG(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,               &
     &   MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,PAGE,PTYP,PGRD,NSORT)

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
      INCLUDE 'DEFSIZE.INC'
!     meteorology file and grid
      INCLUDE 'DEFGRID.INC'

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER   PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8   MASS(MAXDIM,MAXPAR), MASS0(MAXDIM)

!     sort index and temporary sort index
      INTEGER   NSORT(MAXPAR)

!     sort flag flag
      LOGICAL SORTF, EQUAL

      COMMON /GBLGRD/ GRID, DREC, FILE

!     need at least two elements to merge
      IF(KPM.LT.2)RETURN

      KP=1
      INDEX1=NSORT(KP)
!     find pointer to first non-zero (non-particle) index
      DO WHILE (HDWP(INDEX1).LT.1.AND.KP.LE.KPM)
         KP=KP+1
         INDEX1=NSORT(KP)
      END DO
!     none to merge in the array (only particles)
      IF(KP.GE.KPM)RETURN

!     point to first element to test
      DO WHILE (KP.LT.KPM)
         INDEX1=NSORT(KP)
         EQUAL=.TRUE.

!        horizontal grid distance (m) for rounding
         HGD=GRID(PGRD(INDEX1))%SIZE*1000.0

         KN=KP
!        loop indicies until true (no sort required)
         DO WHILE (EQUAL.AND.KN.LT.KPM)
            KN=KN+1
            INDEX2=NSORT(KN)
            CALL PUFRND(ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,               &
     &         XPOS(INDEX1),YPOS(INDEX1),ZPOS(INDEX1),SIGH(INDEX1),     &
     &         SIGV(INDEX1),HDWP(INDEX1),PTYP(INDEX1),                  &
     &         MASS(1,INDEX1),TMASS1,                                   &
     &         XPOS(INDEX2),YPOS(INDEX2),ZPOS(INDEX2),SIGH(INDEX2),     &
     &         SIGV(INDEX2),HDWP(INDEX2),PTYP(INDEX2),                  &
     &         MASS(1,INDEX2),TMASS2,SORTF,EQUAL)
         END DO

!        shift pointer back if last required a sort
         IF(.NOT.EQUAL)KN=KN-1

!        low mass merge option skip remaining high mass particles
         IF(TMASS1.GT.0.0.OR.TMASS2.GT.0.0)RETURN

!        more than one observation in merge set
         IF(KN.GT.KP)THEN
!           zero out summation variables
            XPOS0=0.0
            YPOS0=0.0
            ZPOS0=0.0
            PAGE0=0.0
            SIGH0=0.0
            SIGV0=0.0

            MASS0(1)=0.0
            MM=MAXDIM
            DO WHILE(MM.GT.1)
               MASS0(MM)=0.0
               MM=MM-1
            END DO

!           sum each element weighted by total element mass
            DO KK=KP,KN
               INDEX0=NSORT(KK)
               TMASS=MASS(1,INDEX0)
               MASS0(1)=MASS0(1)+MASS(1,INDEX0)
               MM=MAXDIM
               DO WHILE(MM.GT.1)
                  TMASS=TMASS+MASS(MM,INDEX0)
                  MASS0(MM)=MASS0(MM)+MASS(MM,INDEX0)
                  MM=MM-1
               END DO

               XPOS0=XPOS0+XPOS(INDEX0)*TMASS
               YPOS0=YPOS0+YPOS(INDEX0)*TMASS
               ZPOS0=ZPOS0+ZPOS(INDEX0)*TMASS
               PAGE0=PAGE0+PAGE(INDEX0)*TMASS
               SIGH0=SIGH0+SIGH(INDEX0)*TMASS
               SIGV0=SIGV0+SIGV(INDEX0)*TMASS
            END DO

!           compute total mass from elemental sums
!           and replace sums into first element
            INDEX0=NSORT(KP)
            TMASS=MASS0(1)
            MASS(1,INDEX0)=MASS0(1)
            MM=MAXDIM
            DO WHILE(MM.GT.1)
               MASS(MM,INDEX0)=MASS0(MM)
               TMASS=TMASS+MASS0(MM)
               MM=MM-1
            END DO

!           compute mass weighted real variables
!           integer marker variables are unchanged
            XPOS(INDEX0)=XPOS0/TMASS
            YPOS(INDEX0)=YPOS0/TMASS
            ZPOS(INDEX0)=ZPOS0/TMASS
            PAGE(INDEX0)=PAGE0/TMASS
            SIGH(INDEX0)=SIGH0/TMASS
            SIGV(INDEX0)=SIGV0/TMASS

!           grid identification to zero for deleted elements
            DO KK=(KP+1),KN
               INDEX0=NSORT(KK)
               PGRD(INDEX0)=0
            END DO
         END IF
!        start next search cycle
         KP=KN+1

      END DO
      RETURN
      END
