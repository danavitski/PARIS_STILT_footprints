!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFDEL           PUFf DELete eliminate unused puffs/particles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF DELETE ELIMINATE UNUSED PUFFS/PARTICLES FROM ARRAY BASED
!   UPON EITHER OFF-GRID, EXCEEDING AGE, OR BELOW MINIMUM MASS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,
!              HDWP,PAGE,PTYP,PGRD,NSORT,KHMAX,FMASS)
!   INPUT ARGUMENT LIST:
!     KPM   - int  total number of puffs or particles
!     KHMAX - int  maximum hours limit
!     FMASS - real lower mass limit per particle/puff
!   OUTPUT ARGUMENT LIST:
!     MASS  - real mass of pollutant (arbitrary units)
!     XPOS,YPOS   - real puff center positions (grid units)
!     ZPOS  - real puff center height (sigma)
!     SIGH,SIGV   - real horiz (meters) and vert sigma (sigma)
!     HDWP  - int  Horizontal distribution within pollutant
!     PAGE  - int  pollutant age since release (min)
!     PTYP  - int  pollutant type index number
!     PGRD  - int  meteorological grid of puff position
!     NSORT - int  sort index array by position
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: pufdel.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,              &
     &   HDWP,PAGE,PTYP,PGRD,NSORT,KHMAX,FMASS)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes

!     particle positions
      REAL*8   XPOS(MAXPAR),YPOS(MAXPAR),ZPOS(MAXPAR)
!     horizontal and vertical distributions within puff
      REAL*8   SIGH(MAXPAR),SIGV(MAXPAR)
!     puff age, distribution type, pollutant index, met grid
      INTEGER   PAGE(MAXPAR),HDWP(MAXPAR),PTYP(MAXPAR),PGRD(MAXPAR)
!     pollutant mass array (mass species, number of particles/puffs)
      REAL*8   MASS(MAXDIM,MAXPAR)
!     sort index
      INTEGER   NSORT(MAXPAR)

!     need at least one element
      IF(KPM.LT.1)RETURN

!     pointer to existing puff index for testing
      KP=1
!     pointer to new index
      KN=0

      DO WHILE (KP.LE.KPM)
!        check age against limits
         KHRS=PAGE(KP)/60

!        check mass total
         TMASS=MASS(1,KP)
         MM=MAXDIM
         DO WHILE(MM.GT.1)
            TMASS=TMASS+MASS(MM,KP)
            MM=MM-1
         END DO

         IF(PGRD(KP).EQ.0.OR.KHRS.GE.KHMAX.OR.TMASS.LE.FMASS)THEN
!           puff to delete
            KP=KP+1
         ELSE
!           puff to keep
            KN=KN+1
            IF(KP.NE.KN)THEN
               XPOS(KN)=XPOS(KP)
               YPOS(KN)=YPOS(KP)
               ZPOS(KN)=ZPOS(KP)
               PAGE(KN)=PAGE(KP)
               SIGH(KN)=SIGH(KP)
               SIGV(KN)=SIGV(KP)
               HDWP(KN)=HDWP(KP)
               PTYP(KN)=PTYP(KP)
               PGRD(KN)=PGRD(KP)
               MASS(1,KN)=MASS(1,KP)
               MM=MAXDIM
               DO WHILE(MM.GT.1)
                  MASS(MM,KN)=MASS(MM,KP)
                  MM=MM-1
               END DO
            END IF
!           save index position in dummy array
            NSORT(KN)=KN

            KP=KP+1
         END IF
      END DO
!     reset counter to new total
      KPM=KN

      RETURN
      END
