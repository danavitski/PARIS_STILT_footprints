!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPSET           DEPosition parameters SET from input data
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION PARAMETERS SET IS THE DATA ENTRY FOR POLLUTANT SPECIES
!   REQUIRED FOR GRAVITATIONAL SETTLING, DRY DEPOSITION, WET REMOVAL
!   AND RADIOACTIVE DECAY COMPUTATIONS. ONE SET OF ENTIRES REQUIRED
!   EACH DEFINED POLLUTANT TYPE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL DEPSET(NUMTYP,IUNIT,CDEP,RDEP,SDEP)
!   INPUT ARGUMENT LIST:
!     NUMTYP      - int number of pollutant types
!     IUNIT - int unit number for input data
!   OUTPUT ARGUMENT LIST:
!     CDEP  - log flag to indicate wet or dry deposition
!     RDEP  - log flag to indicate resistance deposition
!     SDEP  - log flag to indicate resuspension option
!   INPUT FILES:
!     UNIT 5 or UNIT 30 depending if CONTROL file is found
!   OUTPUT FILES:
!     UNIT 31 to file STARTUP if input from unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: depset.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE DEPSET(NUMTYP,IUNIT,CDEP,RDEP,SDEP)

      use module_defconc
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     pollutant and concentration grid
!      INCLUDE 'DEFCONC.INC'

      LOGICAL CDEP, RDEP, SDEP

!      COMMON /GBLCON/ CONC, DIRT

      CDEP=.FALSE.
      RDEP=.FALSE.
      SDEP=.FALSE.

!=>number of pollutants already defined

      IF(IUNIT.EQ.5)WRITE(*,*)'Enter number of pollutants:', NUMTYP
      READ(IUNIT,*)NUMPOL
      IF(NUMPOL.NE.NUMTYP)THEN
         WRITE(*,*)'WARNING depset: # pollutants defined -',NUMTYP
         NUMPOL=NUMTYP
      END IF
      IF(IUNIT.EQ.5)WRITE(31,*)NUMPOL

      DO 100 KK=1,NUMTYP
         IF(IUNIT.EQ.5)                                                 &
     &      WRITE(*,*)'Enter data for Pollutant:',DIRT(KK)%IDENT

!=>define as gas or particle

         DIRT(KK)%DOGAS=.FALSE.
         DIRT(KK)%PDIAM=0.0
         DIRT(KK)%PDENS=0.0
         DIRT(KK)%SHAPE=0.0

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Define as gas (all zero) or particle'
            WRITE(*,*)'Diameter (um), density (g/cc), shape (1-2)'
            WRITE(*,*)DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE
         END IF

         READ(IUNIT,*)DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE
         IF(IUNIT.EQ.5)                                                 &
     &     WRITE(31,*)DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE

!        all three options required to define a particle
         IF(DIRT(KK)%PDIAM.EQ.0.0 .OR. DIRT(KK)%PDENS.EQ.0.0 .OR.       &
     &      DIRT(KK)%SHAPE.EQ.0.0)     DIRT(KK)%DOGAS=.TRUE.

!=>define dry deposition velocity (over-rides grav settling)

         DIRT(KK)%DORES=.FALSE.
         DIRT(KK)%DODRY=.FALSE.
         DIRT(KK)%DRYVL=0.0
         DIRT(KK)%GPMOL=0.0
         DIRT(KK)%ACVTY=0.0
         DIRT(KK)%DIFTY=0.0
         DIRT(KK)%HENRY=0.0

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Dry deposition - '
            WRITE(*,*)'---- explicit ----  --- resistence method --->'
            WRITE(*,*)'Set velocity (m/s), Gram Molecular Wt',          &
     &            'Activity Ratio, Diffusivity Ratio, Henrys Constant'
            WRITE(*,*)DIRT(KK)%DRYVL, DIRT(KK)%GPMOL, DIRT(KK)%ACVTY,   &
     &                DIRT(KK)%DIFTY, DIRT(KK)%HENRY
         END IF

         READ(IUNIT,*)DIRT(KK)%DRYVL, DIRT(KK)%GPMOL, DIRT(KK)%ACVTY,   &
     &                DIRT(KK)%DIFTY, DIRT(KK)%HENRY
         IF(IUNIT.EQ.5)WRITE(31,*)DIRT(KK)%DRYVL, DIRT(KK)%GPMOL,       &
     &            DIRT(KK)%ACVTY, DIRT(KK)%DIFTY, DIRT(KK)%HENRY
         IF(DIRT(KK)%DRYVL.GT.0.0)DIRT(KK)%DODRY=.TRUE.

!        resistence method requires molecular weight for gases
!        and only diameter for particles, however will use molecular
!        weight as the turn-on resistence flag for both particles and gases
         IF(DIRT(KK)%GPMOL.GT.0.0)DIRT(KK)%DORES=.TRUE.

!        check for consistent inputs for gravitational settling
         DIRT(KK)%DOGRV=.FALSE.
         IF(DIRT(KK)%PDIAM.GT.0.0.AND. DIRT(KK)%PDENS.GT.0.0.AND.       &
     &      DIRT(KK)%SHAPE.GT.0.0)     DIRT(KK)%DOGRV=.TRUE.

!        specified dry deposition over-rides gravitational settling
!        and resistance method specification
         IF(DIRT(KK)%DRYVL.GT.0.0)THEN
            DIRT(KK)%DOGRV=.FALSE.
            DIRT(KK)%DORES=.FALSE.
         END IF

!=>wet removal constants

         DIRT(KK)%DOWET=.FALSE.
         DIRT(KK)%WETGAS=0.0
         DIRT(KK)%WETIN=0.0
         DIRT(KK)%WETLO=0.0

         IF(IUNIT.EQ.5)THEN
          WRITE(*,*)'Wet removal constants - '
          WRITE(*,*)'-- gasses ---   ----------- particles -----------'
          WRITE(*,*)'Henrys (M/atm), In-cloud (L/L), Below-cloud (1/s)'
          WRITE(*,*)DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO
         END IF

         READ(IUNIT,*)DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO
         IF(IUNIT.EQ.5)                                                 &
     &     WRITE(31,*)DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO

!        check for consistency of wet removal definitions
         IF(DIRT(KK)%DOGAS)THEN
!           wet removal of gasses only if they are soluable
            DIRT(KK)%WETIN=0.0
            DIRT(KK)%WETLO=0.0
         ELSE
            DIRT(KK)%WETGAS=0.0
         END IF

         IF(DIRT(KK)%WETGAS.GT.0.0 .OR. DIRT(KK)%WETIN.GT.0.0 .OR.      &
     &      DIRT(KK)%WETLO .GT.0.0)     DIRT(KK)%DOWET=.TRUE.

!=>radioactive decay options (one year default)

         DIRT(KK)%DORAD=.FALSE.
         DIRT(KK)%RHALF=0.0

         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Radioactive decay half-life (days / 0=none)'
            WRITE(*,*)DIRT(KK)%RHALF
         END IF

         READ(IUNIT,*)DIRT(KK)%RHALF
         IF(IUNIT.EQ.5)WRITE(31,*)DIRT(KK)%RHALF

         IF(DIRT(KK)%RHALF.GT.0.0)DIRT(KK)%DORAD=.TRUE.

!=>deposition resuspension option

         DIRT(KK)%DOSUS=.FALSE.
         DIRT(KK)%SRATE=0.0
         IF(IUNIT.EQ.5)THEN
            WRITE(*,*)'Deposition resuspension constant (1/m)'
            WRITE(*,*)DIRT(KK)%SRATE
         END IF

         READ(IUNIT,*)DIRT(KK)%SRATE
         IF(IUNIT.EQ.5)WRITE(31,*)DIRT(KK)%SRATE
         IF(DIRT(KK)%SRATE.GT.0.0)DIRT(KK)%DOSUS=.TRUE.

!=>set species independent flags to see if any deposition options enabled

         IF(DIRT(KK)%DODRY.OR.DIRT(KK)%DORES.OR.                        &
     &      DIRT(KK)%DOGRV.OR.DIRT(KK)%DOWET.OR.DIRT(KK)%DORAD)         &
     &      CDEP=.TRUE.

         IF(DIRT(KK)%DORES)RDEP=.TRUE.
         IF(DIRT(KK)%DOSUS)SDEP=.TRUE.

  100 END DO

      RETURN
      END
