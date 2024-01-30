!======================================================================
!     module_defconc - defines concentration grid  structure and the
!     pollutant identification and release information.  It is required
!     in most CON??? routines.
!----------------------------------------------------------------------
!     $Id: module_defconc.f90,v 1.2 2007-05-03 13:09:13 skoerner Exp $
!======================================================================

MODULE module_defconc

  USE module_defsize

!     time identification structure:
!     year, month, day, hour, minute, forecast, accumulated minutes
       TYPE SSET

         INTEGER      YR
         INTEGER      MO
         INTEGER      DA
         INTEGER      HR
         INTEGER      MN
         INTEGER      IC
         INTEGER      MACC

      END TYPE SSET

!     define the geographic GRID parameters
       TYPE CSET

!        file location, identification, and unit number
                CHARACTER       DIR*100
                CHARACTER       FILE*100
                INTEGER         UNIT
!        horizontal number of grid points 
                INTEGER         NUMB_LAT
                INTEGER         NUMB_LON
!        grid spacing in fractional degrees
                REAL*8            DELT_LAT
                REAL*8            DELT_LON
!        minimum grid size (km)
                REAL*8          SIZE
!        position of lower left corner of grid
                REAL*8            X1Y1_LAT
                REAL*8            X1Y1_LON
!        type of averaging (0:std 1:snapshot)
                INTEGER         SNAP
!        vertical grid number of levels and heights (m AGL)
                INTEGER         LEVELS
                INTEGER         HEIGHT(MAXZP)
!        time to start, stop, integration interval, current time
                TYPE (SSET)     START
                TYPE (SSET)     STOP
                TYPE (SSET)     DELTA
                TYPE (SSET)     NOW
       
      END TYPE CSET

!     define the pollutant parameters
      TYPE PSET

!        identifcation label string (array element = type index)
                CHARACTER       IDENT*4
!        emission rate in units per hour
                REAL*8          QRATE
!        hours of emission
                REAL*8          QHRS
!        emission starting time
                TYPE (SSET)     START
!        logical removal flags (wet, dry, radioactive)
                LOGICAL         DOWET, DODRY, DORAD
!                               resistence, gasses, settling, resuspen
                LOGICAL         DORES, DOGAS, DOGRV, DOSUS
!        wet removal scavenging ratio within-cloud (by volume)
                REAL*8          WETIN
!        wet removal scavenging coefficient below-cloud (1/sec)
                REAL*8          WETLO
!        soluability (henry's consant) for wet removal of gases
                REAL*8          WETGAS
!        effective henry's constant for dry deposition
                REAL*8          HENRY
!        dry removal deposition velocity (m/s)
                REAL*8          DRYVL
!        gravitational settling: density (g/m3), diameter (um), shape
                REAL*8          PDENS, PDIAM, SHAPE
!        radioactive half life (days)
                REAL*8          RHALF
!        gram molecular weight
                REAL*8          GPMOL
!        pollutant activity ratio for gaseous dry deposition
                REAL*8          ACVTY
!        ratio of pollutant diffusivity to that of water vapor
                REAL*8          DIFTY
!        deposition resuspension constant
                REAL*8          SRATE

      END TYPE PSET

!======================================================================

      TYPE (CSET)        CONC(MAXGRD)
      TYPE (PSET)        DIRT(MAXTYP)

END MODULE module_defconc
