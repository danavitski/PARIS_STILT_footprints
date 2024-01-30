!======================================================================
!     module_defgrid - defines meteorological grid location, variables
!     contained during each time, and file structure. It is required
!     in most MET??? routines.
!----------------------------------------------------------------------
!     $Id: module_defgrid.f90,v 1.3 2008-03-27 14:21:32 tnehrkor Exp $
!======================================================================

MODULE module_defgrid

  USE module_defsize
  USE map_utils, only : proj_info

  real*8, parameter :: vmiss=-999.99, vmissle=-999.

!     define the geographic GRID parameters
      TYPE GSET

!        real array that defines grid conversions
         REAL*8                 GBASE(15)

         REAL*8                 POLE_LAT
         REAL*8                 POLE_LON
         REAL*8                 REF_LAT
         REAL*8                 REF_LON
         REAL*8                 SIZE
         REAL*8                 ORIENT
         REAL*8                 TANG_LAT
         REAL*8                 SYNC_XP
         REAL*8                 SYNC_YP
         REAL*8                 SYNC_LAT
         REAL*8                 SYNC_LON
         REAL*8                 DUMMY

!        grid-points in x, in y, and data levels (incl sfc)
         INTEGER        NX
         INTEGER        NY
         INTEGER        NZ

!        character model data identification
         INTEGER        NUMBER
         CHARACTER      MODEL_ID*4

! JCL:(07/06/2004) needed for global grid
         LOGICAL        LATLON   ! flag if input grid is latlon
         LOGICAL        GLOBAL   ! flag if latlon subgrid is global

         type (proj_info) :: proj  !WPS mapping routines data structure
         real, dimension(:,:), pointer :: mapfactor !mapfactor array for WPS

      END TYPE GSET

!======================================================================

!     define structure of Data RECords within a file
      TYPE MSET   

!        checksum for each field
                INTEGER                 CHK_SUM(MVAR,MLVL)  
!        height of each data level
                REAL*8                  HEIGHT(MLVL)
!        variables per data level
                INTEGER                 NUM_VARB(MLVL)
!        data type flag 1:old 2:new
                INTEGER                 TYPE
!        1:sigma 2:pressure 3:terrain
                INTEGER                 Z_FLAG
!        records per time period
                INTEGER                 REC_PER
!        minutes between times
                INTEGER                 DELTA
!        number of records offset
                INTEGER                 OFFSET
!        accumulation cycle (minutes)
                INTEGER                 ACYCLE
! CHG:(11/20/01) convective precipitation flag
                LOGICAL                 CFLG
! CHG:(11/20/01) total cloud flag
                LOGICAL                 TCLF
! CHG:(12/04/01) low cloud flag
                LOGICAL                 LCLF
! CHG:(11/20/01) shortw. downw. radiation flag
                LOGICAL                 RADF
! CHG:(22/01/03) soil moisture flag
                LOGICAL                 SLMF
!        upper humidity, w-velocity 
                LOGICAL                 QFLG,RFLG(MLVL),WFLG(MLVL)
!        surface flux flags 
!               LOGICAL                 SFLG,EFLG,STAR
!        surface wind and temp flags
                LOGICAL                 UFLG,TFLG
! JCL:(5/7/01)  new surface flags that Draxler implemented
!        surface flux flags
                LOGICAL                 EFLX,HFLX,UFLX,USTR,TSTR
! CHG:(23/11/06) TKE flag
                LOGICAL                 TKEF
                !        surface terrain height, pressure
                LOGICAL                 SHGT,PRSS
!        character identification
                CHARACTER               VARB_ID(MVAR,MLVL)*4

      END TYPE MSET

!======================================================================

!     defines a time structure
      TYPE TSET

         INTEGER        YR
         INTEGER        MO
         INTEGER        DA
         INTEGER        HR
         INTEGER        MN
         INTEGER        IC
!        accumulated minutes from 1 jan
         INTEGER        MACC

      END TYPE TSET

!======================================================================

!     defines structure of meteorological data FILE
      TYPE FSET

!        beginning time period of file
                TYPE (TSET)             FIRST
!        ending time period of file
                TYPE (TSET)             LAST
!        unit number assigned 
                INTEGER                 KUNIT          
!        pointer to time period 1 or 2
                INTEGER                 PERIOD
!        record length (bytes) 
                INTEGER                 REC_LEN
!        number of last record in file
                INTEGER                 ENDREC
!        directory of meteorology
                CHARACTER               DIR*100
!        name of meteorology file
                CHARACTER               METEO*100

      END TYPE FSET

!======================================================================

      TYPE (GSET)     GRID(MGRD)
      TYPE (MSET)     DREC(MGRD)
      TYPE (FSET)     FILE(MGRD,2)

END MODULE module_defgrid
