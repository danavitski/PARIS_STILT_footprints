!======================================================================
!     module_defsize - sets dimensions for many of the customized 
!     include files as well as other program variables.  Should be
!     optimized for concentration, trajectory, utility applications.
!----------------------------------------------------------------------
!     $Id: module_defsize.f90,v 1.2 2007-05-03 13:09:13 skoerner Exp $
!======================================================================

MODULE module_defsize

!     maximum grid dimension for meteorological sub-grid
   INTEGER, PARAMETER :: NZM=90, NXM=281, NYM=221

!     maximum levels, variables, and grids for structure definition
!     of meteorological data files
   INTEGER, PARAMETER :: MLVL=90, MVAR=35, MGRD=17

!======================================================================

!     pollutant maximum limits:
!     MAXPAR - total # of particles in calculation
!     MAXDIM - # of species associated with each particle
!     MAXTYP - # of species associated with different particles
   INTEGER, PARAMETER :: MAXPAR=10000, MAXDIM=1, MAXTYP=2

!     concentration grid dimensions limits:
!     MAXGRD - # number of simultaneous concentration grids
!     MAXXP  - # of grid points in longitude
!     MAXYP  - # of grid points in latitude
!     MAXZP  - # of levels in the vertical
   INTEGER, PARAMETER :: MAXGRD=2, MAXXP=300, MAXYP=300, MAXZP=10

!======================================================================

!     maximum number of sources per simulation           
   INTEGER, PARAMETER :: MLOC=1000

END MODULE module_defsize
