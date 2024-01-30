!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBL2XY           GloBaL position to XY cordinates
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONVERTS LAT LON POSITION TO X Y GRID COORDINATES BASED UPON
!   UPON THE GRID SPACING SPECIFIED.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Mar 2001 (RRD) - initial version
!
! USAGE:  CALL GBL2XY(KG,CLAT,CLON,X,Y)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
! $Id: gbl2xy.f90,v 1.5 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

SUBROUTINE GBL2XY(KG,CLAT,CLON,X,Y)

  use module_defgrid
  IMPLICIT NONE

!  INCLUDE 'DEFGRID.INC'                                     ! meteorology grid and file

  INTEGER, INTENT(IN)  :: kg                                ! active grid number
  REAL*8,  INTENT(IN)  :: clat,clon                         ! latlon location
  REAL*8,  INTENT(OUT) :: x,y                               ! grid position

  REAL*8               :: tlat,tlon

!  COMMON /GBLGRD/ GRID, DREC, FILE

  IF(.NOT.GRID(KG)%LATLON) RETURN

! Grid system is simply defined as the number of grid points
! from the corner point at 1,1 using an even lat-lon increment
! for the x and y directions. Grid distances are computed
! where needed according to the latitude of the grid point

  TLAT=CLAT
  IF(TLAT.GT. 90.0)TLAT= 180.0-TLAT
  IF(TLAT.LT.-90.0)TLAT=-180.0-TLAT
  Y=1.0+(TLAT-GRID(KG)%SYNC_LAT)/GRID(KG)%REF_LAT

  TLON=CLON
  IF(TLON.LT.0.0)  TLON=360.0+TLON
  IF(TLON.GT.360.0)TLON=TLON-360.0
  TLON=TLON-GRID(KG)%SYNC_LON
  IF(TLON.LT.0.0)TLON=TLON+360.0
  X=1.0+TLON/GRID(KG)%REF_LON

END SUBROUTINE GBL2XY
