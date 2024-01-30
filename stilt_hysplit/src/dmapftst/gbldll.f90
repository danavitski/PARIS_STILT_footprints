!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBLDLL           GloBaL grid size by LAT LON coord
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DETERMINES THE GRID SPACING IN KM FOR A LAT LON POSITION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 20 Apr 2001 (RRD) - initial version
!
! USAGE:  CALL GBLDLL(KG,CLAT,CLON,GSX,GSY)
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
!$$$

!JCL:(07/12/2004) started out with HYSPLIT Ver. 45 Code

SUBROUTINE GBLDLL(KG,CLAT,CLON,GSX,GSY)

  ! IMPLICIT NONE  !JCL: (07/09/2004) compatibility with old code
  IMPLICIT REAL*8 (A-H,O-Z)

  INCLUDE 'DEFSIZE.INC' ! JCL: (07/09/2004) need array size info
  INCLUDE 'DEFGRID.INC'               ! meteorology grid and file

  INTEGER, INTENT(IN)  :: kg          ! active grid number
  REAL*8,    INTENT(IN)  :: clat,clon   ! grid position
  REAL*8,    INTENT(OUT) :: gsx,gsy     ! grid size at that location

  REAL*8, PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL*8, PARAMETER :: DEGPRD = 57.29578  ! deg per radian
  REAL*8, PARAMETER :: PI     =  3.14159

  COMMON /GBLGRD/ GRID, DREC, FILE

!-------------------------------------------------------------------------------

  IF(.NOT.GRID(KG)%LATLON) RETURN

! latitude grid spacing
  GSY = (2.0*PI*REARTH) / (360.0/GRID(KG)%REF_LAT)

! longitude grid spacing
! GSX = COS(CLAT/DEGPRD)*GSY*GRID(KG)%REF_LON/GRID(KG)%REF_LAT
  GSX = DCOS(CLAT/DEGPRD)*GSY*GRID(KG)%REF_LON/GRID(KG)%REF_LAT

END SUBROUTINE gbldll
