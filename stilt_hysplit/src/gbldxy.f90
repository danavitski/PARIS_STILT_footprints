!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBLDXY           GloBaL grid size by x,y coordinate
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DETERMINES THE GRID SPACING IN KM FOR A X,Y COORDINATE POSITION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Mar 2001 (RRD) - initial version
!
! USAGE:  CALL GBLDXY(KG,X,Y,GSX,GSY)
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
! $Id: gbldxy.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

!JCL:(07/12/2004) started out with HYSPLIT Ver. 45 Code

SUBROUTINE GBLDXY(KG,X,Y,GSX,GSY)

  use module_defgrid
  !IMPLICIT NONE  !JCL: (07/09/2004) compatibility with old code
  IMPLICIT REAL*8 (A-H,O-Z)

!  INCLUDE 'DEFGRID.INC'               ! meteorology grid and file

  INTEGER, INTENT(IN)  :: kg          ! active grid number
  REAL*8,    INTENT(IN)  :: x,y         ! grid position
  REAL*8,    INTENT(OUT) :: gsx,gsy     ! grid size at that location

  REAL*8            :: CLAT,CLON
  REAL*8, PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL*8, PARAMETER :: DEGPRD = 57.29578  ! deg per radian
  REAL*8, PARAMETER :: PI     =  3.14159

!  COMMON /GBLGRD/ GRID, DREC, FILE

!-------------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE GBL2LL(KG,X,Y,CLAT,CLON)
  ! IMPLICIT NONE  !JCL: (07/09/2004) compatibility with old code
    INTEGER, INTENT(IN)  :: kg          ! active grid number
    REAL*8,    INTENT(IN)  :: x,y         ! grid position
    REAL*8,    INTENT(OUT) :: clat,clon   ! latlon location
    END SUBROUTINE gbl2ll
  END INTERFACE
!-------------------------------------------------------------------------------

  IF(.NOT.GRID(KG)%LATLON) RETURN

  CALL GBL2LL(KG,X,Y,CLAT,CLON)

! latitude grid spacing
  GSY = (2.0*PI*REARTH) / (360.0/GRID(KG)%REF_LAT)

! longitude grid spacing
! GSX = COS(CLAT/DEGPRD)*GSY*GRID(KG)%REF_LON/GRID(KG)%REF_LAT
  GSX = DCOS(CLAT/DEGPRD)*GSY*GRID(KG)%REF_LON/GRID(KG)%REF_LAT

END SUBROUTINE gbldxy
