!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBL2LL           GloBaL position to Lat Lon
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONVERTS X,Y GRID POSITION TO LAT LON COORDINATES BASED UPON
!   UPON THE GRID SPACING SPECIFIED.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 09 Mar 2001 (RRD) - initial version
!
! USAGE:  CALL GBL2LL(KG,X,Y,CLAT,CLON)
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
! $Id: gbl2ll.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

!JCL:(07/12/2004) started out with HYSPLIT Ver. 45 Code

SUBROUTINE GBL2LL(KG,X,Y,CLAT,CLON)

  ! IMPLICIT NONE  !JCL: (07/09/2004) compatibility with old code
  use module_defgrid
  IMPLICIT REAL*8 (A-H,O-Z)

!  INCLUDE 'DEFGRID.INC'               ! meteorology grid and file

  INTEGER, INTENT(IN)  :: kg          ! active grid number
  REAL*8,    INTENT(IN)  :: x,y         ! grid position
  REAL*8,    INTENT(OUT) :: clat,clon   ! latlon location
  REAL*8 :: TLAT,TLON

!  COMMON /GBLGRD/ GRID, DREC, FILE

  IF(.NOT.GRID(KG)%LATLON) RETURN

  CLAT=GRID(KG)%SYNC_LAT+(Y-1.0)*GRID(KG)%REF_LAT
  IF(CLAT.GT. 90.0)CLAT= 180.0-CLAT
  IF(CLAT.LT.-90.0)CLAT=-180.0-CLAT

  CLON=GRID(KG)%SYNC_LON+(X-1.0)*GRID(KG)%REF_LON
  !CLON=AMOD(CLON,360.0)
  CLON=DMOD(CLON,DBLE(360.0))   !JCL: (07/12/2004) double precision
  IF(CLON.GT.180.0)CLON=CLON-360.0

END SUBROUTINE gbl2ll
