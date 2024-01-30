!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METGRD           computes METeorological GRiD size
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL GRID COMPUTES THE GRID SIZE FOR THE CURRENT GRID
!   IN ADDITION LOADS OTHER GRID SENSITIVE PARAMETERS SUCH AS LANDUSE
!   AND ROUGHNESS LENGTH.  VARIABLES PASSED THROUGH AS REQUIRED.
!   THE ROUTINE MUST BE CALLED EACH TIME THE SUBGRID CHANGES
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 11 May 1998 (RRD)
!                  14 Apr 1999 (RRD) - refined grid size = 0 test
!
! USAGE:  CALL METGRD(KG,LX1,LY1,NXS,NYS,GD,Z0,LU)
!   INPUT ARGUMENT LIST:
!     KG    - int grid selection index
!     LX1,LY1   - int   subgrid lower left position
!     NXS,NYS     - int subgrid dimensions
!   OUTPUT ARGUMENT LIST:
!     GD    - real      grid size array (m)
!     Z0    - real      aerodynamic roughness length (m)
!     LU    - int land-use category (1-11)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT 30 - diagnostic MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metgrd.f90,v 1.8 2008-03-26 19:16:05 tnehrkor Exp $
!
!$$$
      SUBROUTINE METGRD(KG,LX1,LY1,NXS,NYS,GX,GY,Z0,LU)

      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8    GX(NXS,NYS),GY(NXS,NYS),Z0(NXS,NYS)
      INTEGER LU(NXS,NYS)


!---------------------------------------------------------------------------------------------------
      DO JJ=1,NYS
      DO II=1,NXS

!        convert index from subgrid to full grid
         XI=II+LX1-1
         YJ=JJ+LY1-1

! JCL:(07/12/2004) added global grid code from HYSPLIT Vers. 45
      IF(GRID(KG)%LATLON)THEN
!       find position of grid node
        CALL GBL2LL(KG,XI,YJ,CLAT,CLON)

!       meters per grid cell
        CALL GBLDXY(KG,XI,YJ,GDX,GDY)

!       at pole set spacing to one grid point below pole
        IF(CLAT.EQ. 90.0) CALL GBLDXY(KG,XI,(YJ-1.0),GDX,GDY)
        IF(CLAT.EQ.-90.0) CALL GBLDXY(KG,XI,(YJ+1.0),GDX,GDY)

!       convert to meters
        GX(II,JJ)=GDX*1000.0
        GY(II,JJ)=GDY*1000.0

      ELSE
!       meters per grid cell
        GX(II,JJ)=CGSZXY(GRID(KG)%GBASE,XI,YJ,GRID(KG))*1000.0

!       accounts for error in grid conversion at pole points
        IF(GX(II,JJ).LT.1.0)                                            &
     &     GX(II,JJ)=CGSZXY(GRID(KG)%GBASE,XI+0.5,YJ+0.5,GRID(KG))*1000.0

!       conformal projection
        GY(II,JJ)=GX(II,JJ)

!       find position of grid node
        CALL CXY2LL(GRID(KG)%GBASE,XI,YJ,CLAT,CLON, GRID(KG)%proj)
      END IF

! JCL:
!        WRITE(45,*) 'GX,GY in METGRD:',II,JJ,GX(II,JJ),GY(II,JJ)

!        uses units 60 (land-use) and 62 (roughness length)
!        to load data into each node from lat/lon based input file
         CALL SFCINP(CLAT,CLON,Z0(II,JJ),LU(II,JJ))

! JCL:   occasionally, Z0=0, and model crashes b/c dividing by 0
!            so reset Z0 to a very small value (0.01)
         IF(Z0(II,JJ).EQ.0) Z0(II,JJ)=0.01

      END DO
      END DO

      WRITE(30,'(A,2I3)')                                               &
     &   ' NOTICE metgrd: subgrid buffers loaded at - ',LX1,LY1

      RETURN
      END
