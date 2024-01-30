!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBHOR           STaBility HORizonal
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:98-12-17
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   COMPUTES THE HORIZONTAL MIXING COEFFICIENT BASED UPON THE
!   SMAGORINSKY FORMULATION OF THE HORIZONTAL WIND FIELD DEFORMATION
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 12 Dec 1998 - RRD
!
! USAGE:  CALL STBHOR(NXS,NYS,NLVL,GD,U,V,H)
!   INPUT ARGUMENT LIST:
!     NXS,NYS     - int horizontal subgrid dimensions
!     NLVL  - int number of output levels
!     GD    - real      horizontal grid spacing
!     U,V         - real      horizontal wind components
!   OUTPUT ARGUMENT LIST:
!     H           - real      3d horizontal mixing
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: stbhor.f90,v 1.5 2007-02-16 17:53:48 tnehrkor Exp $
!
!$$$

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!     SUBROUTINE STBHOR(NXS,NYS,NLVL,GD,U,V,H)
      SUBROUTINE STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)

!     dimension information

!     3D meteorological variables
      REAL*8 U(NXS,NYS,NZM), V(NXS,NYS,NZM), H(NXS,NYS,NZM)

!     2D meteorological variables - grid size array
!     REAL*8 GD(NXS,NYS)
      REAL*8 GX(NXS,NYS),GY(NXS,NYS)

!     process each node on subgrid
      DO J=1,NYS-1
      DO I=1,NXS-1
      DO K=1,NLVL

!        to obtain the horizontal diffusivity (from Smagorinsky, 1963)
!        from velocity deformation (shearing and tension stress)
!        centered differences and all horizontal velocities in grid/min
!        (u2-u1)/(2 delta-x), where x = 1 grid length

!        check edge limits (P-plus and M-minus)
         IP=MIN(I+1,NXS)
         JP=MIN(J+1,NYS)
         IM=MAX(I-1,1)
         JM=MAX(J-1,1)

!        deformation (1/min)
         DUDY=U(I+1,J+1,K)-U(I+1,J,K)
         DVDX=V(I+1,J+1,K)-V(I,J+1,K)
         DUDX=U(I+1,J+1,K)-U(I,J+1,K)
         DVDY=V(I+1,J+1,K)-V(I+1,J,K)
         DEFT=0.5*SQRT((DUDY+DVDX)*(DUDY+DVDX)+(DUDX-DVDY)*(DUDX-DVDY))

!        deformation in 1/min => mixing m2/sec
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        H(I+1,J+1,K)=0.014*GD(I+1,J+1)*GD(I+1,J+1)*DEFT/60.0
         H(I+1,J+1,K)=0.014*GX(I+1,J+1)*GY(I+1,J+1)*DEFT/60.0

      END DO
      END DO
      END DO

      RETURN
      END
