!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METWND           METeorologically derived WiND
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICALLY DERIVED WIND DETERMINES WHICH VARIABLE IS TO BE
!   ANALZED FOR ITS SLOPE DEPENDING UPON THE VERTICAL VELOCITY OPTION
!   SELECTED ON INPUT - THESE OPTIONS ASSUME VARIOUS PROPERTIES ARE
!   CONSTANT ALONG THE TRAJECTORY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,
!              P1,P2,T1,T2,D1,D2)
!   INPUT ARGUMENT LIST:
!     KVEL  - int vertical motion flag (0-leave alone)
!     NXS,NYS     - int dimensions of sub-grid
!     NZM   - int vertical grid dimension
!     NLVL  - int number of data levels
!     DM    - real      time between obs
!     ZSG   - real      sigma levels
!     U,V,W,etc - real  meteorological variables (see advpnt for list)
!   OUTPUT ARGUMENT LIST:
!     W           - real      vertical velocity
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metwnd.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

      SUBROUTINE METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,             &
     &   P1,P2,T1,T2,D1,D2)

      IMPLICIT REAL*8 (A-H,O-Z)

!     velocity components
      REAL*8   U(NXS,NYS,NZM), V(NXS,NYS,NZM), W(NXS,NYS,NZM)

!     variables at old (1) and new (2) times to maintain on surface
      REAL*8   P1(NXS,NYS,NZM),P2(NXS,NYS,NZM),T1(NXS,NYS,NZM),         &
     &       T2(NXS,NYS,NZM),D1(NXS,NYS,NZM),D2(NXS,NYS,NZM)

!     internal model sigma levels
      REAL*8   ZSG(NLVL)

      IF(KVEL.EQ.0)THEN
         RETURN
      ELSEIF(KVEL.EQ.1)THEN
!        constant pressure
         CALL METSLP(NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,P1,P2)
      ELSEIF(KVEL.EQ.2)THEN
!        constant temperature (isentropic)
         CALL METSLP(NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,T1,T2)
      ELSEIF(KVEL.EQ.3)THEN
!        constant density
         CALL METSLP(NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,D1,D2)
      ELSEIF(KVEL.EQ.4)THEN
!        constant sigma level
         DO K=1,NLVL
         DO J=1,NYS
         DO I=1,NXS
            W(I,J,K)=0
         END DO
         END DO
         END DO
      END IF

      RETURN
      END
