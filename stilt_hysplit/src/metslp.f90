!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METSLP           METeorological SLoPe of a surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL SLOPE COMPUTES THE SLOPE OF METEOROLOGICAL VARIABLE
!   AND THEIR LOCAL TIME DERIVATIVE AND DETERMINES THE VERTICAL VELOCITY
!   REQUIRED TO STAY ON THE SURFACE FOR THAT VARIABLE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 - RRD
!
! USAGE:  CALL METSLP(NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,X1,X2)
!   INPUT ARGUMENT LIST:
!     NXS,NYS     - int dimensions of sub-grid
!     NZM   - int vertical grid dimension
!     NLVL  - int number of data levels
!     DM    - real      time between obs
!     ZSG   - real      sigma levels
!     U,V   - real      meteorological wind components
!     X1    - real      data mapping variable (last time)
!     X2    - real      data mapping variable (next time)
!   OUTPUT ARGUMENT LIST:
!     W           - real      vertical motion (ds/dt)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metslp.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

      SUBROUTINE METSLP(NXS,NYS,NZM,NLVL,DM,ZSG,U,V,W,X1,X2)

      IMPLICIT REAL*8 (A-H,O-Z)

!     velocity components, variable to maintain surface
      REAL*8    U(NXS,NYS,NZM), V(NXS,NYS,NZM), W(NXS,NYS,NZM),         &
     &       X1(NXS,NYS,NZM),X2(NXS,NYS,NZM)

!     internal model sigma levels
      REAL*8   ZSG(NLVL)

      DO I=2,NXS-1
      DO J=2,NYS-1

         DO K=1,NLVL

!              compute horizontal gradients according to level
               IF(K.EQ.1)THEN
                  DELX=X2(I,J,K+1)-X2(I,J,K)
                  DELZ=ZSG(K+1)-ZSG(K)
               ELSEIF(K.EQ.NLVL)THEN
                  DELX=X2(I,J,K)-X2(I,J,K-1)
                  DELZ=ZSG(K)-ZSG(K-1)
               ELSE
                  DELX=X2(I,J,K+1)-X2(I,J,K-1)
                  DELZ=ZSG(K+1)-ZSG(K-1)
               END IF

!              vertical slope of surface (per sigma)
               VERT=DELX/DELZ

!              horizontal slope of surface (grid/min)
               DXD=0.5*U(I,J,K)*(X2(I+1,J,K)-X2(I-1,J,K))
               DYD=0.5*V(I,J,K)*(X2(I,J+1,K)-X2(I,J-1,K))
               HORZ=DXD+DYD

!              local derivative (per min)
               TIME=(X2(I,J,K)-X1(I,J,K))/DM

!              w velocity in sigma/min to maintain surface
               IF(VERT.NE.0.0)THEN
                  WVEL=-(TIME+HORZ)/VERT
!                 restrict maximum to 0.005 sig/min (about 1 m/s)
                  WNEW=DMIN1(DBLE(0.005),ABS(WVEL))
                  W(I,J,K)=SIGN(WNEW,WVEL)
               END IF
         END DO

      END DO
      END DO

      RETURN
      END
