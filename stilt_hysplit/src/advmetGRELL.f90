!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONCALL           calls Convective redistribution using GRELL scheme
!   PRGMMR:    CHRISTOPH GERBIG   DATE:09-23-03
!
! ABSTRACT:
!   extracts convective fluxes for gridcell and calls CGRELL convective
!   redistribution of particles using updraft & downdraft mass
!   Comparable to ADVMET
! USAGE:  CALL ADVMETGRELL(X,Y,CFXUP1, CFXUP2, CFXDN1, DFXUP1, DFXUP2, EFXUP1,
!                          EFXUP2, DFXDN1, EFXDN1, RAUP1, RAUP2, RADN1)
!   INPUT ARGUMENT LIST:
!     x,y       - real  horizontal particle position (grid units)
!     ?FX*      - real  Cloud/Entrainment/Detrainment Mass Fluxes
!                       (Updraft, downdraft; 1:deep; 2: shallow)
!                       UNITS: kg/m^2s; convention: upward=positive
!                       TO GET MASSFLUX (kg/s):
!                       MASSFLUX=?FX* * DX*DY (up/dndraft & hor. fluxes)
!                       numbers represent ro*w*a w/ a=fractional coverage
!                       to get mass flux[kg/s]: ?FX* * area (gridcell)
!     RA*       - real  updraft area coverage (fract.),
!                       (Updraft, downdraft; 1:deep; 2: shallow)
!     NXS,NYS   - int   horizontal grid dimensions
!     NLVL      - int   number of vertical levels
!   OUTPUT ARGUMENT LIST:
!     GBLMET COMMON WITH VARIABLES DEFINED IN DEFMETO.INC
!
! $Id: advmetGRELL.f90,v 1.6 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$
! CHG (09/23/03) Radius up/downdraft (RAUP1,RAUP2,RADN1) not yet from RAMS, so assign small value
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
      SUBROUTINE ADVMETGRELL(X,Y,CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,    &
     &  EFXUP1,EFXUP2,DFXDN1,EFXDN1,TKEN,NXS,NYS,NLVL,GLOBAL,NXP,NYP)

      use module_defmeto
      IMPLICIT REAL*8 (A-H,O-Z)
!     dimension information
!     contains meteorological summary at last advection point
!      INCLUDE 'DEFMETO.INC'
! JCL:(07/12/2004) added cyclic boundary condition flag
      LOGICAL GLOBAL
! CHG (09/23/03) Radius up/downdraft not yet from RAMS, so assign small value
      REAL*8 CFXUP1(NXS,NYS,NZM),CFXUP2(NXS,NYS,NZM),                   &
     &     CFXDN1(NXS,NYS,NZM),DFXUP1(NXS,NYS,NZM),DFXUP2(NXS,NYS,NZM), &
     &     EFXUP1(NXS,NYS,NZM),EFXUP2(NXS,NYS,NZM),DFXDN1(NXS,NYS,NZM), &
     &     EFXDN1(NXS,NYS,NZM)
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
      REAL*8 TKEN(NXS,NYS,NZM)

!      COMMON /GBLMET/ METO

      XX=DNINT(X)
      YY=DNINT(Y)
      DO KL=1,NLVL
!     set vertical interpolation point to index position
        ZK=KL
        CALL ADVINT(CFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%CFXUP1(KL)=VAR1
        CALL ADVINT(CFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%CFXUP2(KL)=VAR1
        CALL ADVINT(CFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%CFXDN1(KL)=VAR1
        CALL ADVINT(DFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%DFXUP1(KL)=VAR1
        CALL ADVINT(DFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%DFXUP2(KL)=VAR1
        CALL ADVINT(EFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
! CHG(09/24/03) need to shift?
                             !no shift
        METO%EFXUP1(KL)=VAR1
!        IF(KL.LT.NLVL)METO%EFXUP1(KL+1)=VAR1 !leaves EFXUP1(1) zero...
        CALL ADVINT(EFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%EFXUP2(KL)=VAR1
        CALL ADVINT(DFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%DFXDN1(KL)=VAR1
        CALL ADVINT(EFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
     &              VAR1)
        METO%EFXDN1(KL)=VAR1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
        CALL ADVINT(TKEN,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,    &
     &              VAR1)
        METO%TKEN(KL)=VAR1
!      IF(ZK.EQ.2)Write(45,*)'advmetGRELL CFXDN1',CFXDN1(7,7,1:10)
!      IF(ZK.EQ.2)Write(45,*)'advmetGRELL EFXDN1',EFXDN1(7,7,1:10)
!      IF(ZK.EQ.2)Write(45,*)'advmetGRELL XX,YY',XX,YY
!      IF(ZK.EQ.2)Write(45,*)'advmetGRELL M%CFXDN1',METO%CFXDN1(1:10)
!      IF(ZK.EQ.2)Write(45,*)'advmetGRELL M%EFXDN1',METO%EFXDN1(1:10)
      END DO

      RETURN
      END
