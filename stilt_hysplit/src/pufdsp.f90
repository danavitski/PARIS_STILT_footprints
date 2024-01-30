!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFDSP           PUFf DiSPersion subgrid puff dispersion
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF DISPERSION ASSUMES LINEAR GROWTH OF HORIZONTAL DISTRIBUTION
!   SQUARE ROOT VERTICAL GROWTH. GROWTH RATES ARE COMPUTED FROM THE
!   STANDARD DEVIATION OF THE TURBULENT VELOCITY.  TURBULENT VELOCITY
!   VARIANCES MAY BE OBTAINED FROM THE METEOROLOGICAL MODEL DIRECTLY
!   PARAMETERIZED FROM THE HORIZONTAL AND VERTICAL DIFFUSIVITY,
!   SIGMA^2 = H / TL.  VERTICAL GROWTH MAY REQUIRE FINER TIME STEPS
!   MAINTAIN COMPUTATIONAL STABILITY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 18 Aug 1998 (RRD) - isotroptic turbulence option
!                 20 Apr 1999 (RRD) - added terrain compression
!                 26 Aug 1999 (RRD) - pass vertical grid in common block
!
! USAGE:  CALL PUFDSP(HMIX,DT,ZMDL,ZSFC,NLVL,VMIX,ZSG,MASS,ZPOS,SIGH,
!              SIGV,HDWP,ZNDX,ISOT)
!   INPUT ARGUMENT LIST:
!     HMIX  - real horizontal diffusivity (m2/s)
!     DT    - real horizontal time step from advection (min)
!     ZMDL  - real top of computational domain (meters)
!     ZSFC  - real surface terrain elevation (m)
!     NLVL  - int  number of levels in subgrid
!     VMIX  - real vertical mixing profile (m2/s)
!     ZSG   - real internal model sigma levels
!     MASS  - real mass of pollutant (arbitrary units)
!     HDWP  - int  Horizontal distribution within pollutant
!     ZNDX  - real fraction vertical index position
!     ISOT  - int  flag to set isotropic turbulence (1=yes 0=no)
!   OUTPUT ARGUMENT LIST:
!     ZPOS  - real puff center height (sigma)
!     SIGH,SIGV   - real horiz (meters) and vert sigma (sigma)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: pufdsp.f90,v 1.4 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

      SUBROUTINE PUFDSP(HMIX,DT,ZMDL,ZSFC,NLVL,VMIX,ZSG,MASS,           &
     &   ZPOS,SIGH,SIGV,HDWP,ZNDX,ISOT)

      use module_defsize
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes

!     puff age, distribution type
      INTEGER   HDWP
!     pollutant mass array (mass species)
      REAL*8   MASS(MAXDIM)
!     vertical mixing, sigma levels
      REAL*8   VMIX(NLVL), ZSG(NLVL)

!     common block to pass parameters for vertical grid
      COMMON /ZZTOKK/ AA,BB,CC

!     Lagrangian time scales vertical and horizonatl (sec)
      DATA VSCALE/ 100.0 /, HSCALE/ 10800.0 /

!     vertical extent,    and sqrt(2)
      DATA SIGR/ 1.54 /,  root2/ 1.414 /

!     full particles not calculated in this routine
      IF(HDWP.EQ.0)RETURN

!=>horizontal puff diffusion

!     horizontal turbulent velocity sigma (m/s) from deformation
      DSH=DMAX1(DBLE(0.2),DMIN1(DBLE(2.), SQRT(HMIX/HSCALE) ))

!     replace horizontal turbulence with vertical for isotropic option
      IF(ISOT.EQ.1)THEN

!        find the puff vertical index number for center position
         KZ=MAX0(1, MIN0(NLVL, NINT(ZNDX) ))
         KT=MIN0(KZ+1, NLVL)
         KB=KT-1

         AMIX=SQRT(VMIX(KT)*VMIX(KB))
         DSH=SQRT(AMIX/VSCALE)
      END IF

!     linear horizontal diffusion summation (m/sec => meters)
! JCL:take ABS of DT, since SIGH should be a positive number
      SIGH=SIGH+ROOT2*DSH*ABS(DT)*60.0
!      SIGH=SIGH+ROOT2*DSH*DT*60.0

!     particle in vertical return now
      IF(HDWP.EQ.3.OR.HDWP.EQ.4)RETURN

!=>vertical puff diffusion

!     find the puff vertical index number for center position
      KZ=MAX0(1, MIN0(NLVL, NINT(ZNDX) ))
      KT=MIN0(KZ+1, NLVL)
      KB=KT-1

!     vertical turbulent velocity variance (m2/s2)
      WVVT=VMIX(KT)/VSCALE
      WVVB=VMIX(KB)/VSCALE

!     determine the required vertical dispersion time step (min)
      DELZ=(ZMDL-ZSFC)*(ZSG(KB)-ZSG(KT))
! JCL:take ABS of DT to make sure that DELT > 0
      DELT=DMIN1(ABS(DT),  0.125*DELZ*DELZ/WVVT/VSCALE/60.0)
!      DELT=AMIN1(DT,   0.125*DELZ*DELZ/WVVT/VSCALE/60.0)
      DELT=DMIN1(DELT,0.125*DELZ*DELZ/WVVB/VSCALE/60.0)
      DELT=AMAX0(1,NINT(DELT))

!     round down time step to even multiple of DT
! JCL:take ABS of DT
      DO WHILE (MOD(INT(ABS(DT)),INT(DELT)).NE.0.AND.INT(DELT).GT.1)
!      DO WHILE (MOD(INT(DT),INT(DELT)).NE.0.AND.INT(DELT).GT.1)
         DELT=DELT-1.0
      END DO

!     go through iterations to match external DT
! JCL: take ABS of DT
      DO KNUM=1,NINT(ABS(DT)/DELT)
!      DO KNUM=1,NINT(DT/DELT)

!        index at bottom and top of puff
         SGT=DMAX1(ZPOS-SIGR*SIGV,DBLE(0.0))
         SGB=DMIN1(ZPOS+SIGR*SIGV,DBLE(1.0))

!        compute vertical index (at point just below and just above)
!        from integer index based upon quadratic relation between
!        height (or sigma) and array index position

!        height agl based on zsfc=0
!        index equation based on zsfc=0
! use the ind_zsg routine instead to allow for externally specified zsg levels:
!!$         ZZ=ZMDL*(1.0-DMIN1(DBLE(1.),SGB))
!!$         ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
         call ind_zsg(zmdl,zsg,nlvl,sgb,zx,aa,bb,cc)
         KB=MIN(MAX(1,INT(ZX)),NLVL)

!        height agl based on zsfc=0
!        index equation based on zsfc=0
! use the ind_zsg routine instead to allow for externally specified zsg levels:
!!$         ZZ=ZMDL*(1.0-DMIN1(DBLE(1.0),SGT))
!!$         ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
         call ind_zsg(zmdl,zsg,nlvl,sgt,zx,aa,bb,cc)
         KT=MIN(MAX(1,NINT(ZX)),NLVL)

!        vertical turbulent velocity variance (m2/s2)
         WVVT=VMIX(KT)/VSCALE
         WVVB=VMIX(KB)/VSCALE

!        square of vertical diffusion rate of at bottom and top
         DS2T=2.0*WVVT*VSCALE
         DS2B=2.0*WVVB*VSCALE

!        compute new puff sigma's normalized vertical coordinate
         VSIGT=SQRT(SIGV*SIGV+DS2T*DELT*60.0/(ZMDL-ZSFC)/(ZMDL-ZSFC))
         VSIGB=SQRT(SIGV*SIGV+DS2B*DELT*60.0/(ZMDL-ZSFC)/(ZMDL-ZSFC))

!        updated bottom and top of puff
         SGB=DMIN1(ZPOS+SIGR*VSIGB,DBLE(1.0))
         SGT=DMAX1(ZPOS-SIGR*VSIGT,DBLE(0.0))

!        fraction lost through upper boundary
!        OUT=AMAX1(ZSG(NLVL)-SGT,0.0)/(SGB-SGT)
         SGT=DMAX1(SGT,ZSG(NLVL))

!        update sigma and vertical position
         SIGV=(SGB-SGT)/SIGR/2.0
         ZPOS=0.5*(SGB+SGT)

!        mass lost out of upper lid
!        IF(OUT.GT.0.0)THEN
!           MASS(1)=MASS(1)*(1.0-OUT)
!           MM=MAXDIM
!           DO WHILE(MM.GT.1)
!              MASS(MM)=MASS(MM)*(1.0-OUT)
!              MM=MM-1
!           END DO
!        END IF

!     time step loop
      END DO

      RETURN
      END
