!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBSND           STaBility SouNDing computes vertical mixing
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY SOUNDING COMPUTES THE VERTICAL MIXING COEFFICIENT
!   PROFILE FROM SURFACE STABILITY PARAMETERS AND BULK RICHARDSON
!   NUMBER GIVEN MONIN-OBUKHOV LENGTH AND FRICTION VELOCITY. ROUTINE
!   ANALYZES SOUNDING TO DETERMINE MIXED LAYER DEPTH.  USING
!   SIMILARITY FLUX PROFILE RELATIONSHIPS, THE VERTICAL DIFFUSIVITY
!   IS COMPUTED WITHIN THE PBL AND IN THE REMAINDER OF THE ATMOSPHERE
!   USING STANDARD MIXING LENGTH THEORY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 13 Mar 1997 - RRD
!                  18 Aug 1998 - RRD: isotropic turbulence no BL average
!
! USAGE:  CALL STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,
!                     NL,UU,VV,TT,ZZ,XX)
!   INPUT ARGUMENT LIST:
!     ISOT  - int flag to indicate isotropic turbulence (0-no 1-yes)
!     KSFC  - int index that reprsents top of surface layer
!     TSTR  - real      friction temperature (K)
!     USTR  - real      friction velocity (m/s)
!     WSTR  - real      convective velocity scale (m/s)
!     SLEN  - real      Obukhov length scale (m)
!     ZMIX  - real      boundary layer depth (m)
!     NL    - int number of sigma levels
!     UU,VV - real      array horizontal wind components
!     TT    - real      array virtual potential temperature (pot K)
!     ZZ    - real      array height at levels (m)
!   OUTPUT ARGUMENT LIST:
!     XX    - real      array vertical mixing coefficient (m2/s)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: stbsnd.f90,v 1.3 2005-12-14 17:05:59 tnehrkor Exp $
!
!$$$

! JCL:add Lagrangian timescale TTL & stddev of vertical velocity SIGMAW
!     as arguments, so can be output to PRFCOM, which calls STBSND
! JCL:add the roughness length as argument
      SUBROUTINE STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,             &
     &                  NL,UU,VV,TT,ZZ,XX,TTL,SIGMAW,Z0)
!      SUBROUTINE STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,
!     :                  NL,UU,VV,TT,ZZ,XX)

      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8   UU(NL),VV(NL),TT(NL),ZZ(NL),XX(NL)
! JCL:
      REAL*8   TTL(NL),SIGMAW(NL)

! JCL:z/h, the alt divided by mixed-layer height--will be used often
      REAL*8   ZH

! JCL:roughness length (m)
      REAL*8   Z0

!          gravity (m/s2)  Karman's,   Neut Prandtl
      DATA GRAV/9.80616/,  VONK/0.4/,  PRN/0.923/

!     Mixing minimum,    maximum
      DATA   VMIN/0.01/, VMAX/200.0/

!     coefficients for Beljaars-Holtslag data (1991)
      DATA A/1.0/, B/0.66667/, C/5.0/, D/0.35/

!     coefficients for mixing length scales
      DATA A1/0.2828E-03/,  A2/0.8049/,     A3/1.6583/,                 &
     &     A4/0.5090E-02/,  A5/-1.0063E-03/

! JCL:typical value for Coriolis parameter in mid-lats (s-1)
      DATA FF/0.0001/

!     summation for average mixing
      VSUM=0.0
      WGHT=0.0

!     total flux for diabatic Prandtl number calculation
      WM=(USTR*USTR*USTR+0.6*WSTR*WSTR*WSTR)**0.3333

!     top of surface layer (minimum defined for data)
      ZSFC=DMAX1(ZZ(KSFC),0.1*ZMIX)

!=>develop vertical mixing profile for top of each layer

      DO K=1,NL
!=>level is within the PBL

         IF(ZZ(K).LE.ZMIX)THEN
!           recompute stablility for level height
            SL=ZZ(K)/SLEN
            SL=DMAX1(DBLE(-2.0),DMIN1(DBLE(10.0),SL))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! JCL:      use more sophisticated parameterization for Lagrangian
!           timescale, and directly calculate sigmaw, instead of
!           deriving it from eddy diffusivity
!           based on Hanna(1982) & Stohl's FLEXPART Ver3.1
! JCL:z/h, the alt divided by mixed-layer height--will be used often
            ZH=ZZ(K)/ZMIX
!***********STABLE***********
            IF(SL.GE.0.0)THEN
! JCL:        parameterization for stable conditions
              SIGMAW(K)=1.3*USTR*(1.0-ZH)
! CHG(10/30/02): avoid SIGMAW becoming zero
              SIGMAW(K)=DMAX1(SIGMAW(K),DBLE(0.00001))
! JCL:(10/30/02) wrong exponent used previously from FLEXPART V3.1 Manual
!             TTL(K)=0.1*(ZMIX/SIGMAW(K))*(ZH**0.5)
              TTL(K)=0.1*(ZMIX/SIGMAW(K))*(ZH**0.8)
! JCL:(11/1/02) TL has singularity at ZZ(K)==ZMIX that pumps particles out of ML, so to prevent this assign TL of 100.0 at this leve
              IF(ZH.EQ.1.0)TTL(K)=100.0

!***********UNSTABLE*********
            ELSE

! JCL:      parameterization for unstable conditions
! CHG(10/30/02): use real Hanna, not the modified one found in FLEXPART
!              SIGMAW(K)=WSTR*(((1.2*(1.0-0.9*ZH)*(ZH**0.6666))
!     &           +((1.8-1.4*ZH)*(USTR/WSTR)*(USTR/WSTR)))**0.5)
              IF(ZH.LT.0.03)THEN
                 SIGMAW(K)=WSTR*0.96*(3.0*ZH-SLEN/ZMIX)**0.3333
              ELSEIF(ZH.LT.0.4)THEN
                 SIGMAW(K)=WSTR*                                        &
     &                DMIN1(0.96*(3.0*ZH-SLEN/ZMIX)**0.3333,            &
     &                0.763*ZH**0.175)
              ELSEIF(ZH.LT.0.96)THEN
                 SIGMAW(K)=WSTR*0.722*(1-ZH)**0.207
              ELSE
                 SIGMAW(K)=WSTR*0.37
              END IF

! CHG(10/30/02): avoid SIGMAW becoming zero
              SIGMAW(K)=DMAX1(SIGMAW(K),DBLE(0.00001))

              IF((ZH.LT.0.1).AND.((ZZ(K)-Z0).LT.ABS(SLEN)))THEN
                TTL(K)=0.1*ZZ(K)/(SIGMAW(K)*(0.55+0.38*(ZZ(K)-Z0)/SLEN))
              ELSEIF((ZH.LT.0.1).AND.((ZZ(K)-Z0).GT.ABS(SLEN)))THEN
                TTL(K)=0.59*(ZZ(K)/SIGMAW(K))
              ELSE
                TTL(K)=0.15*(ZMIX/SIGMAW(K))*(1-EXP(-5.0*ZH))
              END IF

!           IF(SL.GE.0.0)THEN
            END IF

!          WRITE(45,*)"stbs K: ",K," SIGMAW ",SIGMAW(K)," ZH ",ZH,
!     :                "  TL ",TTL(K)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            IF(SL.GE.0.0)THEN
!              stable surface layer Beljaars-Holtslag
               ETRM=B*EXP(-D*SL)*(1.0-D*SL+C)
               PHIM=1.0+(A+ETRM)*SL
               PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*SL)+ETRM)*SL)

            ELSE
!              unstable surface layer Betchov-Yaglom / Kadar-Perepelkin
               PHIM=((1.0+0.625*SL*SL)/(1.0-7.5*SL))**0.3333
               PHIH=0.64*((3.0-2.5*SL)/(1.0-10.0*SL+50.0*SL*SL))**0.3333
            END IF

!           level within the surface layer
            IF(NINT(ZZ(K)).LE.NINT(ZSFC))THEN
               WH=USTR/PHIH
               VMIX=VONK*WH*ZZ(K)*(1.0-ZZ(K)/ZMIX)*(1.0-ZZ(K)/ZMIX)

!              save the last value of PHI to be valid at z/Zi=0.1
               PHIM0=PHIM
               PHIH0=PHIH

!           level within boundary layer
            ELSE
!              diabatic corrected Prandtl number
               PRD=(PHIH0/PHIM0)+(7.2*VONK*ZZ(K)*WSTR/ZMIX/WM)
               WH=WM/PRD
               VMIX=VONK*WH*ZZ(K)*(1.0-ZZ(K)/ZMIX)*(1.0-ZZ(K)/ZMIX)
            END IF

!           sum values for average (save last index for pbl)
            KPBL=K
            VSUM=VSUM+VMIX
            WGHT=WGHT+1.0

! JCL:(2/27/01)comment out because gives too strong entrainment
!=>compute mixing through the inversion layer (only for convective case)
!         ELSEIF(NINT(ZZ(K)).EQ.NINT(ZMIX).AND.WSTR.GT.0.0)THEN

!           Betts and Beljaars inversion layer jump model
!            DTDZ=(TT(K)-TT(K-1))/(ZZ(K)-ZZ(K-1))
!            DTDZ=DMAX1(DTDZ,DBLE(0.1))
!            VMIX=-0.4*TSTR*USTR/DTDZ

!=>level is within the free troposphere

         ELSE
!           bulk Richardson number
            DELZ=ZZ(K)-ZZ(K-1)
            DELU=(UU(K)-UU(K-1))**2+(VV(K)-VV(K-1))**2
            DELU=DMAX1(DELU,DBLE(0.1))
            RI=GRAV*DELZ*(TT(K)-TT(K-1))/DELU/TT(K-1)
            RI=DMAX1(DBLE(0.0),DMIN1(DBLE(20.0),RI))

!           vertical length scale (l)
            VL=1.0/(1.0/VONK/ZZ(K)+1.0/150.0)

!           stability scale (l/lo)
            IF(RI.GE.0.0.AND.RI.LE.0.001)THEN
               XL=1.0893*RI
            ELSE
               XL=A1+RI*(A2+RI*(A3+RI*(A4+A5*RI)))
            END IF

!           local stability function
            ETRM=B*EXP(-D*XL)*(1.0+C-D*XL)
            PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*XL)+ETRM)*XL)
            VMIX=VL*VL*ABS(SQRT(DELU)/DELZ)/PHIH
         END IF

!        check limits
         XX(K)=DMAX1(VMIN, DMIN1(VMAX, VMIX))

! JCL:   if particle in free troposphere or in inversion layer,
!           then use old parameterization, with constant TL & SIGMAW
!           determined from eddy diffusitivity
         IF(NINT(ZZ(K)).GT.NINT(ZMIX))THEN
            TTL(K)=100.0
            SIGMAW(K)=SQRT(XX(K)/TTL(K))
         END IF
      END DO
!     pbl mixing profile replaced by average value
!     up to the last level below the mixed layer depth
!     only average BL mixing if isotropic turbulence turned off
      IF(ISOT.EQ.0)THEN
         DO K=1,KPBL
            XX(K)=VSUM/WGHT
         END DO
      END IF

      RETURN
      END
