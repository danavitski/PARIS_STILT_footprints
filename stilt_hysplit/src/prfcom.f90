!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFCOM           PRoFile COMmon driver for input meteo
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE COMMON IS THE COMMON DRIVER FOR THE FOUR METEOROLGICAL
!   PROFILE ANALYSIS PROGRAMS.  ROUTINE EXAMINES THE PROFILE AT EACH
!   SUB-GRID NODE AND CONVERTS THE INPUT DATA TO COMMON UNITS AND
!   INTERPOLATES THE SOUNDING TO THE INTERNAL MODEL SIGMA SYSTEM.
!   THE DIFFERENT INPUT DATA COORDINATE SYSTEMS ARE SUPPORTED FOR
!   CONVERSION: ABSOLUTE PRESSURE (PRFPRS), SIGMA PRESSURE (PRFSIG),
!   TERRAIN FOLLOWING SIGMA (PRFTER), AND ECMWF HYBRID.  AFTER UNITS
!   CONVERSION THE STABILITY ANALYSIS ROUTINES ARE CALLED AT EACH POINT.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 01 Apr 1998 (RRD)
!                  18 Aug 1998 (RRD) - added isotropic turbulence variable
!                  04 Mar 1999 (RRD) - test for proper RH limits
!                  19 Apr 1999 (RRD) - added terrain array
!
! USAGE:  CALL PRFCOM(ISOT,KG,KSFC,GD,Z0,ZT,NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,
!              P0,T0,U0,V0,UF,VF,SF, U,V,W,T,Q,P,D,X)
!   INPUT ARGUMENT LIST:
!     ISOT  - int flag to indicate isotropic turbuence (0-no; 1-yes)
!     KG    - int grid indicator
!     KSFC  - int index that represents the top of the sfc layer
!     GD    - real      grid size (m)
!     Z0    - real      aerodynamic roughness length (m)
!     ZT    - real      terrain height elevations (m)
!     NXS,NYS     - int subgrid dimensions
!     NZS   - int number of input data levels
!     ZMDL  - real      internal model top (meters)
!     ZSG   - real      array internal sigma levels
!     NLVL  - int number of output levels
!     VMIX  - log indicator for mixing computation
!     P0,T0,etc - real  meteorological values at the surface (see advpnt)
!     U,V,T,Q,P,D,X - meteorological variables on input surface
!   OUTPUT ARGUMENT LIST:
!     U,V,T,Q,P,D,X - meteorological variables on model surface
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: prfcom.f90,v 1.13 2007-02-16 17:53:47 tnehrkor Exp $
!
!$$$

! JCL:add the array of Lagrangian timescale 'TL' &
!     stddev of vertical velocity 'SIGW' to argument list
! JCL:(4/28/00)add array of mixed-layer heights 'ZML' as input & output
! CHG:(12/04/01)add array of lim. of conv. heights 'ZLOC' as input & output
! JCL:(4/3/02)add mass violation array as argument
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer ht--used to prescribe mixed-layer ht
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
      SUBROUTINE PRFCOM(ISOT,KG,KSFC,GX,GY,Z0,ZT,NXS,NYS,NZS,ZMDL,ZSG,  &
     &   NLVL,VMIX,ICONVECT,P0,T0,U0,V0,W0,hpbl,UF,VF,SF,U,V,W,T,Q,P,D,X,TL,    &
     &   SIGW,ZML,ZLOC,DELMASS,ZISCALE,ALT0,mu,muu,muv,msfu,msfv,msft,fluxflg, deepflg, shallflg, &
     &      CFXUP1,       &
     &      CFXUP2,CFXDN1,DFXUP1,         &
     &      DFXUP2,EFXUP1,EFXUP2,         &
     &      DFXDN1,EFXDN1,TKEN )
!      SUBROUTINE PRFCOM(ISOT,KG,KSFC,GD,Z0,ZT,NXS,NYS,NZS,ZMDL,ZSG,
!     :   NLVL,VMIX,P0,T0,U0,V0,UF,VF,SF, U,V,W,T,Q,P,D,X)

      use module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorological array defined in this structure
!      INCLUDE 'DEFGRID.INC'

! JCL(03/27/03):add flag specifying whether data from RAMS or not
      LOGICAL RAMSFLG, AWRFFLG
!     3D meteorological variables
      REAL*8 U(NXS,NYS,NZM), V(NXS,NYS,NZM), W(NXS,NYS,NZM),            &
     &     T(NXS,NYS,NZM), Q(NXS,NYS,NZM), P(NXS,NYS,NZM),              &
     &     D(NXS,NYS,NZM), X(NXS,NYS,NZM), &
     &      CFXUP1(nxs,nys,nzm),       &
     &      CFXUP2(nxs,nys,nzm),CFXDN1(nxs,nys,nzm),DFXUP1(nxs,nys,nzm),         &
     &      DFXUP2(nxs,nys,nzm),EFXUP1(nxs,nys,nzm),EFXUP2(nxs,nys,nzm),         &
     &      DFXDN1(nxs,nys,nzm),EFXDN1(nxs,nys,nzm),TKEN(nxs,nys,nzm)

! CHG(09/18/03) add mass flux fields for DMASS calculation
      REAL*8 UMF(NXS,NYS,NZM), VMF(NXS,NYS,NZM), WMF(NXS,NYS,NZM)

! JCL:Lagrangian timescale & sigmaW
      REAL*8 TL(NXS,NYS,NZM),SIGW(NXS,NYS,NZM)
      REAL*8 ALT0(NXS,NYS,NZM)  !Dry inverse density for WRF
      REAL*8 mu(NXS,NYS),msft(NXS,NYS) !total mu and map scale factor at WRF mass grid points
      real*8 muu(NXS,NYS),muv(NXS,NYS),msfu(NXS,NYS),msfv(NXS,NYS) !mu and map-scale factors for u,v
      logical :: fluxflg !flag for WRF momentum flux input
      logical :: deepflg, shallflg

! JCL:(4/3/02)mass violation grid [fraction of gridcell/min]
      REAL*8 DELMASS(NXS,NYS,NZM)
      REAL*8 TESTMASS(NXS,NYS,NZM)
      REAL*8 TEST2(5,5)

!     2D meteorological variables
      REAL*8 P0(NXS,NYS), U0(NXS,NYS), V0(NXS,NYS), T0(NXS,NYS),        &
     &     UF(NXS,NYS), VF(NXS,NYS), SF(NXS,NYS), Z0(NXS,NYS),          &
     &     ZT(NXS,NYS), W0(NXS,NYS), hpbl(NXS,NYS)

      character (len=80) :: hpbl_msg=' ', &
           & hpbl_have='Using HPBL from file, replacing internally computed ZMIX', &
           & hpbl_miss='Missing HPBL in file, using internally computed ZMIX'

! JCL:(4/28/2000)array to store mixed-layer depth
      REAL*8 ZML(NXM,NYM)

! CHG:(12/04/2001)array to store lim. of conv. heights
      REAL*8 ZLOC(NXM,NYM)

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!     grid size array
!     REAL*8 GD(NXS,NYS)
      REAL*8 GX(NXS,NYS), GY(NXS,NYS)

!     input level data, output internal sigma levels
      REAL*8 PSG(MLVL), ZSG(NLVL), ESG(MLVL)

!     temporary holding arrays (note output arrays can be larger
!     than input data due to finer interpolation, although
!     they cannot exceed max dimension nzm

      REAL*8 U1(NZM),V1(NZM),W1(NZM),T1(NZM),Q1(NZM),P1(NZM)
      REAL*8 U2(NZM),V2(NZM),W2(NZM),T2(NZM),Q2(NZM),P2(NZM),           &
     &     Z2(NZM),D2(NZM),XV(NZM)
      real*8 alt1(NZM) !Dry inverse density for WRF
      real*8 CFXUP1_1(nzm),       &
     &      CFXUP2_1(nzm),CFXDN1_1(nzm),DFXUP1_1(nzm),         &
     &      DFXUP2_1(nzm),EFXUP1_1(nzm),EFXUP2_1(nzm),         &
     &      DFXDN1_1(nzm),EFXDN1_1(nzm),TKEN_1(nzm)
      real*8 CFXUP1_2(nzm),       &
     &      CFXUP2_2(nzm),CFXDN1_2(nzm),DFXUP1_2(nzm),         &
     &      DFXUP2_2(nzm),EFXUP1_2(nzm),EFXUP2_2(nzm),         &
     &      DFXDN1_2(nzm),EFXDN1_2(nzm),TKEN_2(nzm)


! JCL:temporary holding arrays for TL & SIGW
      REAL*8 TL2(NZM),SIGW2(NZM)

!     flag to indicate if vertical mixing information required
      LOGICAL VMIX

      save hpbl_msg, hpbl_have, hpbl_miss

! JCL:(4/5/02) parameters to convert P=>Density for mass violation claculation
!     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
      DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/

!      COMMON /GBLGRD/ GRID, DREC, FILE

!      WRITE(45,*)'In PRFCOM:'

! JCL(03/27/03): flag specifying whether data from RAMS or not
      RAMSFLG=GRID(KG)%MODEL_ID.EQ.'RAMS'

! JCL(050725): flag specifying whether data from WRFFLG or not
      AWRFFLG=GRID(KG)%MODEL_ID(2:4).EQ.'WRF'

!     remap vertical input data heights into input array
! CHG(09/08/03)replace loop
!      DO K=1,NZS
!         PSG(K)=DREC(KG)%HEIGHT(K+1)
!      END DO
       PSG(1:NZS)=DREC(KG)%HEIGHT(2:(NZS+1))

!     process each node on subgrid
      DO J=1,NYS
      DO I=1,NXS

!        place subgrid profile into temporary variables
         DO K=1,NZS
!           each level should have u,v,t
            U1(K)=U(I,J,K)
            V1(K)=V(I,J,K)
            T1(K)=T(I,J,K)
!           sigma coordinates - pressure is computed from hyspometric eq
!           pressure coordinates - contains height the replaced by pressure
            P1(K)=P(I,J,K)

!           vertical motion and humidity may not be present at all levels
            W1(K)=0.0
            Q1(K)=0.0
            if (awrfflg .and. (deepflg .or. shallflg)) tken_1(k) = tken(i,j,k)
            if (awrfflg .and. deepflg) then
               CFXUP1_1(k) = CFXUP1(i,j,k)
               CFXDN1_1(k) = CFXDN1(i,j,k)
               DFXUP1_1(k) = DFXUP1(i,j,k)
               EFXUP1_1(k) = EFXUP1(i,j,k)
               DFXDN1_1(k) = DFXDN1(i,j,k)
               EFXDN1_1(k) = EFXDN1(i,j,k)
            endif
            if (awrfflg .and. shallflg) then
               CFXUP2_1(k) = CFXUP2(i,j,k)
               DFXUP2_1(k) = DFXUP2(i,j,k)
               EFXUP2_1(k) = EFXUP2(i,j,k)
            endif
            IF(DREC(KG)%WFLG(K)) W1(K)=W(I,J,K)
            IF(DREC(KG)%RFLG(K)) Q1(K)=Q(I,J,K)
!           only needed for WRF:
            alt1(k)=0.0
            if (awrfflg) alt1(k)=alt0(i,j,k)
         END DO

!        remap data according to coordinate system
! TN (20050824): add wrf vertical coordinate system
         IF( AWRFFLG )THEN
!              input data on WRF pressure-sigma system
               CALL PRFwrf(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
     &            DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
     &            P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),w0(i,j),alt1, &
     &            mu(i,j),muu(i,j),muv(i,j),msfu(i,j),msfv(i,j),msft(i,j),fluxflg,&
     &            NZS,PSG,P1,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,              &
     &            P2,U2,V2,W2,T2,Z2,Q2,D2, &
     &            deepflg, shallflg, &
     &      CFXUP1_1,       &
     &      CFXUP2_1,CFXDN1_1,DFXUP1_1,         &
     &      DFXUP2_1,EFXUP1_1,EFXUP2_1,         &
     &      DFXDN1_1,EFXDN1_1,TKEN_1,           &
     &      CFXUP1_2,       &
     &      CFXUP2_2,CFXDN1_2,DFXUP1_2,         &
     &      DFXUP2_2,EFXUP1_2,EFXUP2_2,         &
     &      DFXDN1_2,EFXDN1_2,TKEN_2 )

         ELSEIF( (DREC(KG)%Z_FLAG).EQ.1)THEN

!              input data on pressure-sigma system
               CALL PRFSIG(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
     &            DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
     &            P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),              &
     &            NZS,PSG,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,                 &
     &            P2,U2,V2,W2,T2,Z2,Q2,D2)

         ELSEIF( (DREC(KG)%Z_FLAG).EQ.2)THEN

!              input data on pressure-absolute system
               CALL PRFPRS(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
     &            DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
     &            P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),              &
     &            NZS,PSG,P1,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,              &
     &            P2,U2,V2,W2,T2,Z2,Q2,D2)

         ELSEIF( (DREC(KG)%Z_FLAG).EQ.3)THEN

! JCL(03/27/03):add flag specifying whether data from RAMS or not
!              input data on terrain-sigma system
               CALL PRFTER(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
     &            DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
     &            P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),NZS,PSG,      &
     &            P1,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,                      &
     &            P2,U2,V2,W2,T2,Z2,Q2,D2,RAMSFLG)
!       WRITE(*,*)'prfcom',DREC(KG)%QFLG,DREC(KG)%UFLG,DREC(KG)%TFLG,
!     :   DREC(KG)%PRSS,DREC(KG)%SHGT,Q1

!      WRITE(45,*)'all levels:',Z2
!      WRITE(45,*)'p external:',i,j,P1
!      WRITE(45,*)'p internal:',i,j,P2

         ELSEIF( (DREC(KG)%Z_FLAG).EQ.4)THEN
!              input data on ecmwf hybrid sigma system
               CALL PRFECM(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
     &            DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
     &            P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),              &
     &            NZS,PSG,ESG,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,             &
     &            P2,U2,V2,W2,T2,Z2,Q2,D2)
         ELSE
               WRITE(*,*)'ERROR prfcom: no vertical data remapping'
               WRITE(*,*)'Unknown vertical type:', DREC(KG)%Z_FLAG
               STOP
         END IF


         IF (VMIX) THEN

! CHG:(11/19/01) pass on humidity (RH fraction) to subroutine STBANL as well
! CHG:(12/04/01) get lim. of convection height (ZLOCN) from stbanl
!           analyze surface stability (z/L) from fluxes  if available
!            CALL STBANL(KSFC,DREC(KG)%SFLG,DREC(KG)%EFLG,DREC(KG)%STAR,
!     :         Z0(I,J),UF(I,J),VF(I,J),SF(I,J),NLVL,U2,V2,T2,Z2,D2,
!     :         USTR,TSTR,WSTR,ZMIX,SLEN,PSI)
!           analyze surface stability (z/L) from fluxes  if available
! CHG:(9/17/02) add 'ICONVECT' as convection flag
            CALL STBANL(KSFC,DREC(KG)%EFLX,DREC(KG)%HFLX,DREC(KG)%UFLX, &
     &         DREC(KG)%USTR,DREC(KG)%TSTR,                             &
     &         Z0(I,J),UF(I,J),VF(I,J),SF(I,J),NLVL,U2,V2,T2,Z2,D2,Q2,  &
     &         USTR,TSTR,WSTR,ZMIX,SLEN,PSI,ZLOCN,ICONVECT)


            if (ziscale .lt. 0.0) then
! Special value to denote use of hpbl from met file, if available
               if (hpbl(i,j) .gt. -90.0) then
                  if (hpbl_msg .ne. hpbl_have) then
                     hpbl_msg = hpbl_have
                     write(30,'(1x,a)') hpbl_msg
                     write(45,'(1x,a)') hpbl_msg
                     write(*,'(1x,a)') hpbl_msg
                  endif
                  zmix = hpbl(i,j)
               else
                  if (hpbl_msg .ne. hpbl_miss) then
                     hpbl_msg = hpbl_miss
                     write(30,'(1x,a)') hpbl_msg
                     write(45,'(1x,a)') hpbl_msg
                     write(*,'(1x,a)') hpbl_msg
                  endif
               endif
            else
! JCL:(9/16/02) scaling factor used to prescribe mixed-layer height from file (ZICONTROL)
               ZMIX=ZMIX*ZISCALE
            end if
            
! JCL:(10/31/02) make sure that ZMIX is at NEAREST model level after applying scaling factor
            IF(ZISCALE.NE.1.0)THEN
            ZMIXNEW=Z2(2)
            DO KK=3,NLVL
               IF(DABS(Z2(KK)-ZMIX).LT.DABS(ZMIXNEW-ZMIX))THEN
                  ZMIXNEW=Z2(KK)
               END IF
            END DO
! JCL:(10/31/02) min mixed layer height set to second model level (~75 m)
            ZMIX=DMAX1(ZMIXNEW,Z2(2))
            END IF

! JCL:      add temporary holding arrays TL2&SIGW2 to hold results of
!             calculations for Lagrangian timescale & stddev of W
! JCL:      add the roughness length as argument as well
            CALL STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,             &
     &         NLVL,U2,V2,T2,Z2,XV,TL2,SIGW2,Z0(I,J))
!           compute vertical mixing profile (diffusivity)
!            CALL STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,
!     :         NLVL,U2,V2,T2,Z2,XV)

! JCL:(4/28/00)store mixed-layer height in array
            ZML(I,J)=ZMIX

! CHG:(12/04/01)store lim. of conv. height in array
            ZLOC(I,J)=ZLOCN

!           save friction velocity in u-momentum variable
            UF(I,J)=USTR
!           save integrated stability function in v-momentum variable
            VF(I,J)=PSI
         END IF

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
! JCL:   need following line to prevent GD(I,J)=0, causing model to crash
!        because would be dividing by 0!
!        IF(GD(I,J).EQ.0)GD(I,J)=GD(I,J-1)
         IF(GX(I,J).EQ.0)GX(I,J)=GX(I,J-1)
         IF(GY(I,J).EQ.0)GY(I,J)=GY(I,J-1)

!        wind speed conversion (m/s) to (grid/min)
! CHG(09/18/03) use interpolated grid size for RAMS
         IF(.NOT.RAMSFLG)THEN
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!          GSPD=60.0/GD(I,J)
           GSPD_X=60.0/GX(I,J)
           GSPD_Y=60.0/GY(I,J)
           if (awrfflg .and. fluxflg) then
! GX is valid at the mass point; GX*msfu/msft is the gx at the staggered u-point
              gspd_x = gspd_x * msft(i,j)/msfu(i,j)
              gspd_y = gspd_y * msft(i,j)/msfv(i,j)
           end if
         ELSE
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!          IF(I.LT.NXS)GSPDX=60.0*2.0/(GD(I+1,J)+GD(I,J)) !at u flux position i
!          IF(J.LT.NYS)GSPDY=60.0*2.0/(GD(I,J+1)+GD(I,J)) !at v flux position j
!          IF(I.EQ.NXS)GSPDX=60.0/GD(I,J)
!          IF(J.EQ.NYS)GSPDY=60.0/GD(I,J)
           !at u flux position i
           IF(I.LT.NXS)GSPDX=60.0*2.0/(GX(I+1,J)+GX(I,J))
           !at v flux position j
           IF(J.LT.NYS)GSPDY=60.0*2.0/(GY(I,J+1)+GY(I,J))
           IF(I.EQ.NXS)GSPDX=60.0/GX(I,J)
           IF(J.EQ.NYS)GSPDY=60.0/GY(I,J)
         END IF

         DO K=1,NLVL
           IF(.NOT.RAMSFLG)THEN
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!           U(I,J,K)=U2(K)*GSPD
!           V(I,J,K)=V2(K)*GSPD
            U(I,J,K)=U2(K)*GSPD_X
            V(I,J,K)=V2(K)*GSPD_Y
           ELSE
            UMF(I,J,K)=U2(K)
            VMF(I,J,K)=V2(K)
            WMF(I,J,K)=W2(K)
            U(I,J,K)=U2(K)*GSPDX
            V(I,J,K)=V2(K)*GSPDY
           END IF
            W(I,J,K)=W2(K)*60.0
!           temperature has been converted to potential
            T(I,J,K)=T2(K)
!           humidity is now relative humidity fraction
            Q(I,J,K)=DMAX1(DBLE(0.0),DMIN1(DBLE(1.0),Q2(K)))
!           local air density (kg/m3)
            D(I,J,K)=D2(K)
!           local pressure (mb)
            P(I,J,K)=P2(K)
!           vertical mixing coefficient (m2/s)
            X(I,J,K)=XV(K)
! JCL:      Lagrangian timescale (s)
            TL(I,J,K)=TL2(K)
! JCL:      stddev of vertical velocity (m/s)
            SIGW(I,J,K)=SIGW2(K)
            if (awrfflg .and. (deepflg .or. shallflg)) tken(i,j,k) = tken_2(k)
            if (awrfflg .and. deepflg) then
               CFXUP1(i,j,k) = CFXUP1_2(k)
               CFXDN1(i,j,k) = CFXDN1_2(k)
               DFXUP1(i,j,k) = DFXUP1_2(k)
               EFXUP1(i,j,k) = EFXUP1_2(k)
               DFXDN1(i,j,k) = DFXDN1_2(k)
               EFXDN1(i,j,k) = EFXDN1_2(k)
            endif
            if (awrfflg .and. shallflg) then
               CFXUP2(i,j,k) = CFXUP2_2(k)
               DFXUP2(i,j,k) = DFXUP2_2(k)
               EFXUP2(i,j,k) = EFXUP2_2(k)
            endif
         END DO

      END DO
      END DO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! JCL(4/3/02): following lines quantify MASS VIOLATION
!     NOTE:  still need to add the (Drho/Dt)s term--will be done in ADVPNT
!     loop over each grid position
!     mass violation array has one less element in each dimension,
!         except for the VERTICAL (b/c impose no-slip b.c. at the surface)
! CHG(09/18/03) use mass fluxes
      Model: IF (RAMSFLG) THEN


! CHG TEST scaling with dens(k)/dens0(k)
!      DO KK=1,NLVL
!        IF(KK.EQ.1)DENS0=1.1929761
!        IF(KK.EQ.2)DENS0=1.1794773
!        IF(KK.EQ.3)DENS0=1.1637981
!        IF(KK.EQ.4)DENS0=1.1456604
!        IF(KK.EQ.5)DENS0=1.1230835
!        IF(KK.GT.5)DENS0=1.0
!        IF(KK.EQ.1)DENS0U=1.1794773
!        IF(KK.EQ.2)DENS0U=1.1637981
!        IF(KK.EQ.3)DENS0U=1.1456604
!        IF(KK.EQ.4)DENS0U=1.1230835
!        IF(KK.EQ.5)DENS0U=1.1
!        IF(KK.GT.5)DENS0U=1.0
!          UMF(1:(NXS-1),1:NYS,KK)=UMF(1:(NXS-1),1:NYS,KK)*0.5*
!     &           (D(1:(NXS-1),1:NYS,KK)+D(2:NXS,1:NYS,KK))/DENS0
!          VMF(1:NXS,1:(NYS-1),KK)=VMF(1:NXS,1:(NYS-1),KK)*0.5*
!     &           (D(1:NXS,1:(NYS-1),KK)+D(1:NXS,2:NYS,KK))/DENS0
!          IF(KK.LT.NLVL)WMF(1:NXS,1:NYS,KK)=WMF(1:NXS,1:NYS,KK)*
!     &           (D(1:NXS,1:NYS,KK)/DENS0+D(1:NXS,1:NYS,KK+1)/DENS0U)
!          IF(KK.EQ.NLVL)WMF(1:NXS,1:NYS,KK)=WMF(1:NXS,1:NYS,KK)*
!     &           D(1:NXS,1:NYS,KK)/DENS0
!       WRITE(*,*)'prfcom:',dens0,dens0l,KK
!      ENDDO

      DO K=0,NLVL-1
      DO J=1,NYS-1
      DO I=1,NXS-1
! to get from mass flux density to mass flux (kg/s), need proper area for hor. fluxes
! need to interpolate terrain height to locations
         !at u-flux location I+1
         ZTT1=ZMDL-(ZT(I+1,J+1)+ZT(MIN0(I+2,NXS),J+1))/2
         !at u-flux location I
         ZTT2=ZMDL-(ZT(I,J+1)+ZT(I+1,J+1))/2
         !at v-flux location J+1
         ZTT3=ZMDL-(ZT(I+1,J+1)+ZT(I+1,MIN0(J+2,NYS)))/2
         !at v-flux location J
         ZTT4=ZMDL-(ZT(I+1,J)+ZT(I+1,J+1))/2
         ZTT=ZMDL-ZT(I+1,J+1)
! TEST
!        ZTT1=ZTT
!        ZTT2=ZTT
!        ZTT3=ZTT
!        ZTT4=ZTT

! CHG need to apply correct GD(i,j) ("delta-x" etc.) and correct area (through ZTT weighting)
!        HDIVU=(UMF(I+1,J+1,K+1)*ZTT1-UMF(I,J+1,K+1)*ZTT2)/
!     &        ZTT/GD(I+1,J+1)
!        HDIVV=(VMF(I+1,J+1,K+1)*ZTT3-VMF(I+1,J,K+1)*ZTT4)/
!     &        ZTT/GD(I+1,J+1)
! following works nicer..., but no idea what it means
!        HDIVU=(UMF(I+1,J+1,K+1)/ZTT1-UMF(I,J+1,K+1)/ZTT2)*
!     &        ZTT/GD(I+1,J+1)
!        HDIVV=(VMF(I+1,J+1,K+1)/ZTT3-VMF(I+1,J,K+1)/ZTT4)*
!     &        ZTT/GD(I+1,J+1)
!        HDIV=HDIVU+HDIVV
! horizontal divergence
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
! CHG need to apply correct GD(i,j) ("delta-x" etc.)
!       FMAP=45000.0/GD(I+1,J+1)
! WJS (4/1/05): Make grid size dynamic; multiply by 1000 to convert to meters
!        FMAP_X=45000.0/GX(I+1,J+1)
!        FMAP_Y=45000.0/GY(I+1,J+1)
        FMAP_X=(GRID(KG)%SIZE * 1000.0)/GX(I+1,J+1)
        FMAP_Y=(GRID(KG)%SIZE * 1000.0)/GY(I+1,J+1)
        HDIVU=(UMF(I+1,J+1,K+1)/ZTT-UMF(I,J+1,K+1)/ZTT)*ZMDL
        HDIVV=(VMF(I+1,J+1,K+1)/ZTT-VMF(I+1,J,K+1)/ZTT)*ZMDL
!       HDIVU=(UMF(I+1,J+1,K+1)*ZTT1-UMF(I,J+1,K+1)*ZTT2)*ZMDL/ZTT**2.0
!       HDIVV=(VMF(I+1,J+1,K+1)*ZTT3-VMF(I+1,J,K+1)*ZTT4)*ZMDL/ZTT**2.0
!       HDIVU=(UMF(I+1,J+1,K+1)/ZTT1-UMF(I,J+1,K+1)/ZTT2)*ZMDL
!       HDIVV=(VMF(I+1,J+1,K+1)/ZTT3-VMF(I+1,J,K+1)/ZTT4)*ZMDL
!       HDIVU=HDIVU/GD(I+1,J+1)
!       HDIVV=HDIVV/GD(I+1,J+1)

! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!       HDIVU=HDIVU*FMAP*FMAP/45000.0
!       HDIVV=HDIVV*FMAP*FMAP/45000.0
! WJS (4/1/05): Make grid size dynamic; multiply by 1000 to convert to meters
!        HDIVU=HDIVU*FMAP_X*FMAP_X/45000.0
!        HDIVV=HDIVV*FMAP_Y*FMAP_Y/45000.0
        HDIVU=HDIVU*FMAP_X*FMAP_X/(GRID(KG)%SIZE * 1000.0)
        HDIVV=HDIVV*FMAP_Y*FMAP_Y/(GRID(KG)%SIZE * 1000.0)

        HDIV=HDIVU+HDIVV
! CHG use different map factors, specific for each flux location
! ->only minimal changes, not significant...
!        FMAP=45000.0/GD(I+1,J+1)
!        FMAPU1=2.0*45000.0/(GD(MIN0(I+2,NXS),J+1)+GD(I+1,J+1))
!        FMAPU2=2.0*45000.0/(GD(I+1,J+1)+GD(I+1,J+1))
!        FMAPV1=2.0*45000.0/(GD(I+1,MIN0(J+2,NYS))+GD(I+1,J+1))
!        FMAPV2=2.0*45000.0/(GD(I+1,J+1)+GD(I+1,J))
!        HDIVU=UMF(I+1,J+1,K+1)*FMAPU1**2.0*ZMDL/ZTT1
!     :       -UMF(I,J+1,K+1)*FMAPU2**2.0*ZMDL/ZTT2
!        HDIVV=VMF(I+1,J+1,K+1)*FMAPV1**2.0*ZMDL/ZTT3
!     :       -VMF(I+1,J,K+1)*FMAPV2**2.0*ZMDL/ZTT4
!        HDIVU=HDIVU/45000.0
!        HDIVV=HDIVV/45000.0
!        HDIV=HDIVU+HDIVV
! CHG also need vertical map factor
!        HDIVU=(UMF(I+1,J+1,K+1)/ZTT1-UMF(I,J+1,K+1)/ZTT2)*
!     &        ZMDL/GD(I+1,J+1)
!        HDIVV=(VMF(I+1,J+1,K+1)/ZTT3-VMF(I+1,J,K+1)/ZTT4)*
!     &        ZMDL/GD(I+1,J+1)
!        HDIV=HDIVU+HDIVV
! vertical divergence
! back to w* fluxes
                                     !w* fluxes
        WMFS1=-WMF(I+1,J+1,K+1)*ZMDL
        IF(K.GT.0)WMFS2=-WMF(I+1,J+1,K)*ZMDL
        IF(K.EQ.0)WMFS2=0.0
                                                 !z* altitude difference
        IF(K.GT.0)DZ=-1.0*(ZSG(K+1)-ZSG(K))*ZMDL
        IF(K.EQ.0)DZ=-1.0*(ZSG(K+1)-1.0)*ZMDL
        VDIV=(WMFS1-WMFS2)/DZ
!        VDIV=VDIV*FMAP*FMAP
!     &      * (ZMDL-ZT(I+1,J+1))/ZMDL

!        IF(K.LT.NLVL-1)THEN
!          VDIV=(WMF(I,J,K+2)-WMF(I,J,K+1))/(ZSG(K+2)-ZSG(K+1))
!        ELSE !near top
!          VDIV=(WMF(I,J,K+1)-WMF(I,J,K))/(ZSG(K+1)-ZSG(K))
!        END IF
!        WRITE(45,*)'prfcom:',I,J,K,HDIVU,HDIVV,VDIV
!        WRITE(45,*)I,ZTT,HDIVU,HDIVV,VDIV,D(I+1,J+1,K+1)
!        WRITE(45,*)I,UMF(I+1,J+1,K+1),HDIVU,HDIVV,VDIV,D(I+1,J+1,K+1)
!     &      ,FMAPU1,FMAPU2,FMAPV1,FMAPV2,FMAP
!     &      ,ZTT1,ZTT
! CHG&JCL (10/07/03) normalize by density in advpnt
!        DELMASS(I+1,J+1,K+1)=(HDIV+VDIV)/D(I+1,J+1,K+1)*60.0 !also convert to 1/minute
                                              !also convert to 1/minute
        DELMASS(I+1,J+1,K+1)=(HDIV+VDIV)*60.0
!        WRITE(45,*)I,J,K,HDIVU,HDIVV,VDIV,DELMASS(I+1,J+1,K+1)
             !of DO I=1,NXS-1
      END DO
             !of DO J=1,NYS-1
      END DO
             !of DO K=0,NLVL-1
      END DO
           ! for awrf:
      ELSE IF (awrfflg) THEN Model
! C-grid setup in WRF: for x-staggered u and unstaggered m grid points:
!                        i-1  i   i  i+1 i+1
!                         m   u   m   u   m
!
!                      for z-staggered w and unstaggered m (and u,v) grid points:
!                        k-1  k   k  k+1 k+1 k+2
!                         m   w   m   w   m   w
      DO K=1,NLVL
      DO J=1,NYS
      DO I=1,NXS
! couple decoupled velocities at staggered u,v-grid points with density
        UMF(i,j,k)=0.5*(D(max(i-1,1),j,k)+D(i,j,k))*u(i,j,k)
        VMF(i,j,k)=0.5*(D(i,max(j-1,1),k)+D(i,j,k))*v(i,j,k)
! couple decoupled w at staggered w-grid points with density
        WMF(i,j,k)=0.5*(D(i,j,max(1,k-1))+D(i,j,k))*w(i,j,k)
      enddo
      enddo
      enddo
      DO K=0,NLVL-1
! Compute depth of layer bounded by w-levels:
      if (k .le. 0) then
! kk=k+1=1: layer bounded by zsigw(1) and zsigw(2)
         zsigw_dn=1.
         zsigw_up=2*zsg(k+1) - 1. !zs of second w-level
      else
! kk=k+1: layer bounded by zsigw(kk) and zsigw(kk+1)
         zsigw_up=2.*zsg(k+1)-zsigw_dn  !zs of w-level kk+1
      endif
      DZSIG=zsigw_up-zsigw_dn
      zsigw_dn=zsigw_up
      DO J=1,NYS-1
      DO I=1,NXS-1
! compute zmdl-terrain height
         !at u-flux location I+1
         ZTT1=ZMDL-(ZT(I,J)+ZT(MIN(I+1,NXS),J))/2
         !at u-flux location I
         ZTT2=ZMDL-(ZT(max(I-1,1),J)+ZT(I,J))/2
         !at v-flux location J+1
         ZTT3=ZMDL-(ZT(I,J)+ZT(I,MIN(J+1,NYS)))/2
         !at v-flux location J
         ZTT4=ZMDL-(ZT(I,max(J-1,1))+ZT(I,J))/2
         ZTT=ZMDL-ZT(I,J)
!        hor divergence term (not have to divide by GD b/c already done so in wind speed conversion above)
!        (computed at mass level,
         DUDXUP=UMF(I+1,J,K+1)*ztt1 - UMF(I,J,K+1)*ztt2
         DVDYUP=VMF(I,J+1,K+1)*ztt3 - VMF(I,J,K+1)*ztt4
         rhobar=d(i,j,k+1)
         if (k+1 .ge. NLVL) then
            wup = 0         !upper b.c. for w
         else
            wup=wmf(i,j,k+2)
         endif
         wdn=wmf(i,j,k+1)
! same as for other non-RAMS models (see below)         
         HDIV=(DUDXUP+DVDYUP)/ZTT
         VDIV=(WUP-WDN)/DZSIG
         DELMASS(I,J,K+1)=(HDIV+VDIV)/RHOBAR
      enddo
      enddo
      enddo
           ! for other than RAMS or awrf
      ELSE Model
        UMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
     &                          U(1:NXS,1:NYS,1:NLVL)
        VMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
     &                          V(1:NXS,1:NYS,1:NLVL)
        WMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
     &                          W(1:NXS,1:NYS,1:NLVL)
      DO K=0,NLVL-1
      DO J=1,NYS-1
      DO I=1,NXS-1
         ZTT1=ZMDL-ZT(I,J)
         ZTT2=ZMDL-ZT(I+1,J)
         ZTT3=ZMDL-ZT(I,J+1)
         ZTT4=ZMDL-ZT(I+1,J+1)
         ZTT=(ZTT1+ZTT2+ZTT3+ZTT4)/4.0
!        hor divergence term(not have to divide by GD b/c already done so in wind speed conversion above
         DUDXUP=(UMF(I+1,J+1,K+1)*ZTT4+                                 &
     &           UMF(I+1,J,K+1)*ZTT2)/2.0-                              &
     &          (UMF(I,J+1,K+1)*ZTT3+                                   &
     &           UMF(I,J,K+1)*ZTT1)/2.0
         DVDYUP=(VMF(I,J+1,K+1)*ZTT3+                                   &
     &           VMF(I+1,J+1,K+1)*ZTT4)/2.0-                            &
     &          (VMF(I+1,J,K+1)*ZTT2+                                   &
     &           VMF(I,J,K+1)*ZTT1)/2.0
         WUP=(WMF(I,J,K+1)+WMF(I+1,J,K+1)                               &
     &       +WMF(I,J+1,K+1)+WMF(I+1,J+1,K+1))/4.0

!        NOT near ground
         IF(K.GT.0)THEN
!           horizontal divergence in lower model layer
            DUDXDN=(UMF(I+1,J+1,K)*ZTT4+                                &
     &           UMF(I+1,J,K)*ZTT2)/2.0-                                &
     &          (UMF(I,J+1,K)*ZTT3+                                     &
     &           UMF(I,J,K)*ZTT1)/2.0
            DVDYDN=(VMF(I,J+1,K)*ZTT3+                                  &
     &           VMF(I+1,J+1,K)*ZTT4)/2.0-                              &
     &          (VMF(I+1,J,K)*ZTT2+                                     &
     &           VMF(I,J,K)*ZTT1)/2.0
!           vertical velocity at lower model layer
            WDN=(WMF(I,J,K)+WMF(I+1,J,K)                                &
     &       +WMF(I,J+1,K)+WMF(I+1,J+1,K))/4.0
            RHODN=(D(I,J,K)+D(I+1,J,K)+D(I,J+1,K)+D(I+1,J+1,K))/4.0
            RHOUP=(D(I,J,K+1)+D(I+1,J,K+1)+                             &
     &             D(I,J+1,K+1)+D(I+1,J+1,K+1))/4.0
            RHOBAR=(RHODN+RHOUP)/2.0
            DZSIG=ZSG(K+1)-ZSG(K)
!        near ground
         ELSE
!           no-slip b.c. at ground surface
            DUDXDN=0.0
            DVDYDN=0.0
            WDN=0.0
            RHOUP=(D(I,J,K+1)+D(I+1,J,K+1)+                             &
     &             D(I,J+1,K+1)+D(I+1,J+1,K+1))/4.0
! JCL:(5/2/02) NGM data doesn't have surface temp data, so simply assign RHOUP to RHODN
            IF(GRID(KG)%MODEL_ID.EQ.'NGM')THEN
               RHODN=RHOUP
            ELSE
!              density at the surface--from Ideal Gas Law
               RHODN=(P2JM/RDRY)*(P0(I,J)/T0(I,J)+P0(I+1,J)/T0(I+1,J)+  &
     &            P0(I,J+1)/T0(I,J+1)+P0(I+1,J+1)/T0(I+1,J+1))/4.0
            END IF
            RHOBAR=(RHODN+RHOUP)/2.0
            DZSIG=ZSG(K+1)-1.0
         END IF

!        horizontal divergence term--average over upper and lower model layers
! fixed bug: was (DVDYDN+DVDYDN)
         HDIV=((DUDXDN+DUDXUP)/2.0+(DVDYDN+DVDYUP)/2.0)/ZTT
!        vertical divergence term
         VDIV=(WUP-WDN)/DZSIG

!        assign element to mass violation array; still need to add the (Drho/Dt)s term
!        note convention: mass viol betw levels 1 & 2, e.g., gets stored in level 2 (K+1) of mass viol grid
!        mass violation is in units of [fraction of mass in gridcell/min]
         DELMASS(I,J,K+1)=(HDIV+VDIV)/RHOBAR

! JCL:   see what happens when adjust vertical wind to conserve mass
!         TMPDELMASS=(HDIV+VDIV)/RHOBAR
!         W(I+1,J+1,K+1)=W(I+1,J+1,K+1)-4*(RHOBAR*(DZSIG)/
!     &                  D(I+1,J+1,K+1))*DELMASS(I,J,K+1)
!        recalculate mass violation after adjusted vertical wind
!         WUP=(W(I,J,K+1)*D(I,J,K+1)+W(I+1,J,K+1)*D(I+1,J,K+1)
!     &      +W(I,J+1,K+1)*D(I,J+1,K+1)+W(I+1,J+1,K+1)*D(I+1,J+1,K+1))/4.0
!         VDIV=(WUP-WDN)/DZSIG
!         DELMASS(I,J,K+1)=(HDIV+VDIV)/RHOBAR
!         WRITE(45,*)'DMASS:',I,J,K,TMPDELMASS,DELMASS(I,J,K+1)
              !of DO I=1,NXS-1
       END DO
              !of DO J=1,NYS-1
       END DO
              !of DO K=0,NLVL-1
       END DO
             !of IF(RAMSFLAG) ... ELSE ...
      END IF Model


!       WRITE(45,*)'In PRFCOM:'
!       TESTMASS(1:(NXS-1),1:NYS,1:10)=(U(2:NXS,1:NYS,1:10)-
!     &          U(1:(NXS-1),1:NYS,1:10))/4.0
!       TESTMASS(1:(NXS-1),1:NYS,11:12)=0.0
!       DO K=1,14
!       DO J=1,NYS
!       DO I=1,NXS-1
!         WRITE(45,*)I,J,K,(U(I+1,J,K)-U(I,J,K))/4.0,TESTMASS(I,J,K)
!       END DO
!       END DO
!       END DO
!       WRITE(45,*)'DU:',TESTMASS(1:6,1:6,1)
!       WRITE(45,*)'DU2:',TEST2(1:5,1:5)
!       WRITE(45,*)'UGRID:',U(1:6,1:6,1)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      END
