!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METINP           METeorological INPut reads meteo data file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL INPUT READS METEOROLOGICAL DATA FILE STARTING AT
!   AND FILLS DATA ARRAYS FOR ONLY ONE TIME PERIOD.  MISSING DATA ARE
!   SKIPPED AND OLD VALUES REMAIN IN THE ARRAY.  HOWEVER THIS MAY CAUSE
!   PROBLEMS IN OTHER ROUTINES WHEN DATA ARE REMAPPED TO NEW COORDINATES
!   OR UNIT CONVERSIONS. NO OTHER CHECKS ARE PERFORMED AND DATA FILL
!   ARRAY ACCORDING TO THE LEVEL NUMBER CONTAINED IN THE RECORD LABEL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Mar 1998 (RRD)
!                  16 Dec 1998 (RRD) - added end-of-file test KEND
!                  19 Apr 1999 (RRD) - added terrain height
!                  17 Jun 1999 (RRD) - correction to handle minutes
!                  02 Nov 1999 (RRD) - flag to test for data with sfc=0
!
! USAGE:  CALL METINP(BACK,KG,KUNIT,KREC,NXY,LX1,LY1,NXS,NYS,NZS,
!              CDATA,MC,KEND,IFHR,ZT,P0,T0,U0,V0,UF,VF,SF,RT,U,V,W,T,Q,P)
!   INPUT ARGUMENT LIST:
!     BACK      - log   flag to indicate integration direction
!     KG    - int number of active grid
!     KUNIT - int input device unit number
!     KREC  - int record number of index record
!     NXY   - int product of full grid dimensions
!     LX1,LY1     - int lower left corner of sub-grid in FG units
!     NXS,NYS     - int dimensions of sub-grid
!     NZS   - int number of data levels to read
!     CDATA - char      packed data array of length NXY
!     MC    - int accumulated minutes of data read
!     KEND  - int last valid record number of input file
!   OUTPUT ARGUMENT LIST:
!     IFHR  - int current forecast hour
!     ZT    - real      terrain height array (m)
!     P,T,U,V,W,Q - meteorological variables (see advpnt for complete list)
!   INPUT FILES:
!     UNIT 10,12,14,etc - according to how many input files defined
!   OUTPUT FILES:
!     UNIT 30 for diagnostic MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: metinp.f90,v 1.17 2008-08-12 09:54:44 skoerner Exp $
!
!$$$

 ! CHG:(11/20/01) add RC conv precip, TC total cld, SW radiation
 ! CHG:(11/20/01) add LF latent heat flux
 ! CHG:(12/04/01) add LC low cloud cover
 ! CHG:(01/22/03) add SM soil moisture
 ! JCL:(03/27/03) add arrays that get assigned if reading RAMS fields
 ! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
 ! CHG(09/25/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
 ! CHG(09/25/03) add RAMS turb. kin. energy TKEN
      SUBROUTINE METINP (BACK,KG,KUNIT,KREC,NXY,LX1,LY1,NXS,NYS,NZS,    &
         MC,KEND,IFHR,ZT,P0,T0,U0,V0,w0,hpbl,                           &
         muu, muv, mu, msfu, msfv, msft, fluxflg, deepflg, shallflg,    &
         UF,VF,SF,RT,U,V,W,T,Q,P,                                       &
         RC,LF,TC,LC,SW,SM,TLRAMS,SIGWRAMS,CFU1,CFU2,CFD1,DFU1,DFU2,    &
         EFU1,EFU2,DFD1,EFD1,TKEN)

      USE module_defgrid
      IMPLICIT REAL*8 (A-H,O-Z)

      ! missing times, variables
      INTEGER MISST, MISSV
      LOGICAL BACK, SFC1

      ! 3D meteorological variables
      REAL*8 U(NXS,NYS,NZM), V(NXS,NYS,NZM), W(NXS,NYS,NZM),          &
             T(NXS,NYS,NZM), Q(NXS,NYS,NZM), P(NXS,NYS,NZM)

      ! JCL:(03/27/03) add arrays that get assigned if reading RAMS fields
      REAL*8   TLRAMS(NXS,NYS,NZM), SIGWRAMS(NXS,NYS,NZM)

      ! CHG(09/23/03) RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
      REAL*8   CFU1(NXS,NYS,NZM),CFU2(NXS,NYS,NZM),CFD1(NXS,NYS,NZM),   &
               DFU1(NXS,NYS,NZM),DFU2(NXS,NYS,NZM),EFU1(NXS,NYS,NZM),   &
               EFU2(NXS,NYS,NZM),DFD1(NXS,NYS,NZM),EFD1(NXS,NYS,NZM)

      ! CHG(09/25/03) add RAMS turb. kin. energy TKEN
      REAL*8   TKEN(NXS,NYS,NZM)

      ! 2D meteorological variables
      ! CHG:(11/20/01) add RC conv precip, TC total cld, SW radiation
      ! CHG:(11/20/01) add LF latent heat flux
      ! CHG:(12/04/01) add LC low cloud cover
      ! CHG:(01/22/03) add SM soil moisture
      REAL*8 P0(NXS,NYS), U0(NXS,NYS), V0(NXS,NYS), T0(NXS,NYS), hpbl(NXS,NYS),      &
             UF(NXS,NYS), VF(NXS,NYS), SF(NXS,NYS), RT(NXS,NYS),      &
             ZT(NXS,NYS), RC(NXS,NYS), LF(NXS,NYS),                   &
             TC(NXS,NYS), LC(NXS,NYS), SW(NXS,NYS), SM(NXS,NYS), w0(nxs,nys)
      REAL*8 muu(nxs,nys),muv(nxs,nys),mu(nxs,nys), &
             msfu(nxs,nys),msfv(nxs,nys),msft(nxs,nys)
      LOGICAL, INTENT(inout) :: fluxflg, deepflg, shallflg
      CHARACTER*4 :: fluxvars(2) = (/'MUU0','MUV0'/), &
             wrfvars(6) = (/'MUBA','MUPE','ALT0','MSFU','MSFV','MSFT'/), &
             deepvars(6) = (/'CFU1','CFD1','DFU1','DFD1','EFU1','EFD1'/), &
             shallvars(3) = (/'CFU2','DFU2','EFU2'/)
      LOGICAL :: haveflux(2),havewrf(6),havedeep(6),haveshall(3),havehpbl
      ! data strings for unpacking
      CHARACTER LABEL*50, CDATA(NXY)*1, VARB*4, HEADER*108
      INTEGER :: CHKSUM=-1

      ! CHG(10/02/03) add ramsflag, need for 16 bit...
      LOGICAL :: RAMSFLG, ECMFLG, SXTNBIT, awrfflg, ppertflg

      REAL*8 ppert(NXS,NYS,NZM), mupert(nxs,nys)


!---------------------------------------------------------------------------------------------------
      ! set flag default to indicate that surface fields have the level
      ! index set to 0 (new style) rather than 1 (old style)
      SFC1 = .FALSE.
      ppertflg = .FALSE.

      ! see defgrid.inc for following definitions
      NXG = GRID(KG)%NX
      NYG = GRID(KG)%NY

      ! set starting record number at index record position
      JREC = KREC

      ! CHG(10/02/03) add ramsflag, need for 16 bit...
      RAMSFLG = GRID(KG)%MODEL_ID == 'RAMS'
      ECMFLG  = GRID(KG)%MODEL_ID(1:2) == 'EC'
      AWRFFLG = GRID(KG)%MODEL_ID(2:4) .EQ. 'WRF'
      SXTNBIT = RAMSFLG .OR. GRID(KG)%MODEL_ID(1:3) == 'ECX' .OR. GRID(KG)%MODEL_ID == 'DWRF'
      haveflux = .FALSE.
      havewrf = .FALSE.
      havedeep = .FALSE.
      haveshall = .FALSE.
      DREC(KG)%TKEF = .FALSE.
      havehpbl = .FALSE.

      ! check for end-of-file
      IF (JREC > KEND) GOTO 900

      ! loop to this position when missing data found
      MISST = 0
  200 CONTINUE
      MISSV = 0

      ! read index record (or first sfc field) for data time
      READ (KUNIT,REC=JREC,ERR=910) LABEL,HEADER
      READ (LABEL,100) IY,IM,ID,IH,IFHR
      CALL TM2MIN (IY,IM,ID,IH,0,MC)

      IF ((DREC(KG)%TYPE) == 2) THEN
         READ (HEADER(1:108),'(4X,I3,97X,I4) ')IFHR,LENH
         ! the number of index records per time period
         IF (SXTNBIT) THEN
            NNDX = LENH/(NXG*NYG*2)+1
         ELSE
            NNDX = LENH/(NXG*NYG)+1
         ENDIF
         ! new format data skip index record for subsequent reads
         JREC = JREC+NNDX
         ! decode extended portion for forecast hour, and minutes
         READ (HEADER,'(4X,I3,I2) ')IFHR,MNTS
         ! add minutes field to accumulated time (6/17/99)
         MC = MC+MNTS
      END IF

      ! adjust read for data offset (old style SH)
      JREC = JREC+DREC(KG)%OFFSET

      ! loop through only number of levels in sub-grid
      ppert = 0  !in case ppert is available at only some levels
      w0 = 0  !use zero w at sfc unless available from input file
      DO KK=0,NZS

         ! start with level 0 surface parameters
         NVAR = DREC(KG)%NUM_VARB(KK+1)

         DO NN=1,NVAR

            READ (KUNIT,REC=JREC,ERR=920) LABEL,CDATA
            READ (LABEL,100) IY,IM,ID,IH,IC,LL,IG,VARB,NEXP,PREC,VAR1

            ! test for missing data
            IF (IC < 0 .OR. VARB == 'NULL') THEN
               ! eta missing flux correction
               ! JCL: (5/30/00) for unknown reason, get compiler errmsg if
               ! leave the following lines in lower-case
               ! JCL:(5/7/01) implement Draxler's changes for generalized input of flux fields made in July 28, 2000
               ! IF((.NOT.DREC(KB)%SFLG).AND.
               IF ((.NOT.DREC(KG)%UFLX) .AND. (VARB == 'UMOF' .OR. VARB == 'VMOF')) THEN
                  GO TO 300
               ELSE
               ! if((.not.drec(kg).sflg).and.
               ! :           (varb.eq.'UMOF'.or.varb.eq.'VMOF')) then
               ! go to 300
               ! else
               MISSV = MISSV+1
               WRITE (30,*) 'WARNING metinp: missing data - ',VARB
               WRITE (30,*) 'Time : ',IY,IM,ID,IH,IC
               WRITE (30,*) 'Level: ',LL,'   Variable: ',NN
               GO TO 300
               END IF
            END IF

            ! All old format data files will always have the surface data
            ! at the beginning of each time period.  Prior to 1992 the surface
            ! index was set to #1 and after that the surface started at #0.
            IF (KK == 0 .AND. NN == 1 .AND. LL == 1) SFC1 = .TRUE.
            ! old style data that starts at level #1 adjust index by one
            IF (SFC1) LL = LL-1

            IF (LL == 0) THEN
               ! surface level fields
               IF (VARB == 'PRSS')                                       &
                  CALL PAKINP (SXTNBIT,P0, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'TMPS' .OR. VARB == 'T02M')                     &
                  CALL PAKINP (SXTNBIT,T0, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'U10M')                                       &
                  CALL PAKINP (SXTNBIT,U0, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'V10M')                                       &
                  CALL PAKINP (SXTNBIT,V0, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'HPBL') THEN
                  havehpbl = .TRUE.
                  CALL PAKINP (SXTNBIT,HPBL, CDATA,NXG,NYG,NXY,LX1,LY1,    &
       	             NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ENDIF
               IF (VARB == 'WWND' .OR. VARB == 'DZDT')                                       &
                  CALL PAKINP (SXTNBIT,W0, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! JCL:(5/7/01) implement Draxler's changes for generalized input of flux fields made in July 28, 2000
               ! friction velocity or momentum fluxes
               IF (DREC(KG)%USTR) THEN
                  IF (VARB == 'USTR')                                    &
                     CALL PAKINP (SXTNBIT,UF, CDATA,NXG,NYG,NXY,LX1,LY1, &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ELSE
                  IF (VARB == 'UMOF' .OR. VARB == 'EXCO')                  &
                     CALL PAKINP (SXTNBIT,UF, CDATA,NXG,NYG,NXY,LX1,LY1, &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                  IF (VARB == 'VMOF')                                    &
                     CALL PAKINP (SXTNBIT,VF, CDATA,NXG,NYG,NXY,LX1,LY1, &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               END IF
               ! friction temperature or sensible heat flux
               IF (DREC(KG)%TSTR) THEN
                  IF (VARB == 'TSTR')                                    &
                     CALL PAKINP (SXTNBIT,SF, CDATA,NXG,NYG,NXY,LX1,LY1, &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ELSE
                  IF (VARB == 'SHTF' .OR. VARB == 'HFLX') THEN
                     CALL PAKINP (SXTNBIT,SF, CDATA,NXG,NYG,NXY,LX1,LY1, &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                     !  correction for sign reversal in ECMWF/ETA/EDAS sensible heat flux
                     IF (GRID(KG)%MODEL_ID(1:2) == 'EC'  .OR.                 &
                         GRID(KG)%MODEL_ID      == ' ETA' .OR. GRID(KG)%MODEL_ID == 'EDAS') SF = -SF
                  ENDIF
               END IF

               ! latent heat flux
               IF (VARB == 'LHTF' .OR. VARB == 'WFLX')                     &
                  CALL PAKINP (SXTNBIT,LF, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                  NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               !  correction for sign reversal in ECMWF/ETA/EDAS latent heat flux
               IF (GRID(KG)%MODEL_ID(1:2) == 'EC'  .OR.                 &
                   GRID(KG)%MODEL_ID      == ' ETA' .OR. GRID(KG)%MODEL_ID == 'EDAS') LF = -LF

               ! terrain height and precipitation
               ! IF(VARB == 'SHGT')                                       &
               IF (VARB == 'SHGT' .OR. VARB == 'TERR')                     &
                  CALL PAKINP (SXTNBIT,ZT, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB(1:3) == 'TPP')                                   &
                  CALL PAKINP (SXTNBIT,RT, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! CHG:(11/20/01) conv. precip.
               IF (VARB(1:3) == 'CPP' .OR. VARB == 'CPRC')                                   &
                  CALL PAKINP (SXTNBIT,RC, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! CHG:(11/20/01) total cloud cover
               IF (VARB == 'TCLD')                                       &
                  CALL PAKINP (SXTNBIT,TC, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! CHG:(12/04/01) low cloud cover
               IF (VARB == 'LCLD')                                       &
                  CALL PAKINP (SXTNBIT,LC, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! CHG:(11/20/01) shortw. radiative flux
               IF (VARB == 'DSWF') THEN
                  CALL PAKINP (SXTNBIT,SW, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                  SW = MAX(SW,0D0)
               ENDIF
               ! CHG:(22/01/03) soil moisture
               IF (VARB == 'SOLW')                                       &
                  CALL PAKINP (SXTNBIT,SM, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)

               IF (awrfflg) THEN
                 ! check for additional variables needed for momentum fluxes
                 where (VARB .EQ. fluxvars) haveflux = .TRUE.
                 where (VARB .EQ. wrfvars) havewrf = .TRUE.
                 IF (VARB == 'MUU0') &
                  CALL PAKINP (SXTNBIT,muu, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MUV0') &
                  CALL PAKINP (SXTNBIT,muv, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MUBA') &
                  CALL PAKINP (SXTNBIT,mu, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MUPE') &
                  CALL PAKINP (SXTNBIT,mupert, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MSFU') &
                  CALL PAKINP (SXTNBIT,msfu, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MSFV') &
                  CALL PAKINP (SXTNBIT,msfv, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
                 IF (VARB == 'MSFT') &
                  CALL PAKINP (SXTNBIT,msft, CDATA,NXG,NYG,NXY,LX1,LY1,    &
                     NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              ENDIF
            ELSE
               ! upper levels
               IF (VARB == 'UWND')                                       &
                  CALL PAKINP (SXTNBIT,U(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'VWND')                                       &
                  CALL PAKINP (SXTNBIT,V(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'WWND' .OR. VARB == 'DZDT')                                       &
                  CALL PAKINP (SXTNBIT,W(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'TEMP' .OR. VARB .EQ. 'THET')                                       &
                  CALL PAKINP (SXTNBIT,T(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'SPHU' .OR. VARB == 'RELH')                     &
                  CALL PAKINP (SXTNBIT,Q(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               ! pressure level data requires heights of pressure surface
               ! height will be replaced by pressure in prfprs subroutine
               IF (VARB == 'HGTS')                                       &
                  CALL PAKINP (SXTNBIT,P(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB == 'PRES')                                       &
                  CALL PAKINP (SXTNBIT,P(:,:,LL), CDATA,NXG,NYG,NXY,     &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
               IF (VARB .EQ. 'PPRE') THEN
                  ppertflg = .TRUE.
                  CALL PAKINP (SXTNBIT,Ppert(:,:,LL), CDATA,NXG,NYG,NXY, &
                     LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              ENDIF
              ! JCL(03/27/03):new RAMS variables
              ! JCL(03/27/03):turbulence variables are directly calculated by RAMS, based on TKE
              IF (VARB == 'TLGR' .AND. GRID(KG)%MODEL_ID == 'RAMS')        &
                CALL PAKINP (SXTNBIT,TLRAMS(:,:,LL), CDATA,NXG,NYG,NXY,  &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'SIGW' .AND. GRID(KG)%MODEL_ID == 'RAMS')        &
                CALL PAKINP (SXTNBIT,SIGWRAMS(:,:,LL), CDATA,NXG,NYG,NXY,&
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)

              ! convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
              IF (VARB == 'CFU1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,CFU1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'CFU2' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,CFU2(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'CFD1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,CFD1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              ! IF(VARB.EQ.'CFD1'.AND.(RAMSFLG .OR. ECMFLG .or. awrfflg))
              ! :          WRITE(45,*) 'C: ITS THERE!',LL,LX1,LY1,CFD1(7,7,LL)
              IF (VARB == 'DFU1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,DFU1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'DFU2' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,DFU2(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'EFU1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,EFU1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'EFU2' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,EFU2(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'DFD1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,DFD1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              IF (VARB == 'EFD1' .AND. (RAMSFLG .OR. ECMFLG .OR. awrfflg))        &
                CALL PAKINP (SXTNBIT,EFD1(:,:,LL), CDATA,NXG,NYG,NXY,    &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)

              ! CHG(09/25/03) RAMS turb. kin. energy TKEN
              IF (VARB == 'TKEN' .AND. (ramsflg .OR. awrfflg)) THEN
                 DREC(KG)%TKEF = .TRUE.
                 CALL PAKINP (SXTNBIT,TKEN(:,:,LL), CDATA,NXG,NYG,NXY,    &
                              LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              END IF

               IF (awrfflg) THEN
                 ! check for additional variables needed for momentum and convective fluxes
                 where (VARB .EQ. fluxvars) haveflux = .TRUE.
                 where (VARB .EQ. wrfvars) havewrf = .TRUE.
                 where (VARB .EQ. deepvars) havedeep = .TRUE.
                 where (VARB .EQ. shallvars) haveshall = .TRUE.
                 IF (VARB == 'ALT0')        &
                   CALL PAKINP (SXTNBIT,TLRAMS(:,:,LL), CDATA,NXG,NYG,NXY,  &
                   LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,CHKSUM)
              ENDIF
            END IF

  300       JREC = JREC+1
         END DO

      END DO

      IF (ppertflg) THEN
         p = p + ppert
         pminp = 1e10
         DO ll=1,nzs
            ibad_pmin = 0
            pmin = minval(p(:,:,ll))
            IF (pmin < 0 .OR. pmin > pminp) THEN
               ibad_pmin = ll
               WRITE (45,*) 'metinp: bad pressure at time,ll,pminp,pmin : ',IY,IM,ID,IH,IC,ll,pminp,pmin
               WRITE (*,*) 'metinp: bad pressure at time,ll,pminp,pmin : ',IY,IM,ID,IH,IC,ll,pminp,pmin
            ENDIF
            pminp = pmin
         ENDDO
         IF (ibad_pmin .NE. 0) THEN
            DO ll=1,nzs
               write (45,'(a,i4,a,2g15.6)') 'll=',ll,' Min,max ptot=',minval(p(:,:,ll)),maxval(p(:,:,ll))
               IF (ppertflg) write (45,'(a,i4,a,2g15.6) ') &
                      'll=',ll,' Min,max ppert=',minval(ppert(:,:,ll)),maxval(ppert(:,:,ll))
            ENDDO
            stop 'metinp: bad pressure'
         ENDIF
      ENDIF
      IF (.NOT. havehpbl) THEN
         hpbl(:,:) = -99.
      END IF
      IF (awrfflg) THEN
         fluxflg = any(haveflux)
         IF (fluxflg .AND. .NOT. all(haveflux)) THEN
            WRITE (45,*) 'metinp: missing one or more fields for mom flux input:', &
                   ' haveflux= ',haveflux,' fluxvars= ',fluxvars
            WRITE (*,*) 'metinp: missing one or more fields for mom flux input:', &
                   ' haveflux= ',haveflux,' fluxvars= ',fluxvars
            stop 'metinp: missing one or more fields for mom flux input'
         ENDIF
         IF (.NOT. all(havewrf)) THEN
            WRITE (45,*) 'metinp: missing one or more fields for wrf input:', &
                   ' havewrf= ',havewrf,' wrfvars= ',wrfvars
            WRITE (*,*) 'metinp: missing one or more fields for wrf input:', &
                   ' havewrf= ',havewrf,' wrfvars= ',wrfvars
            stop 'metinp: missing one or more fields for wrf input'
         ENDIF
         deepflg = all(havedeep)
         shallflg = all(haveshall)
         IF (.NOT. DREC(KG)%TKEF) THEN
            tken(:,:,:) = 0.
         END IF
         IF (.NOT. deepflg) THEN
            cfu1(:,:,:) = 0.
            cfd1(:,:,:) = 0.
            dfu1(:,:,:) = 0.
            dfd1(:,:,:) = 0.
            efu1(:,:,:) = 0.
            efd1(:,:,:) = 0.
         ENDIF
         IF (.NOT. shallflg) THEN
            cfu2(:,:,:) = 0.
            dfu2(:,:,:) = 0.
            efu2(:,:,:) = 0.
         ENDIF
         mu = mu + mupert
      ENDIF !awrfflg

      ! JCL:(05/01/03) 15km MM5 terrain hts seem to be screwed up, so calculate from surface pressure
      ! JCL:(05/02/03) GFSx (longer forecast) terrain hts seem to be screwed up too
      ! JCL:(04/08/2004) 12km ETA terrain hts also seem to be screwed up
      ! JCL:(04/30/2004) RUC terrain hts may also seem to be screwed up
      IF ((GRID(KG)%MODEL_ID == ' MM5' .AND. IDINT(GRID(KG)%SIZE) == 15)   &
        .OR. (GRID(KG)%MODEL_ID == ' ETA' .AND. IDINT(GRID(KG)%SIZE) == 12)&
        .OR. (GRID(KG)%MODEL_ID == 'GFSN')                               &
        .OR. (GRID(KG)%MODEL_ID == ' RUC')) THEN
         HSCALE = 7700
         PSEA = 1025
         ZT(1:NXS,1:NYS) = HSCALE*DLOG(PSEA/P0(1:NXS,1:NYS))
      END IF

      ! up to five missing variables per time period permitted
      IF (MISSV > 10) THEN
         MISST = MISST+1
         ! only one missing time period permitted
         IF (MISST <= 1) THEN
            ! skip back to previous time period, otherwise continue
            IF (BACK) JREC = JREC-2*DREC(KG)%REC_PER
            GO TO 200
         ELSE
            WRITE (30,*) 'ERROR metinp: too many missing times - ',MISST
            STOP
         END IF
      END IF

      ! optional diagnostic information
      WRITE (30,'(A,3(I0,1X),I12,4X,4(1X,I2)) ') ' NOTICE metinp: ',KREC,LX1,LY1,MC,IY,IM,ID,IH
      RETURN

  900 WRITE(*,*) 'WARNING metinp: trying to read past end-of-file'
      RETURN
  910 WRITE(*,*) 'ERROR metinp: reading index record in meteo file'
      STOP
  920 WRITE(*,*) 'ERROR metinp: reading data record in meteo file'
      STOP

  100 FORMAT(7I2,A4,I4,2E14.7)

      END
