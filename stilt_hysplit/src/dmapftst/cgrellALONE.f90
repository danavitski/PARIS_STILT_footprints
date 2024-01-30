!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CGRELL           Convective redistribution using GRELL scheme
!   PRGMMR:    CHRISTOPH GERBIG   DATE:03-13-03
!
! ABSTRACT:
!   Convective redistribution of particles using updraft & downdraft mass
!   fluxes, entrainment and detrainment fluxes. The vertical advection
!   is done on a fixed vertical grid, with dynamic internal timestep.
!   First, particles are assigned to updraft, downdraft or environment,
!   then they are advected to the next grid location. Then they are advected
!   using the fixed grid, until time is used up. Finally they are advected
!   from there position before time was up to the exact (off-grid) location
!   corresponding to the exact time at the end of the convective
!   redistribution.
!   Initially particles are assigned a cloud index (cloud updraft, downdraft,
!   or environment. When the area coverage changes, there is a probability
!   that they change there cloud index.
!
!
! PROGRAM HISTORY LOG:
!   LAST REVISED:
!
! USAGE:  CALL CGRELL(CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,EFXUP1,EFXUP2,
!              DFXDN1,EFXDN1,RAUP1,RAUP2,RADN1,AREA,AREAPRU,AREAPRD,DENS,NZM,
!              JET,ZPROFM,Z1,Z2,ICNDX1,ICNDX2,ZX,DT,BACK)
!   INPUT ARGUMENT LIST:
!     ?FX*  - real      Cloud/Entrainment/Detrainment Mass Fluxes
!                 (Updraft, downdraft; 1:deep; 2: shallow)
!                 UNITS: kg/m^2s; convention: upward=positive
!                       TO GET MASSFLUX (kg/s):
!                       MASSFLUX=?FX* * DX*DY (up/dndraft & hor. fluxes)
!                       numbers represent ro*w*a w/ a=fractional coverage
!                       to get mass flux[kg/s]: ?FX* * area (gridcell)
!     RA*   - real      updraft area coverage (fract.),
!                 (Updraft, downdraft; 1:deep; 2: shallow)
!     AREA  - real      gridcell area (m^2)
!     AREAPRU     - real      area coverage updr. (fract.) of particle during last call
!     AREAPRD     - real      area coverage downdr. (fract.) of particle during last call
!     DENS  - real      Density profile [kg/m^3]
!     NZM   - int vertical grid dimension
!     JET   - int current elapsed time (minutes)
!     ZPROFM      - real      profile of level altitudes above ground, FLX
!                 START: First level above ground (surface not incl.)
!     Z1    - real      old position at time t
!     ICNDX1      - int cloud index (whether in up/dndrft/environment)
!                 at time t
!     DT    - real      integration step (minutes)
!     BACK      - log   logical flag to indicate direction
!   OUTPUT ARGUMENT LIST:
!     Z2    - real      new position at time t+dt
!     ICNDX2      - real      cloud index (whether in up/dndrft/environment)
!                 at time t+dt
!     AREAPRU     - real      area coverage updr. (fract.) of particle after convection
!     AREAPRD     - real      area coverage downdr. (fract.) of particle after convection
!     ZX    - int last estimate of vertical index
!     JET   - int current elapsed time (minutes)
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  LINUX
!
!$$$

       SUBROUTINE CGRELL(CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,EFXUP1,     &
     &                   EFXUP2,DFXDN1,EFXDN1,RAUP1,RAUP2,RADN1,AREA,   &
     &                   AREAPRU,AREAPRD,DENS,NZM,JET,ZPROFM,Z1,Z2,     &
     &                   ICNDX1,ICNDX2,ZX,DT,BACK,RSEED)

      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BACK
      REAL*8   CFXUP1(NZM),CFXUP2(NZM),CFXDN1(NZM)
      REAL*8   DFXUP1(NZM),DFXUP2(NZM),DFXDN1(NZM)
      REAL*8   EFXUP1(NZM),EFXUP2(NZM),EFXDN1(NZM)
      REAL*8   RAUP1(NZM),RAUP2(NZM),RADN1(NZM)
      REAL*8   DENS(NZM),ZPROFM(NZM),ZPROFT(NZM)
      REAL*8   DZT(NZM)

      INTEGER ICNDX1,ICNDX2,RSEED,ZX

!    grid size arrays for Fluxes
      REAL*8 FUP(3,NZM),FLE(3,NZM),FRI(3,NZM),FDN(3,NZM),FT(3,NZM)
      REAL*8 FDN1(3,NZM),FI(3,NZM),FBUD(3,NZM),FUPS(3,NZM),FDNS(3,NZM)
      REAL*8 FRIS(3,NZM),FLES(3,NZM),FTS(3,NZM)

!    grid size arrays for probabilities
      REAL*8 PUP(3,NZM),PLE(3,NZM),PRI(3,NZM),PDN(3,NZM),PTO(3,NZM)

!    grid size arrays for w's, density and areas
      REAL*8 DENSM(3,NZM),AREAA(3,NZM)
      REAL*8 WUP(3,NZM),WDN(3,NZM),WTT(3,NZM)

!    grid size arrays for MASS and DZ's
      REAL*8 DZTA(3,NZM)

      DATA PI/3.1415926/
      DATA ONE/1.0/

!     get DZT as grid spacing
      ZPROFT(2:NZM)=0.5*(ZPROFM(1:(NZM-1))+ZPROFM(2:NZM))
      ZPROFT(1)=0.5*ZPROFM(1)
      DZT(1:NZM-1)=ZPROFT(2:NZM)-ZPROFT(1:(NZM-1))
      DZT(NZM)=ZPROFT(NZM)-ZPROFT(NZM-1)
      DZTA(1,1:NZM)=DZT
      DZTA(2,1:NZM)=DZT
      DZTA(3,1:NZM)=DZT
!     get areas in array
!     relative area coverage is defined at flux levels
      AREAA(1,1:NZM)=RAUP1+RAUP2
      AREAA(3,1:NZM)=RADN1
      AREAA(2,1:NZM)=1-AREAA(1,1:NZM)-AREAA(3,1:NZM)
!     from relative area to m^2
      AREAA(1:3,1:NZM)=AREAA(1:3,1:NZM)*AREA

!     put upward fluxes in 2D array
!     take into account time direction
!     ALL FLUXES are in MASSFLUX/GRIDAREA (=dxdy)
!     column 1 is updraft, 2 is environment, 3 is downdraft
      IF(.NOT.BACK)THEN
        FUP(1,1:NZM)=DABS(CFXUP1)+DABS(CFXUP2)
        FUP(3,1:NZM)=1.0E-15
!     put downward fluxes in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FDN(1,1:NZM)=1.0E-15
        FDN(3,1:NZM)=-DABS(CFXDN1)
!     put fluxes to "right" in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FRI(1,1:NZM)=DABS(DFXUP1(1:NZM))+DABS(DFXUP2(1:NZM))
        FRI(2,1:NZM)=DABS(EFXDN1(1:NZM))
        FRI(3,1:NZM)=1.0E-15
!     put fluxes to "left" in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FLE(2,1:NZM)=DABS(EFXUP1(1:NZM))+DABS(EFXUP2(1:NZM))
        FLE(3,1:NZM)=DABS(DFXDN1(1:NZM))
        FLE(1,1:NZM)=1.0E-15
!     reverse flows for backward
      ELSEIF(BACK)THEN
        FDN(1,1:NZM)=-DABS(CFXUP1)-DABS(CFXUP2)
        FDN(3,1:NZM)=1.0E-15
!     put downward fluxes in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FUP(1,1:NZM)=1.0E-15
        FUP(3,1:NZM)=DABS(CFXDN1)
        FLE(2,1:NZM)=DABS(DFXUP1(1:NZM))+DABS(DFXUP2(1:NZM))
        FLE(3,1:NZM)=DABS(EFXDN1(1:NZM))
        FLE(1,1:NZM)=1.0E-15
        FRI(1,1:NZM)=DABS(EFXUP1(1:NZM))+DABS(EFXUP2(1:NZM))
        FRI(2,1:NZM)=DABS(DFXDN1(1:NZM))
        FRI(3,1:NZM)=1.0E-15
      END IF
!     get environment fluxes from
!     mass balance; what goes up must come down
      FUP(2,1:NZM)=-FUP(1,1:NZM)-FUP(3,1:NZM)                           &
     &             -FDN(1,1:NZM)-FDN(3,1:NZM)
!     from net flux, seperate up and downward fluxes
!     no heavyside function, so use (sign(x)+1)/2
      FDN(2,1:NZM)=FUP(2,1:NZM)*(DSIGN(ONE,-FUP(2,1:NZM))+1.0)*0.5
      FUP(2,1:NZM)=FUP(2,1:NZM)*(DSIGN(ONE,FUP(2,1:NZM))+1.0)*0.5
!     set small values to 1.0E-15
      FUPS=(-DSIGN(ONE,DABS(FUP)-1.0E-8)+1.0)*0.5
      FUP=1.0E-15*FUPS+FUP*(ONE-FUPS)
      FDNS=(-DSIGN(ONE,DABS(FDN)-1.0E-8)+1.0)*0.5
      FDN=-1.0E-15*FDNS+FDN*(ONE-FDNS)

!     put density in array:
!     DENSM at flux levels (mean of T neighbours)
      DENSM(1,1:NZM-1)=0.5*(DENS(1:NZM-1)+DENS(2:NZM))
      DENSM(1,NZM)=DENS(NZM)
      DENSM(2,1:NZM)=DENSM(1,1:NZM)
      DENSM(3,1:NZM)=DENSM(1,1:NZM)

!     get w (w=F/ro) (WTT from T- to T-level)
!     set small values to 1.0E-15
      WUP=FUP*AREA/AREAA/DENSM
      WUP=1.0E-15*FUPS+WUP*(ONE-FUPS)
      WDN=FDN*AREA/AREAA/DENSM
      WDN=-1.0E-15*FDNS+WDN*(ONE-FDNS)
      WTT=WUP+WDN

!     get probabilities to move to above, below, right and left
!     shift flux down so that FDN1 is flux leaving current level
!      FDN1=FDN
      FDN1(1:3,2:NZM)=FDN(1:3,1:NZM-1)
      FDN1(1:3,1)=1.0E-15
!     set small values to 0.0
!     NEED TO USE CORRECT PROBABILITIES WHEN FLUXES ARE ZERO
      FDNS=(-DSIGN(ONE,DABS(FDN1)-1.0E-8)+1.0)*0.5
      FRIS=(-DSIGN(ONE,DABS(FRI)-1.0E-8)+1.0)*0.5
      FLES=(-DSIGN(ONE,DABS(FLE)-1.0E-8)+1.0)*0.5
      FT=DABS(FUP)*(ONE-FUPS)+DABS(FDN1)*(ONE-FDNS)                     &
     &  +DABS(FRI)*(ONE-FRIS)+DABS(FLE)*(ONE-FLES)
      FTS=(-DSIGN(ONE,DABS(FT)-1.0E-8)+1.0)*0.5
      FT=FT*(ONE-FTS)+1.0E-15*FTS
      PUP=DABS(FUP)*(ONE-FUPS)/FT
      PDN=DABS(FDN1)*(ONE-FDNS)/FT
      PRI=DABS(FRI)*(ONE-FRIS)/FT
      PLE=DABS(FLE)*(ONE-FLES)/FT
!     TOTAL PROBABILITY LEAVING
      PTO=PUP+PDN+PRI+PLE

!     CHECK MASS BUDGET
!     FLUX IN:
!     shift FUP so that use flux entering from below
!      FI(1,2:NZM)=FUP(1,1:NZM-1)-FDN(1,2:NZM)+FLE(2,2:NZM)
!      FI(1,1)=-FDN(1,1)+FLE(2,1)
!      FI(2,2:NZM)=FUP(2,1:NZM-1)-FDN(2,2:NZM)+FLE(3,2:NZM)
!     :           +FRI(1,2:NZM)
!      FI(2,1)=-FDN(2,1)+FLE(3,1)+FRI(1,1)
!      FI(3,2:NZM)=FUP(3,1:NZM-1)-FDN(3,2:NZM)+FRI(2,2:NZM)
!      FI(3,1)=-FDN(3,1)+FRI(2,1)
!     MASS BUDGET
!      FBUD=FT-FI

!      OPEN(35,FILE='flux.out')
!      do k=1,NZM
!      WRITE(35,*)FT(1:3,k),FUP(1:3,k),FDN(1:3,k),FRI(1:3,k),FLE(1:3,k),
!     :  WTT(1:3,k),FI(1:3,k),FBUD(1:3,k),PUP(1:3,k),PDN(1:3,k),
!     :  PRI(1:3,k),PLE(1:3,k)
!      enddo
!      close(35)


      ETIME=0.0
      DTN=DABS(DT)
      IF(BACK)DTN=-DABS(DT)

!     position before conv.: Z1
      IZNDXT=1
      DO K=1,NZM
       IF(ZPROFM(K).LT.Z1)IZNDXT=K+1
      ENDDO

      IF(IZNDXT.EQ.1.AND.Z1.LT.ZPROFT(1))THEN
!     BELOW FIRST LEVEL
!     NOT MOVING, LEAVE
        Z2=Z1
        ICNDX2=ICNDX1
        ZX=IZNDXT
        JET=JET+IDINT(DTN)
        RETURN
      ENDIF

!     Need to allow for change in ICNDX if area
!     changed (previous area cov.: AREAPRU, AREAPRD)
!     Compare previous with current area
!     When area now is smaller, need to detrain
!     Need random variable for this
      AREANOW=AREAA(ICNDX1,IZNDXT)/AREA
      IF(ICNDX1.EQ.1)AREAPR=AREAPRU
      IF(ICNDX1.EQ.3)AREAPR=AREAPRD
      IF(ICNDX1.EQ.2)AREAPR=1.0-AREAPRD-AREAPRU
!      WRITE(*,*)"before",ICNDX1,AREANOW,AREAPR,RSEED,IZNDXT
      IF(AREANOW.GT.0.0001.AND.AREANOW.LT.AREAPR)THEN
        RAND=RAN3(RSEED)
        IF(RAND.GT.AREANOW/AREAPR)THEN
!     detrain from current ICNDX1 with probability
!     p(detr)=(areapr-areanow)/areapr=1-areanow/areapr
!           from up- and downdraft to environment
          IF(ICNDX1.EQ.1.OR.ICNDX1.EQ.3)THEN
            ICNDX1=2
          ELSEIF(ICNDX1.EQ.2)THEN
!         from environment to updraft or downdraft
!         depending on random number and prev. updr. and dndr. coverage
!         transform remaining RAND to random values between 0,1
            RANDN=(RAND-AREANOW/AREAPR)/((AREAPR-AREANOW)/AREAPR)
!         probability to go to updraft: fraction of "new" updraft
!         from total "new" area in up+downdraft
            PDTRUP=(AREAA(1,IZNDXT)/AREA-AREAPRU)/                      &
     &              ((1-AREANOW)-(1-AREAPR))
            IF(AREAA(1,IZNDXT)/AREA.LT.0.0001)PDTRUP=0.0
            IF(RANDN.LT.PDTRUP)ICNDX1=1
            IF(RANDN.GE.PDTRUP)ICNDX1=3
          END IF
        END IF
      END IF
!      WRITE(*,*)"after",ICNDX1,RAND,RANDN,RSEED


!     Advect particle to next grid position (T-grid)
!     need to figure out wether moving up or down

!     if zero flux leaving
!     current level:
!     leave routine, since nothing will happen
      IF(PTO(ICNDX1,IZNDXT).LT.1.1E-14)THEN
      Z2=Z1
      ICNDX2=ICNDX1
      ZX=IZNDXT
      JET=JET+IDINT(DTN)
      RETURN
      ENDIF
!      WRITE(*,*)"2nd: ",z1,IZNDXT,ICNDX1

      IF(Z1.LE.ZPROFT(IZNDXT))THEN
!     BELOW CURRENT LEVEL
        IF(WTT(ICNDX1,IZNDXT-1).GT.1.1E-15)THEN
!     MOVING UP
          ETIME=(ZPROFT(IZNDXT)-Z1)/DABS(WTT(ICNDX1,IZNDXT-1))
          IF(ETIME.GE.DABS(DT)*60.0)THEN
            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*DABS(DT)*60.0
          ENDIF
        ELSEIF(WTT(ICNDX1,IZNDXT-1).LT.-1.1E-15)THEN
!     MOVING DOWN
          IZNDXT=IZNDXT-1
          ETIME=(Z1-ZPROFT(IZNDXT))/DABS(WTT(ICNDX1,IZNDXT))
          IF(ETIME.GE.DABS(DT)*60.0)THEN
            Z2=Z1-DABS(WTT(ICNDX1,IZNDXT))*DABS(DT)*60.0
          ENDIF
        ELSE
!     NOT MOVING, LEAVE
          Z2=Z1
          ICNDX2=ICNDX1
          ZX=IZNDXT
          JET=JET+IDINT(DTN)
          RETURN
        ENDIF
      ELSE
!     ABOVE CURRENT LEVEL
        IF(WTT(ICNDX1,IZNDXT).GT.1.1E-15)THEN
!     MOVING UP
          ETIME=(ZPROFT(IZNDXT+1)-Z1)/DABS(WTT(ICNDX1,IZNDXT))
          IZNDXT=IZNDXT+1
          IF(ETIME.GE.DABS(DT)*60.0)THEN
            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*DABS(DT)*60.0
          ENDIF
        ELSEIF(WTT(ICNDX1,IZNDXT).LT.-1.1E-15)THEN
!     MOVING DOWN
          ETIME=(Z1-ZPROFT(IZNDXT))/DABS(WTT(ICNDX1,IZNDXT))
          IF(ETIME.GE.DABS(DT)*60.0)THEN
            Z2=Z1-DABS(WTT(ICNDX1,IZNDXT))*DABS(DT)*60.0
          ENDIF
        ELSE
!     NOT MOVING, LEAVE
          Z2=Z1
          ICNDX2=ICNDX1
          ZX=IZNDXT
          JET=JET+IDINT(DTN)
          RETURN
        ENDIF
      ENDIF
!      WRITE(*,*)"WTT: ",DABS(WTT(ICNDX1,IZNDXT)),ICNDX1,IZNDXT
!      WRITE(*,*)"P UDRL:",PUP(ICNDX1,IZNDXT),PDN(ICNDX1,IZNDXT),
!     : PRI(ICNDX1,IZNDXT),PLE(ICNDX1,IZNDXT)

      OPEN(34,FILE='mass.out')

!      IF(DABS(WTT(ICNDX1,IZNDXT)).GT.1.1E-15) THEN
        DO WHILE(DABS(DT)*60.0.GE.ETIME)
          ICNDX2=icndx1
          IZNDXT2=IZNDXT
          ETIME2=ETIME
      Z2=ZPROFT(IZNDXT)
!      WRITE(*,*)"WTT2: ",DABS(WTT(ICNDX1,IZNDXT)),ICNDX1,IZNDXT,Z2,ETIME
!     Reassign Cloud index and vert. index
!     get random number
!          RSEED=30985
          RAND=RAN3(RSEED)

!     check if moves to left or right, do immediately, no change in ETIME
          IF(RAND.LT.PLE(ICNDX1,IZNDXT))THEN
            ICNDX1=ICNDX1-1
          ELSEIF(RAND.LT.PLE(ICNDX1,IZNDXT)+PRI(ICNDX1,IZNDXT))THEN
            ICNDX1=ICNDX1+1
!     check if moves up, do it, change ETIME
          ELSEIF(RAND.LT.PLE(ICNDX1,IZNDXT)+PRI(ICNDX1,IZNDXT)          &
     &      +PUP(ICNDX1,IZNDXT)) THEN
              ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/DABS(WTT(ICNDX1,IZNDXT))
              IZNDXT=IZNDXT+1
!       WRITE(*,*)"MOVED UP"
            IF(ETIME.GT.DABS(DT)*60.0)THEN
             Z2=ZPROFT(IZNDXT2)+DABS(WTT(ICNDX1,IZNDXT2))*              &
     &          (DABS(DT)*60.0-ETIME2)
            ENDIF

!     moves down, do it, change ETIME
          ELSE
            IZNDXT=IZNDXT-1
            ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/DABS(WTT(ICNDX1,IZNDXT))
!       WRITE(*,*)"MOVED DN"
            IF(ETIME.GT.DABS(DT)*60.0)THEN
             Z2=ZPROFT(IZNDXT2)-DABS(WTT(ICNDX1,IZNDXT))*               &
     &          (DABS(DT)*60.0-ETIME2)
            ENDIF
          ENDIF

      WRITE(34,*)izndxt,icndx1,ETIME,Z2
!      WRITE(*,*)izndxt,icndx1,ETIME,Z2
!       WRITE(*,*)"RSEED",RSEED,RAND
!      WRITE(*,*)"WTT3: ",DABS(WTT(ICNDX1,IZNDXT)),ICNDX1,IZNDXT,Z2
!     END OF WHILE LOOP
        END DO
!     ASSIGN VALUES AFTER CONVECTION
      ICNDX2=ICNDX1
      ZX=IZNDXT
      JET=JET+IDINT(DTN)

!      ENDIF
      close(34)

      RETURN
      END
