!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVRNG           ADVection RaNGe to determine subgrid
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:99-03-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION RANGE COMPUTES THE RANGE OF PARTICLE POSITIONS ON THE
!   COMPUTATIONAL GRID (METEO GRID) WHICH IS USED TO DETERMINE THE
!   OPTIMUM LOCATION AND SIZE OF THE METEOROLOGICAL SUBGRID. THE
!   SUBGRID SIZE IS ALSO DETERMINED BY THE FREQUENCY OF THE METEO DATA
!   IN THAT IT IS NECESSARY TO CONSIDER HOW LONG A PARTICLE REMAINS ON
!   THE SUBGRID BEFORE MORE DATA IS REQUIRED TO BE LOADED.  IT IS BETTER
!   TO AVOID LOADING DATA MORE THAN ONCE FOR ANY TIME PERIOD.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Mar 1999 (RRD) - Initial version of the code
!
! USAGE:  CALL ADVRNG(KG,UMAX,KPM,XPOS,YPOS)
!   INPUT ARGUMENT LIST:
!     KG    - int current meteorological grid
!     UMAX  - real      maximum wind speed (grid pts / min)
!     KPM   - int number of points
!     XPOS,YPOS   - real      puff center positions (grid units)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: advrng.f90,v 1.6 2007-02-16 17:53:46 tnehrkor Exp $
!
!$$$

      SUBROUTINE ADVRNG(KG,UMAX,KPM,XPOS,YPOS)

      use module_defmeto
      use module_defgrid

      IMPLICIT REAL*8 (A-H,O-Z)

!     array sizes
!     meteorology grid and file
!      INCLUDE 'DEFGRID.INC'
!     meteo variables returned after advection
!      INCLUDE 'DEFMETO.INC'

!     particle positions
      REAL*8 XPOS(*),YPOS(*)

!     associated with previous include files
!      COMMON /GBLMET/ METO
!      COMMON /GBLGRD/ GRID, DREC, FILE

!     check for sufficient number
      IF(KPM.LE.0)RETURN

!     zero out summation
      XSUM=0.0
      YSUM=0.0
      CNTR=0.0
      LXMAX=1
      LXMIN=NXM
      LYMAX=1
      LYMIN=NYM

!     sum positions and deterime max/min
      DO KP=1,KPM
         XSUM=XSUM+XPOS(KP)
         YSUM=YSUM+YPOS(KP)
         LXMAX=MAX0(LXMAX,IDNINT(XPOS(KP)))
         LXMIN=MIN0(LXMIN,IDNINT(XPOS(KP)))
         LYMAX=MAX0(LYMAX,IDNINT(YPOS(KP)))
         LYMIN=MIN0(LYMIN,IDNINT(YPOS(KP)))
         CNTR=CNTR+1.0
      END DO

!     compute mean particle position
! CHG(10/03/03) change center position for RAMS to account for stagger
      IF(GRID(KG)%MODEL_ID.EQ.'RAMS')THEN
        METO%LXC=IDNINT(XSUM/CNTR)
        METO%LYC=IDNINT(YSUM/CNTR)
      ELSE
        METO%LXC=XSUM/CNTR
        METO%LYC=YSUM/CNTR
      ENDIF

!     transport distance = (grid units / min ) (meteo time interval)
! CHG(09/26/03) some metfiles have intervals shorter than 60 minutes (e.g. RAMS)
!                60 min is chosen because hymodelc.f has 60-min timestep loop
!     DIST=UMAX*DREC(KG)%DELTA
      DIST=UMAX*MAX0(DREC(KG)%DELTA,60)
!     minimum subgrid size should be twice the distance since the subgrid
!     will be centered over the mean particle position
      NGRD=NINT(2.0*DIST)

!     the subgrid should be at least 10 grid units or at least the size
!     determined from the advection distance plus particle distribution
! CHG&JCL (03/10/2004) change minimum subgrid dimension from 10
!     a value depending on grid size
!     at 40 km and larger still keep 10, at smaller grid sizes keep
!     physical dimension constant
      MING=MAX(10,IDNINT(10.0*40.0/GRID(KG)%SIZE))
      METO%LXR=MIN0(MAX0(MING,NGRD+3*(LXMAX-LXMIN)),NXM)
      METO%LYR=MIN0(MAX0(MING,NGRD+3*(LYMAX-LYMIN)),NYM)

! JCL:(07/12/2004) add lines from HYSPLIT Vers 45 to deal with global grids
!     when subgrid reaches 75% then set limits to maximum
      IF(METO%LXR.GT.NINT(0.75*GRID(KG)%NX).OR.GRID(KG)%GLOBAL)THEN
        METO%LXR=GRID(KG)%NX
        METO%LXC=NINT((GRID(KG)%NX+1)/2.0)
      ELSE
!       compute grid center position
        METO%LXC=NINT((LXMAX+LXMIN)/2.0)
      END IF

!     when subgrid reaches 75% then set limits to maximum
      IF(METO%LYR.GT.NINT(0.75*GRID(KG)%NY).OR.GRID(KG)%GLOBAL)THEN
        METO%LYR=GRID(KG)%NY
        METO%LYC=NINT((GRID(KG)%NY+1)/2.0)
      ELSE
!       compute grid center position
        METO%LYC=NINT((LYMAX+LYMIN)/2.0)
      END IF


!     optional diagnostic
!     WRITE(30,*)'Subgrid optimization: ',NGRD,METO.LXR,METO.LYR
!      WRITE(*,*)'Subgrid optimization: ',NGRD,METO%LYR,
!     &   METO%LYC,DIST,UMAX,YSUM/CNTR

      RETURN
      END
