!======================================================================
!     module_defmeto - defines meteorological varaibles that are returned
!     from the advection step at the last particle or puff position.
!----------------------------------------------------------------------
!     $Id: module_defmeto.f90,v 1.3 2009-03-05 08:26:22 gerbig Exp $
!======================================================================

MODULE module_defmeto

  USE module_defsize


!     define structure advection variables for a particle
      TYPE ASET

!       vertical fractional index position
                REAL*8                  ZNDX               
!       height of ground surface
                REAL*8                  ZTER
!       grid spacing (m)
                REAL*8                  GDISX, GDISY
!       trajectory marker variable (press, temp, density)
                REAL*8                  TMRK
!       precipitation rate (meters/min)
                REAL*8                  RAIN
! CHG:(11/20/01) add convective precip rate (m/min)
                REAL*8                  CRAI
! CHG:(11/20/01) add total cloud cover (%)
                REAL*8                  TCLD
! CHG:(12/04/01) add low cloud cover (%)
                REAL*8                  LCLD
! CHG:(11/20/01) add donwn. shortw. radiation (w/m2)
                REAL*8                  DSWF
! CHG:(11/20/01) add sensible heat flux (w/m2)
                REAL*8                  SHTF
! CHG:(11/20/01) add latent heat flux (w/m2)
                REAL*8                  LHTF
! CHG:(22/01/03) add soil moisture content (rel)
                REAL*8                  SOLW
! add surface temperature (T02M)
                REAL*8                  T02M
!       horizontal diffusivity (m2/s)
                REAL*8                  HMIX 
!       aerodynamic roughness length (m)
                REAL*8                  AERO
!       friction velocity (m/s)
                REAL*8                  USTR
!       integrated stability profile
                REAL*8                  PSI
!       pressure profile (mb)
                REAL*8                  PRES(NZM)
!       temperature profile (deg K)
                REAL*8                  TEMP(NZM)
!       relative humidity fraction (0-1)
                REAL*8                  RHFR(NZM)
! JCL:  density profile @ t
                REAL*8                  DENSPREV(NZM)
! CHG:(03/17/2004)      temperature profile @ t
                REAL*8                  TEMPPREV(NZM)
! CHG:(03/17/2004)      rel. hum. profile @ t
                REAL*8                  RHFRPREV(NZM)
! jcl:  this is density profile @ t+dt
!       density profile (kg/m3)
                REAL*8                  DENS(NZM)
! JCL:(6/28/02) local air density (kg/m3)
                REAL*8                  DENSLOCAL
!       vertical mixing profile (m2/s)
                REAL*8                  VMIX(NZM)
! JCL:  Lagrangian timescale (s) @ t+dt
                REAL*8                  TL(NZM)
! JCL:  std dev of vertical velocity (m/s) @ t+dt
                REAL*8                  SIGW(NZM)
! JCL:  Lagrangian timescale (s) @ t
                REAL*8                  TLPREV(NZM) 
! JCL:  std dev of vertical velocity (m/s) @ t
                REAL*8                  SIGWPREV(NZM)
! JCL:  mixed layer heights (m) @ t
                REAL*8                  ZMLPREV
! JCL:  mixed layer heights (m) @ t+dt
                REAL*8                  ZMLNEXT
! CHG:(12/04/01)  lim. of convection heights (m) @ t
                REAL*8                  ZLOCPREV
! CHG:(12/04/01)  lim. of convection heights (m) @ t+dt
                REAL*8                  ZLOCNEXT
! JCL:(5/9/01) horizontal velocity (grid/min) @ t
                REAL*8                  UUPREV(NZM)
! JCL:(5/9/01) horizontal velocity (grid/min) @ t+dt
                REAL*8                  UUNEXT(NZM)
! JCL:(5/9/01) horizontal velocity (grid/min) @ t
                REAL*8                  VVPREV(NZM)
! JCL:(5/9/01) horizontal velocity (grid/min) @ t+dt
                REAL*8                  VVNEXT(NZM)
! JCL:(4/3/02) mass violation (grid/min) @ t
                REAL*8                DMASSPREV(NZM)
! JCL:(4/3/02) mass violation (grid/min) @ t+dt
                REAL*8                DMASSNEXT(NZM)
! CHG:(09/23/03) convective fluxes (RAMS, ECMWF, ...)
                REAL*8                CFXUP1(NZM)
                REAL*8                CFXUP2(NZM)
                REAL*8                CFXDN1(NZM)
                REAL*8                DFXUP1(NZM)
                REAL*8                DFXUP2(NZM)
                REAL*8                EFXUP1(NZM)
                REAL*8                EFXUP2(NZM)
                REAL*8                DFXDN1(NZM)
                REAL*8                EFXDN1(NZM)
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
                REAL*8                TKEN(NZM)
! CHG:(09/23/03) relative area coverages from RAMS
                REAL*8                RAUP1(NZM)
                REAL*8                RAUP2(NZM)
                REAL*8                RADN1(NZM)
!       true position of endpoint
                REAL*8                  PLAT,PLON

!       current sub-grid corner position 
                INTEGER                 LX1,LY1
!       potential sub-grid corner position 
                INTEGER                 LXC,LYC
!       potential sub-grid size or range  
                INTEGER                 LXR,LYR
!       land-use category
                INTEGER                 LAND

      END TYPE ASET

!======================================================================

      TYPE (ASET) METO

END MODULE module_defmeto
