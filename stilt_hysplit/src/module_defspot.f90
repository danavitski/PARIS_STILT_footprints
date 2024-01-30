!======================================================================
!     module_defspot - defines starting location structure
!
!     $Id: module_defspot.f90,v 1.2 2007-05-03 13:09:13 skoerner Exp $
!======================================================================

MODULE module_defspot

   USE module_defsize

!     defines structure of starting location and time
   TYPE RSET

!        starting location            
             REAL*8                  OLAT,OLON
!        starting height (agl)           
             REAL*8                  OLVL
!        emission rate (units/hr)
             REAL*8                  QTRM
!        emission area (m^2)
             REAL*8                  AREA  
!        current x,y,x position
             REAL*8                  XP,YP,ZP
!        previous height for vertical line sources
             REAL*8                  ZV
!        starting time            
             INTEGER                 IBYR,IBMO,IBDA,IBHR 
!        calculation meteo grid number
             INTEGER                 KG

   END TYPE RSET

!======================================================================

   TYPE (RSET)     SPOT(MLOC)

END MODULE module_defspot
