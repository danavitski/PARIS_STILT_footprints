      subroutine map_xy(stcprm, x_map, x, y)
! CHG(03/12/03)
!
! $Id: map_xy.f90,v 1.3 2005-12-14 17:05:58 tnehrkor Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension stcprm(15),x_map(3)
        call map_xe(stcprm, x_map, xi, eta, 'c')
        call xe_xy(stcprm, xi, eta, x, y)
        return
      END
