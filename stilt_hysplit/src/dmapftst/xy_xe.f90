      subroutine xy_xe(stcprm,  x,y, xi,eta)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      PARAMETER (REARTH=6367.47)
        xtmp = x - stcprm(11)
        ytmp = y - stcprm(12)
        xi = stcprm(15) / REARTH *                                      &
     &    (stcprm(13) * xtmp - stcprm(14) * ytmp)
        eta = stcprm(15) / REARTH *                                     &
     &    (stcprm(13) * ytmp + stcprm(14) * xtmp)
        return
      END
