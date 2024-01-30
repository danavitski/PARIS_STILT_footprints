      subroutine ccrvxy(stcprm, x, y, gx, gy)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (REARTH=6367.47)
      real*8 stcprm(15)
      real*8 map(3)
      call xy_map(stcprm, x,y, map)
      if (dabs(map(3)) .ge. 1.) then
        if (stcprm(1) .eq. map(3)) then
          gx = 0.
          gy = 0.
        else
          xi = 0.
          TEMP=.5e21
          eta = dsign(TEMP,map(3)) / REARTH
        endif
      else
        xi0 = (x - stcprm(11)) * stcprm(15) / REARTH
        eta0 = (y - stcprm(12)) * stcprm(15) / REARTH
        xi = xi0 * stcprm(13) - eta0 * stcprm(14)
        eta = xi0 * stcprm(14) + eta0 * stcprm(13)
        xi = - stcprm(1) * xi
        eta = 1. - stcprm(1) * eta
        fact = dsqrt((xi*xi + eta*eta) / (1. - map(3)) / (1. + map(3)))
        xi = xi/fact
        eta = eta/fact
        fact = (stcprm(1) - map(3)) / (1. - map(3)) / (1. + map(3))     &
     &     / REARTH
        xi = xi * fact
        eta = eta * fact
      endif
      gx = xi * stcprm(13) + eta * stcprm(14)
      gy = eta * stcprm(13) - xi * stcprm(14)
      return
      END
