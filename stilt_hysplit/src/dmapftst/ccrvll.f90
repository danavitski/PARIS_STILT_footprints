      subroutine ccrvll(stcprm, alat, along, gx, gy)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (REARTH=6367.47)
      real*8 stcprm(15)
      real*8 map(3),geog(3)
      call ll_geo(alat,along, geog)
      call basg2m(stcprm, geog, map)
      if ( dabs(map(3)) .ge. 1.) then
        if (stcprm(1) .eq. map(3)) then
          gx = 0.
          gy = 0.
          return
        else
          xi = 0.
          TEMP=.5E21
          eta = dsign(TEMP,map(3)) / REARTH
        endif
      else
        fact = (stcprm(1) - map(3)) /                                   &
     &    (1. - map(3)) / (1. + map(3)) / REARTH
        xi = - map(2) * fact
        eta = map(1) * fact
        if (dabs(stcprm(1)) .lt. 1.) then
          glambda = (stcprm(1) - 1.) * datan2(map(2), map(1))
          clambda = dcos(glambda)
          slambda = dsin(glambda)
          fact = xi * clambda - eta * slambda
          eta = xi * slambda + eta * clambda
          xi = fact
        endif
      endif
      gx = xi * stcprm(13) + eta * stcprm(14)
      gy = eta * stcprm(13) - xi * stcprm(14)
      return
      END
