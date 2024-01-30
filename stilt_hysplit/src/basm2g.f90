      subroutine basm2g(stcprm, xin, xout)
!  receives the vector xin, giving the components of a point in the map
!  system; returns xout, the components of the same point in the geo system.
! CHG(03/12/03)
!
! $Id: basm2g.f90,v 1.3 2005-12-14 17:05:57 tnehrkor Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      dimension stcprm(15),xin(3),xout(3)
      ifind(l,k) = -2 + 3*l + k
        do k=1,3
          xout(k) = 0.
        enddo
        do l=1,3
          do k=1,3
            xout(l) = xout(l) + xin(k) * stcprm(ifind(k,l))
          enddo
        enddo
        return
      END
