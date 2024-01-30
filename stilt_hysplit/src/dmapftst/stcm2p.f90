      subroutine stcm2p(stcprm, x1, y1, xlat1, xlong1,                  &
     &  x2, y2, xlat2, xlong2)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 x1a,y1a, x2a,y2a, den,dena
!  stcprm->x0 = stcprm->y0 = stcprm->srotate = 0;
      stcprm(11) = 0.
      stcprm(12) = 0.
      stcprm(14) = 0.
!  stcprm->crotate = stcprm->gridszeq = 1.;
      stcprm(13) = 1.
      stcprm(15) = 1.
      call cll2xy(stcprm, xlat1,xlong1, x1a,y1a)
      call cll2xy(stcprm, xlat2,xlong2, x2a,y2a)
      den = dsqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2))
      dena = dsqrt((x1a - x2a)*(x1a - x2a) + (y1a - y2a)*(y1a - y2a))
      stcprm(13) = ((x1a - x2a)*(x1 - x2) + (y1a - y2a)*(y1 - y2) )     &
     &     /den /dena
      stcprm(14) = ((y1a - y2a)*(x1 - x2) - (x1a - x2a)*(y1 - y2) )     &
     &     /den /dena
!  stcprm->gridszeq *= dena /den;
      stcprm(15) = stcprm(15) * dena / den
      call cll2xy(stcprm, xlat1,xlong1, x1a,y1a)
!  stcprm->x0 += x1 - x1a;
!  stcprm->y0 += y1 - y1a;
      stcprm(11) = stcprm(11) + x1 - x1a
      stcprm(12) = stcprm(12) + y1 - y1a
      return
      END
