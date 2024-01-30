      subroutine cg2cxy(stcprm, x,y, ug,vg, ue,vn)
! CHG(03/12/03)
!
! $Id: cg2cxy.f90,v 1.3 2005-04-13 13:48:35 skoerner Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (sin1dg=.017452406437284d0)
      real*8 stcprm(15)
      real*8 norm
      call cpolxy(stcprm, x,y, enx,eny,enz)
      norm = dsqrt(enx*enx + eny*eny)
      if (norm .le. 0.1e-3) then
        call cgrnxy(stcprm, x,y, enx,eny,enz)
        norm = dsqrt(enx*enx + eny*eny)
      endif
      enx = enx / norm
      eny = eny / norm
      ue = - enx * vg + eny * ug
      vn = eny * vg + enx * ug
      return
      END

      subroutine cg2wxy(stcprm, x,y, ug,vg, ue,vn)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (sin1dg=.017452406437284)
      real*8 stcprm(15)
      real*8 norm
      call cpolxy(stcprm, x,y, enx,eny,enz)
      norm = dsqrt(enx*enx + eny*eny)
      if (norm .le. sin1dg) then
        call cgrnxy(stcprm, x,y, enx,eny,enz)
        norm = dsqrt(enx*enx + eny*eny)
      endif
      enx = enx / norm
      eny = eny / norm
      ue = - enx * vg + eny * ug
      vn = eny * vg + enx * ug
      return
      END
