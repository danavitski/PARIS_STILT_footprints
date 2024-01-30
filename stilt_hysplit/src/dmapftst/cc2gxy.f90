      subroutine cc2gxy(stcprm, x,y, ue,vn, ug,vg)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 norm
      call cpolxy(stcprm, x,y, enx,eny,enz)
      norm = dsqrt(enx*enx + eny*eny)
      if (norm .le. 0.1e-3) then
        call cgrnxy(stcprm, x,y, enx,eny,enz)
        norm = dsqrt(enx*enx + eny*eny)
      endif
      enx = enx/norm
      eny = eny/norm
      ug = enx * vn + eny * ue
      vg = eny * vn - enx * ue
      return
      END

      subroutine cw2gxy(stcprm, x,y, ue,vn, ug,vg)
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
      enx = enx/norm
      eny = eny/norm
      ug = enx * vn + eny * ue
      vg = eny * vn - enx * ue
      return
      END
