      subroutine cgrnxy(stcprm, x,y, enx,eny,enz)
! returns a vector aligned with the direction toward the North Pole.  I.e.
! parallel to the vector from the Earth's center to the North Pole.
! CHG(03/12/03)
!
! $Id: cgrnxy.f90,v 1.2 2005-02-02 19:39:49 skoerner Exp $
!
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      real*8 map(3),pole(3)
      call xy_map(stcprm, x,y, map)
      do k=1,3
        pole(k) = stcprm(3*k - 1)
      enddo
      call proj_3d(stcprm, map, pole, enx,eny,enz)
      return
      END
