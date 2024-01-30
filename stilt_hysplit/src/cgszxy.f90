      real*8 function cgszxy(stcprm, x, y, grid)

      use map_utils
      use module_defgrid, only : gset, vmiss, vmissle

      type(gset) ,intent(in) :: grid

      real*8 stcprm(15)
      real*8 map(3)
      real*8 x,y,ymerc

      integer :: i_mapfactor, j_mapfactor

! CHG(03/12/03)
!
! $Id: cgszxy.f90,v 1.4 2008-03-26 19:16:02 tnehrkor Exp $
!
!      IMPLICIT REAL*8 (A-H,O-Z)

      if (all(stcprm .le. vmissle)) then
         ! WPS mapping routines use precomputed mapfactor array
         i_mapfactor = nint(x)
         i_mapfactor = min(grid%nx,max(1,i_mapfactor))
         j_mapfactor = nint(y)
         j_mapfactor = min(grid%ny,max(1,j_mapfactor))

         cgszxy = grid%proj%dx * grid%mapfactor(i_mapfactor,j_mapfactor) / 1000.

      else
         call xy_map(stcprm, x,y,map)
         if (map(3) .ge. 1.)  then
            if (stcprm(1) .ge. 1.) then
               cgszxy = 2.*stcprm(15)
            else
               cgszxy = 0.
            endif
         else if (map(3) .le. -1.) then
            if (stcprm(1) .le. -1.) then
               cgszxy = 2.*stcprm(15)
            else
               cgszxy = 0.
            endif
         else if (dabs(stcprm(1)) .ge. 1.) then
            cgszxy = stcprm(15) * (1. + dsign(map(3),stcprm(1)))
         else
            ymerc = -.5 * dlog( (1. - map(3))/(1. + map(3)) )
            cgszxy = stcprm(15) * dexp( - (1. - stcprm(1)) * ymerc) *       &
                 &      (1. + map(3))
         endif
      end if
      
      return
      END
