      real*8 function cgszll(stcprm, alat, along, grid)
!
! $Id: cgszll.f90,v 1.4 2008-03-26 19:16:02 tnehrkor Exp $
!
      use map_utils
      use module_defgrid, only : gset, vmiss, vmissle

      IMPLICIT NONE

      type(gset) ,intent(in) :: grid

      real*8 stcprm(15), alat, along

      real*8 map(3)
      real*8 geog(3)
      real*8 ymerc

      integer :: i_mapfactor, j_mapfactor
      real :: x_alat, x_along, xi, xj

      if (all(stcprm .le. vmissle)) then
         ! WPS mapping routines use precomputed mapfactor array
         x_alat = alat
         x_along = along
         call latlon_to_ij(grid%proj,x_alat,x_along,xi,xj)
         i_mapfactor = nint(xi)
         i_mapfactor = min(grid%nx,max(1,i_mapfactor))
         j_mapfactor = nint(xj)
         j_mapfactor = min(grid%ny,max(1,j_mapfactor))

         cgszll = grid%proj%dx * grid%mapfactor(i_mapfactor,j_mapfactor) / 1000.

      else
         call ll_geo(alat, along, geog)
         call basg2m(stcprm, geog, map)
         if (map(3) .ge. 1.) then
            if (stcprm(1) .ge. 1.) then
               cgszll = 2. * stcprm(15)
            else
               cgszll = 0.
            endif
         else if (map(3) .le. -1.) then
            if (stcprm(1) .le. -1.) then
               cgszll = 2. * stcprm(15)
            else
               cgszll = 0.
            endif
         else if (dabs(stcprm(1)) .ge. 1.) then
            cgszll = stcprm(15) * (1. + dsign(map(3),stcprm(1)))
         else
            ymerc = -.5 * dlog( (1. - map(3))/(1. + map(3)) )
            cgszll = stcprm(15) * dexp( - (1.-stcprm(1)) * ymerc) *         &
                 &   (1. + map(3))
         endif
      end if

      return
      END
