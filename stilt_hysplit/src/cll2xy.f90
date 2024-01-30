      subroutine cll2xy(stcprm, alat, along, x, y, proj )
! CHG(03/12/03)
!
! $Id: cll2xy.f90,v 1.4 2008-03-26 19:16:02 tnehrkor Exp $
!
      use map_utils
      use module_defgrid, only : vmiss, vmissle

      IMPLICIT REAL*8 (A-H,O-Z)
      type(proj_info), intent(in) :: proj
      real*4 :: l_x, l_y, l_alat, l_along

      real*8 stcprm(15)
      real*8 geog(3),temp(3)


      if (all(stcprm .le. vmissle)) then
! use WPS mapping routines
         l_along = along
         l_alat = alat
         call latlon_to_ij(proj, l_alat, l_along, l_x, l_y)
         x = l_x
         y = l_y
      else
         call ll_geo(alat,along, geog)
         call basg2m(stcprm, geog, temp)
         call map_xy(stcprm, temp, x,y)
      endif

      return
      END
