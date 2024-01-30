      subroutine cxy2ll(stcprm , x, y, alat, along, proj)
! CHG(03/12/03)
!
! $Id: cxy2ll.f90,v 1.4 2008-03-26 19:16:03 tnehrkor Exp $
!
      use map_utils
      use module_defgrid, only : vmiss, vmissle

      IMPLICIT REAL*8 (A-H,O-Z)
      type(proj_info), intent(in) :: proj
      real*4 :: l_x, l_y, l_alat, l_along

      real*8 stcprm(15)
      real*8 map(3),geog(3)

      if (all(stcprm .le. vmissle)) then
! use WPS mapping routines
         l_x = x
         l_y = y
         call ij_to_latlon(proj, l_x, l_y, l_alat, l_along)
         along = l_along
         alat = l_alat
      else
         call xy_map(stcprm, x,y, map)
         call basm2g(stcprm, map, geog)
         call geo_ll(geog, alat,along)
      endif

      return
      END
