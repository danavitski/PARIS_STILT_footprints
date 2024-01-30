      subroutine stvmrc(stcprm, r_lat, r_long)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
!*
!*  Set Map Parameters for a TransVerse MeRCator Projection
!*  Inputs: r_lat,r_longit - latitude and longitude of the reference
!*            point on the meridian tangent to the cylinder,
!*            and 180 degrees from the cut.
!*  Outputs: stcprm - map parameters
!*/
      real*8 refpnt(3),anchrp(3),projpl(3)
      call ll_geo(r_lat, r_long, refpnt)
      TEMP=r_lat - 90.
      call ll_geo(TEMP,r_long,anchrp)
      anorm = x_prod(anchrp, refpnt, projpl)
      call geo_ll(projpl, p_lat, p_long)
      TEMP2=0.
      call mpstrt(stcprm, TEMP2, p_lat,p_long, r_lat, r_long)
      return
      END
