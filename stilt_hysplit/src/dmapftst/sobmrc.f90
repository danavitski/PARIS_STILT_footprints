      subroutine sobmrc(stcprm, r_lat, r_long, an_lat, an_lng)
!/*
! *  Set Map Parameters for an Oblique MeRCator Projection
! *  Inputs: r_lat,r_long - latitude and longitude of reference
! *            point on the great circle tangent to the cylinder,
! *            and 180 degrees from the cut.
! *       an_lat,an_lng - latitude and longitude of
! *            an anchor point elsewhere on the same great circle.
! *  Outputs: stcprm - map parameters
! */
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 refpnt(3),anchpt(3),prjpol(3)
      real*8 stcprm(15)
      real*8  p_lat,p_long,norm
      call ll_geo(r_lat,r_long,refpnt)
      call ll_geo(an_lat,an_lng,anchpt)
      norm = x_prod(anchpt,refpnt,prjpol)
      if (norm .ne. 0.) then
        call geo_ll(prjpol,p_lat,p_long)
        TEMP=0.
        call mpstrt(stcprm, TEMP, p_lat,p_long, r_lat,r_long)
      else
        call stvmrc(stcprm, r_lat, r_long)
      endif
      return
      END
