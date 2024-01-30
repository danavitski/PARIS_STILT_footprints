      subroutine soblmb(stcprm, r_lat, r_lng, a_lat1, a_lng1,           &
     &    a_lat2, a_lng2)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      parameter (pi=3.14159265358979323846,radpdg=pi/180.)
      parameter (dgprad=180./pi)
      real*8 a_lat1,a_lng1,r_lat,r_lng,a_lat2,a_lng2
      real*8 stcprm(15)
!/*
! *  Set Map Parameters for an Oblique LaMBeRt Conic Conformal
! *          Projection.
! *  Inputs: r_lat,r_lng,
! *          a_lat1,a_lng1,
! *          a_lat2,a_lng2  - latitudes and longitudes of three
! *              points on the circle (great or small) which is
! *              tangent to the projection cone.  Point r_lat,r_lng is
! *              the reference point; 180degrees away from the cut.
! *  Outputs: stcprm - map parameters
! */
      real*8 pnt_a1(3),pnt_rf(3),pnt_a2(3),pol_pt(3)
      real*8 b_a(3),c_b(3),norm,sin_cn,cos_cn,temp,lat_pl,lng_pl
      call ll_geo(a_lat1,a_lng1,pnt_a1)
      call ll_geo(r_lat,r_lng,pnt_rf)
      call ll_geo(a_lat2,a_lng2,pnt_a2)
      do k=1,3
        b_a(k) = pnt_rf(k) - pnt_a1(k)
        c_b(k) = pnt_a2(k) - pnt_rf(k)
      enddo
      norm = x_prod(b_a,c_b,pol_pt)
      if (norm .ne. 0.) then
        call geo_ll(pol_pt,lat_pl,lng_pl)
        sin_cn = 0.
        do k=1,3
          pol_pt(k) = pol_pt(k)/norm
          sin_cn = sin_cn + pnt_rf(k) * pol_pt(k)
        enddo
        cos_cn = 0.
        do k=1,3
          temp = sin_cn *pol_pt(k) - pnt_rf(k)
          cos_cn = cos_cn + temp*temp
        enddo
        call mpstrt(stcprm,dgprad*datan2(sin_cn,cos_cn),lat_pl,lng_pl,  &
     &        r_lat,r_lng)
      else
        temp = 0.
        do k=1,3
          temp = temp + b_a(k) * b_a(k)
        enddo
        if (temp .ne. 0.) then
          call sobmrc(stcprm,r_lat,r_lng,a_lat1,a_lng1)
        else
          do k=1,3
            temp = temp + c_b(k)*c_b(k)
          enddo
          if (temp .ne. 0.) then
            call sobmrc(stcprm,r_lat,r_lng, a_lat1,a_lng1)
          else
            call sobstr(stcprm, r_lat,r_lng)
          endif
        endif
      endif
      return
      END
