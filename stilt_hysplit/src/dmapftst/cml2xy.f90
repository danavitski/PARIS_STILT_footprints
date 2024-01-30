      integer function cml2xy(stcprm, alat, along, lsiz, x, y )
      real*8 stcprm(15)
      real*8 alat(*),along(*),x(*),y(*)
      real*8 geog(3), mapold(3), mapnew(3), mapmid(3), fact, sizmid
      real*8 xi, eta, a, b
      logical cross
      integer kount,npoints,k
! CHG(03/12/03)
!      IMPLICIT REAL*8 (A-H,O-Z)

      cross(a,b) = ((a .le. 0.) .and. (b .gt. 0.))                      &
     & .or. ((b.le.0.) .and. (a .gt. 0.))
!* Given a string of latitudes and longitudes, returns a string
!* of x-y values.  In the event that the *first* segment crosses the
!* 'cut' of the map, from side a to side b, the first x-y returned
!* will be the side b version of the cut-crossing point.  (If size_l<0,
!* then only the first part of the segment, on side a, will be returned.)
!*  If *any other* segment crosses the cut, this ends the transfer, the
!* last x-y value will be the near side of cut, and the returned value
!* will be negative.  Returns an int, the absolute value of which is
!* the number of points transformed, and the sign of which indicates
!* (if negative) that more points remain to be transformed.*/
      if (lsiz .eq. 0) then
!  No points provided; none returned.
        cml2xy = 0
        return
      endif
      kount=1
      npoints = iabs(lsiz)
      call ll_geo(alat(1),along(1), geog)
      call basg2m(stcprm, geog, mapold)
      call map_xy(stcprm, mapold, x(1),y(1))
      if (npoints .eq. 1) then
!* Only one point provided; one point returned.
        cml2xy = 1
        return
      endif
      if (lsiz .gt. 0) then
!* If first segment crossed, truncate first section.
        call ll_geo(alat(2),along(2), geog)
        call basg2m(stcprm, geog, mapnew)
        if (cross(mapold(2),mapnew(2))) then
          fact = mapnew(2) / (mapnew(2) - mapold(2))
          mapmid(2)=0.
          mapmid(1) = mapnew(1) + fact * (mapold(1) - mapnew(1))
          mapmid(3) = mapnew(3) + fact * (mapold(3) - mapnew(3))
          if (mapmid(1) .le. 0. ) then
!* only in this case is there a crossing of the cut.
            sizmid = dsqrt(mapmid(1)**2 + mapmid(3)**2)
            if (sizmid .gt. 0.) then
              mapmid(1) = mapmid(1) / sizmid
              mapmid(3) = mapmid(3) / sizmid
            else
              mapmid(1) = -1.
            endif
            if (mapnew(2) .gt. 0.) then
              call map_xe(stcprm, mapmid, xi,eta, 'e')
            else
              call map_xe(stcprm, mapmid, xi,eta, 'w')
            endif
            call xe_xy(stcprm, xi,eta, x(1),y(1))
          endif
        endif
        call map_xy(stcprm, mapnew, x(2), y(2))
        do k=1,3
          mapold(k) = mapnew(k)
        enddo
        kount=kount+1
      endif
!* geog values indicate v[0] "Toward lat0,lon0", v[1] East @0,0,
!* v[2] North @0,0.  Thus cut half-plane is v[1]=0, v[0]<=0.
!* Interpolate.*/
      do kount=kount+1,npoints
        call ll_geo(alat(kount),along(kount),geog)
        call basg2m(stcprm, geog, mapnew)
        if ( cross(mapold(2),mapnew(2)) ) then
!* Place interpolated value in location count and exit*/
!*If cut crossed, replace xy-point in location 0 with interpolated value*/
          fact = mapnew(2) / (mapnew(2) - mapold(2))
          mapmid(2)=0.
          mapmid(1) = mapnew(1) + fact * (mapold(1) - mapnew(1))
          mapmid(3) = mapnew(3) + fact * (mapold(3) - mapnew(3))
          if (mapmid(1) .le. 0. ) then
!* only in this case is there a crossing of the cut.
            sizmid = dsqrt(mapmid(1)**2 + mapmid(3)**2)
            if (sizmid .gt. 0.) then
              mapmid(1) = mapmid(1) / sizmid
              mapmid(3) = mapmid(3) / sizmid
            else
              mapmid(1) = -1.
            endif
          endif
          if (mapold(2) .gt. 0.) then
            call map_xe(stcprm, mapmid, xi,eta, 'e')
          else
            call map_xe(stcprm, mapmid, xi,eta, 'w')
          endif
          call xe_xy(stcprm, xi,eta, x(kount),y(kount))
! CHG (03/12/03) typo?
!          ml2xy = -kount
          cml2xy = -kount
          return
        endif
        call map_xy(stcprm, mapnew, x(kount),y(kount))
        do k=1,3
          mapold(k) = mapnew(k)
        enddo
      enddo
! CHG (03/12/03) typo?
!      ml2xy = npoints
      cml2xy = npoints
      return
      END
