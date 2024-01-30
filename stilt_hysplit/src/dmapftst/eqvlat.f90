      real*8 function eqvlat(lat1,lat2)
!*  Written 12/21/94 by Dr. Albion Taylor
!*
!*    This function is provided to assist in finding the tangent latitude
!*    equivalent to the 2-reference latitude specification in the legend
!*    of most lambert conformal maps.  If the map specifies "scale
!*    1:xxxxx true at 40N and 60N", then eqvlat(40.,60.) will return the
!*    equivalent tangent latitude.
!*  INPUTS:
!*    lat1, lat2:  The two latitudes specified in the map legend
!*  RETURNS:
!*    the equivalent tangent latitude
!*  EXAMPLE:  stcmap(& strcmp, eqvlat(40.,60.), 90.)
!*/
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 lat1,lat2
      parameter (pi=3.14159265358979323846,radpdg=pi/180.)
      parameter (dgprad=180./pi)
      parameter (fsm = 1.e-3)
        slat1 = dsin(radpdg * lat1)
        slat2 = dsin(radpdg * lat2)
        if (slat1 + slat2 .eq. 0.) then
          eqvlat = 0.
          return
        endif
        if (dabs(slat1) .ge. 1.) then
          TEMP=90.
          eqvlat = dsign(TEMP,slat1)
          return
        endif
        if (dabs(slat2) .ge. 1.) then
          TEMP=90.
          eqvlat = dsign(TEMP,slat2)
        endif
        if ( dabs(slat1 - slat2) .gt. fsm) then
          al1 = dlog( (1. - slat1) / (1. - slat2) )
          al2 = dlog( (1. + slat1) / (1. + slat2) )
        else
          tau = - (slat1 - slat2)/(2. - slat1 - slat2)
          tau = tau * tau
          al1 = 2./(2. - slat1 - slat2) * (1.    + tau *                &
     &                                    (1./3. + tau *                &
     &                                    (1./5. + tau *                &
     &                                    (1./7. + tau))))
          tau = (slat1 - slat2)/(2. + slat1 + slat2)
          tau = tau * tau
          al2 = -2./(2. + slat1 + slat2) * (1.    + tau *               &
     &                                     (1./3. + tau *               &
     &                                     (1./5. + tau *               &
     &                                     (1./7. + tau))))
          endif
          eqvlat = dasin ((al1 + al2) / (al1 - al2)) * DGPRAD
          return
      END
