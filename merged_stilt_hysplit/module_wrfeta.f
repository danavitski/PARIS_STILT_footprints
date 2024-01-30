module module_wrfeta

  ! This file contains subroutines for dealing with WRF sigma-p (eta) vertical coordinates
  !   PRGMMR:    Thomas Nehrkorn  AER, Inc.  August 2012
  ! $Id: module_wrfeta.f,v 1.2 2012/09/05 16:38:27 trn Exp $
  ! PROGRAM HISTORY LOG:
  !   $Log: module_wrfeta.f,v $
  !   Revision 1.2  2012/09/05 16:38:27  trn
  !   Merge changes from stilt_sigmap branch, at tag SLC_hum_tested
  !
  !   Revision 1.1.2.5  2012/09/05 16:23:36  trn
  !   Moved compute_zsgfull calls outside module_wrfeta for speedup
  !
  !   Revision 1.1.2.4  2012/08/31 14:51:28  trn
  !   Removed debugging code
  !
  !   Revision 1.1.2.3  2012/08/31 14:43:31  trn
  !   Bug fixes, replace allocatable by automatic arrays
  !
  !   Revision 1.1.2.2  2012/08/09 20:40:39  trn
  !   Intermediate checkin: successfull compilation and linking
  !
  !   Revision 1.1.2.1  2012/08/09 17:06:59  trn
  !   Initial set of mods for wrf-sigmap support - before compilation and testing
  !

  private

  public :: eta2agl, agl2eta, compute_zsgfull

  real :: GRAV=9.80616

contains


  subroutine eta2agl(etain,zaglout,mu,alt0,zsg,nlvl,zsg_full,dzsg_full,zlvls,zfull,compute_zlvls)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    ! SUBPROGRAM:  eta2agl           compute zagl from eta
    !
    ! ABSTRACT:  
    !   Computes internal model level heights (as m AGL) from the
    !   model level WRF eta (sigma-p) coordinates

    implicit none

    ! Interface:
    real         , INTENT(IN)    :: etain             !input eta coordinate value (dimensionless)
    real         , INTENT(out)    :: zaglout             !output AGL height (m)
    real         , INTENT(IN)    :: mu             !WRF mu=pdh_sfc-pdh_top
    real         , INTENT(IN)    :: alt0(:)             !WRF inverse density alpha
    real         , INTENT(IN)    :: zsg(:)             !model layer eta values (bottom-up)
    real         , INTENT(IN)    :: zsg_full(:)
    real         , INTENT(IN)    :: dzsg_full(:)
    INTEGER      , INTENT(IN) :: nlvl          !number of model layers
    REAL         , INTENT(INOUT) :: zlvls(:)  !model layer AGL values (m)
    REAL         , INTENT(INOUT) :: zfull(:)  !model level interface AGL values (m) (zfull(0)=0, zfull(nlvl)=top)
    LOGICAL      , INTENT(IN) :: compute_zlvls  !flag (T/F) for computing (T) or using (F) the zlvls and zfull

    ! Other local variables
    INTEGER            :: k,l

    if (compute_zlvls) &
         call compute_allagl(mu,alt0,zsg,nlvl,zlvls,zfull,zsg_full,dzsg_full)
    ! find vertical index of layer:
    l=1
    do k=nlvl-1,1,-1
       if (etain .lt. zsg_full(k)) then
          ! point is in layer l=k+1, between zsg_full(k) and zsg_full(k+1)
          l=k+1
          exit
       end if
    end do
    ! compute zagl using WRF hydrostatic eqn: (same as in prfwrf)
    if (l .lt. 2) then
       ! layer 1: between AGL=0,zsg=1 and AGL=zfull(1),zsg=zsg_full(1)
       zaglout=-mu*(etain-1.)*alt0(1)/grav
    else
       ! layer 2-nlvl: between AGL=zfull(l-1),zsg=zsg_full(l-1) and AGL=zfull(l),zsg=zsg_full(l)
       zaglout=zfull(l-1)-mu*(etain-zsg_full(l-1))*alt0(l)/grav
    end if
    return
    !---------------------------------------------------------------------------------------------------
  end subroutine eta2agl

  subroutine agl2eta(zaglin,etaout,mu,alt0,zsg,nlvl,zsg_full,dzsg_full,zlvls,zfull,compute_zlvls)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    ! SUBPROGRAM:  eta2agl           compute zagl from eta
    !
    ! ABSTRACT:  
    !   Computes internal model level heights (as m AGL) from the
    !   model level WRF eta (sigma-p) coordinates

    implicit none

    ! Interface:
    real         , INTENT(IN)    :: zaglin             !input AGL height (m)
    real         , INTENT(out)    :: etaout             !output eta coordinate value (dimensionless)
    real         , INTENT(IN)    :: mu             !WRF mu=pdh_sfc-pdh_top
    real         , INTENT(IN)    :: alt0(:)             !WRF inverse density alpha
    real         , INTENT(IN)    :: zsg(:)             !model layer eta values (bottom-up)
    real         , INTENT(IN)    :: zsg_full(:)
    real         , INTENT(IN)    :: dzsg_full(:)
    INTEGER      , INTENT(IN) :: nlvl          !number of model layers
    REAL         , INTENT(INOUT) :: zlvls(:)  !model layer AGL values (m)
    REAL         , INTENT(INOUT) :: zfull(:)  !model level interface AGL values (m) (zfull(1)=top of first layer, zfull(nlvl)=top)
    LOGICAL      , INTENT(IN) :: compute_zlvls  !flag (T/F) for computing (T) or using (F) the zlvls and zfull

    ! Other local variables
    INTEGER            :: k,l
    real :: zlower, zsglower, tmp_alt0

    if (compute_zlvls) &
         call compute_allagl(mu,alt0,zsg,nlvl,zlvls,zfull,zsg_full,dzsg_full)

    ! find vertical index of layer:
    l=nlvl+1
    do k=1,nlvl
       if (zaglin .lt. zfull(k)) then
          ! point is in layer l=k, between zfull(k-1) and zfull(k)
          l=k
          exit
       end if
    end do
    ! compute eta by applying WRF hydrostatic eqn between layer lower boundary and zaglin:
    ! for lowest layer, lower boundary is at zagl=0:
    zlower=0.
    zsglower=1.
    if (l .gt. 1) then
       zlower=zfull(l-1)
       zsglower=zsg_full(l-1)
    end if
    ! Use top-layer alt0 for zagl > model top:
    tmp_alt0=alt0(min(nlvl,l))
    etaout=zsglower - grav*(zaglin-zlower)/(mu*tmp_alt0)
    return
    !---------------------------------------------------------------------------------------------------
  end subroutine agl2eta
  
  subroutine compute_zsgfull(zsg,nlvl,zsg_full,dzsg_full)

    implicit none

    real         , INTENT(IN)    :: zsg(:)             !model layer eta values (bottom-up)
    INTEGER      , INTENT(IN) :: nlvl          !number of model layers
    REAL         , INTENT(OUT) :: dzsg_full(:)  !model layer thicknesses (eta)
    REAL         , INTENT(OUT) :: zsg_full(:)  !model interface level eta values (starting at top of bottom layer)

    integer :: k

    ! compute layer thicknesses and interface level etas (same as in prfcom)
    dzsg_full(1)=2*(zsg(1)-1.)
    zsg_full(1)=1.+dzsg_full(1)
    do k=2,nlvl
       dzsg_full(k)=2*(zsg(k)-zsg_full(k-1))
       zsg_full(k)=zsg_full(k-1)+dzsg_full(k)
    end do
    return
  end subroutine compute_zsgfull
  
  subroutine compute_allagl(mu,alt0,zsg,nlvl,zlvls,zfull,zsg_full,dzsg_full)

    implicit none

    real         , INTENT(IN)    :: mu             !WRF mu=pdh_sfc-pdh_top
    real         , INTENT(IN)    :: alt0(:)             !WRF inverse density alpha
    real         , INTENT(IN)    :: zsg(:)             !model layer eta values (bottom-up)
    INTEGER      , INTENT(IN) :: nlvl          !number of model layers
    REAL         , INTENT(OUT) :: zlvls(:)  !model layer AGL values (m)
    REAL         , INTENT(OUT) :: zfull(:)  !model layer AGL values (m)
    REAL         , INTENT(IN) :: zsg_full(:)  !model layer thicknesses (eta)
    REAL         , INTENT(IN) :: dzsg_full(:)  !model interface level eta values (starting at top of bottom layer)

    integer :: k

    ! compute zagl at model layers/levels using WRF hydrostatic eqn: (same as in prfwrf)
    zlvls(1)=-mu*(zsg(1)-1.)*alt0(1)/grav
    zfull(1)=-mu*dzsg_full(1)*alt0(1)/grav
    do k=2,nlvl
       zlvls(k)=zlvls(k-1) - mu*0.5*(dzsg_full(k-1)*alt0(k-1)+dzsg_full(k)*alt0(k))/grav
       zfull(k)=zfull(k-1) - mu*dzsg_full(k)*alt0(k)/grav
    end do
    return
  end subroutine compute_allagl
  
end module module_wrfeta
