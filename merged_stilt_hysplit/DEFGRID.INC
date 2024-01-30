! DEFGRID.INC 

! Defines meteorological grid location, variables contained during each time, 
! and file structure. Required as an INCLUDE in most MET subroutines.

! Last Revised: 31 Mar 1998 (RRD)
!               15 Apr 1999 (RRD) - added surface height flag
!               01 May 2000 (RRD) - scaling variable for input data
!               28 Jul 2000 (RRD) - additional surface flux flags
!               14 Aug 2000 (RRD) - zflg for pressure data
!               20 Sep 2000 (RRD) - fortran90 upgrade
!               09 Mar 2001 (RRD) - global grid definitions
!               29 Aug 2001 (RRD) - time periods added to array structure
!               09 Feb 2002 (RRD) - global subgrid switch
!               26 Feb 2002 (RRD) - downward shortwave flux flag
!               21 Jun 2002 (RRD) - maximum levels increased
!               17 Oct 2003 (RRD) - turbulent kinetic energy flag
!               10 Nov 2003 (RRD) - velocity variance flag
!               14 Apr 2006 (RRD) - wrf options in drec
!               22 May 2006 (RRD) - added mixed layer depth flag
!               10 Aug 2006 (RRD) - prime meridian test
!               22 Mar 2007 (RRD) - added zero element to time array
!               15 Aug 2008 (RRD) - split horizontal and vertical mixing

!-------------------------------------------------------------------------------
! meteorological data file limits
!-------------------------------------------------------------------------------

!      USE map_utils 

  USE map_utils,only : proj_info
  real, parameter :: vmiss=-999.99, vmissle=-999.


  INTEGER,  PARAMETER  :: mlvl     = 90    ! maximum levels  
  INTEGER,  PARAMETER  :: mvar     = 35    ! maximum variables

  INTEGER,  PARAMETER  :: mgrd     = 60    ! maximum grids
  INTEGER,  PARAMETER  :: mtim     = 60    ! maximum time periods

!-------------------------------------------------------------------------------
! define the geographic GRID parameters
!-------------------------------------------------------------------------------

  TYPE GSET
     SEQUENCE

     REAL*8         :: gbase(15)   ! array that defines grid conversions
     CHARACTER(4)   :: model_id    ! model data identification
     INTEGER        :: number      ! grid identification number
     REAL           :: pole_lat    ! pole latitude of projection
     REAL           :: pole_lon    ! pole longitude of projection
     REAL           :: ref_lat     ! reference latitude (size valid)
     REAL           :: ref_lon     ! reference longitude 
     REAL           :: size        ! grid size at ref lat lon
     REAL           :: orient      ! orientation
     REAL           :: tang_lat    ! tangent latitude
     REAL           :: sync_xp     ! synchronization point in x
     REAL           :: sync_yp     ! synchronization point in y
     REAL           :: sync_lat    ! latitude of sync point 
     REAL           :: sync_lon    ! longitude of sync point 
     REAL           :: dummy       ! reserved for future use 
     INTEGER        :: nx          ! number of points west to east
     INTEGER        :: ny          ! number of points south to north
     INTEGER        :: nz          ! number of levels 
     LOGICAL        :: latlon      ! flag if input grid is latlon
     LOGICAL        :: prime       ! lat/lon grid crosses the prime meridian
     CHARACTER(2)   :: gbldat      ! flag if latlon input grid is global
     CHARACTER(2)   :: global      ! flag if latlon subgrid is global
     LOGICAL        :: datload     ! flag to load data into subgrid
     INTEGER        :: lx1,ly1     ! potential sub-grid corner position
     INTEGER        :: lxr,lyr     ! potential sub-grid size or range
     REAL           :: lxc,lyc     ! potential sub-grid center position

!dwen(20090818) *********
     type (proj_info) :: proj        !WPS mapping routines data structure
     real,dimension(:,:), pointer :: mapfactor !mapfactor array for WPS 
!*************************

  END TYPE

!-------------------------------------------------------------------------------
! define structure of Data RECords within a file
!-------------------------------------------------------------------------------

  TYPE MSET
     SEQUENCE

     INTEGER        :: type                  ! data type flag 1:old 2:new
     INTEGER        :: z_flag                ! 1:sigma 2:pressure 3:terrain
     REAL           :: zmdlt                 ! scaling variable for input data
     REAL           :: height         (mlvl) ! height of each data level
     INTEGER        :: num_varb       (mlvl) ! variables per data level
     CHARACTER(4)   :: varb_id   (mvar,mlvl) ! character identification
     INTEGER        :: chk_sum   (mvar,mlvl) ! checksum for each field
     INTEGER        :: rec_per               ! records per time period
     INTEGER        :: delta                 ! minutes between times
     INTEGER        :: offset                ! number of records offset
     INTEGER        :: acycle                ! accumulation cycle (minutes)
     INTEGER        :: kdef                  ! deformation horizontal mixing 
     INTEGER        :: kzmix                 ! vertical profile averaging   
     INTEGER        :: kbls                  ! pbl stability derived from
     INTEGER        :: kblt                  ! pbl turbulence method 
     INTEGER        :: wrf                   ! 0:no 1:arw 2:nmm 3:stilt
     REAL           :: tvmix                 ! trop vert mix scale factor
     LOGICAL        :: dzdt                  ! vertical motion flag
     LOGICAL        :: uflg                  ! surface wind flag
     LOGICAL        :: tflg                  ! surface temp flag
     LOGICAL        :: mixd                  ! mixed layer depth available
     LOGICAL        :: dswf                  ! downward shortwave flux flag
     LOGICAL        :: eflx                  ! exchange coefficient flag
     LOGICAL        :: hflx                  ! sensible heat flux flag
     LOGICAL        :: uflx                  ! momentum flux flag
     LOGICAL        :: ustr                  ! friction velocity flag
     LOGICAL        :: tstr                  ! friction temperature flag
     LOGICAL        :: velv                  ! velocity variance flag
     LOGICAL        :: tken                  ! TKE field =0.5*(u'2+v'2+w'2)
     LOGICAL        :: shgt                  ! surface terrain height flag
     LOGICAL        :: prss                  ! surface pressure flag
     LOGICAL        :: zflg                  ! upper level pressure
     LOGICAL        :: qflg                  ! humidity expressed as specific
     logical        :: tkef
     LOGICAL        :: rflg           (mlvl) ! upper level humidity
     LOGICAL        :: wflg           (mlvl) ! upper level w-velocity

!dwen(20090818) ********
! CHG:(11/20/01) convective precipitation flag
     LOGICAL        :: cflg
! CHG:(11/20/01) total cloud flag
     LOGICAL        :: tclf
! CHG:(12/04/01) low cloud flag
     LOGICAL        :: lclf
! CHG:(11/20/01) shortw. downw. radiation flag
     LOGICAL        :: radf
! CHG:(22/01/03) soil moisture flag
     LOGICAL        :: slmf
! *******************************
     REAL :: wrfptop
     REAL, ALLOCATABLE :: wrfvcoords(:,:)
     
  END TYPE

!-------------------------------------------------------------------------------
! defines a time structure
!-------------------------------------------------------------------------------

  TYPE TSET
     SEQUENCE

     INTEGER        :: yr          ! year
     INTEGER        :: mo          ! month
     INTEGER        :: da          ! day
     INTEGER        :: hr          ! hour
     INTEGER        :: mn          ! minute
     INTEGER        :: ic          ! forecast hour
     INTEGER        :: macc        ! minutes accumulated

  END TYPE

!-------------------------------------------------------------------------------
! defines structure of meteorological data FILE
!-------------------------------------------------------------------------------

  TYPE FSET
     SEQUENCE

     CHARACTER(80)  :: dir         ! directory of meteorology
     CHARACTER(80)  :: meteo       ! name of meteorology file
     INTEGER        :: kunit       ! unit number assigned
     INTEGER        :: rec_len     ! record length (bytes)
     INTEGER        :: endrec      ! number of last record in file
     TYPE(tset)     :: first       ! beginning time period of file
     TYPE(tset)     :: last        ! ending time period of file
!dwen(20090818)**************
     integer        :: period      ! pointer to time period 1 or 2 
!****************************

  END TYPE

!-------------------------------------------------------------------------------
! Define file characteristics for each meteorological grid
!-------------------------------------------------------------------------------

  TYPE(gset)        :: grid  (0:mgrd,0:mtim)
  TYPE(mset)        :: drec  (0:mgrd,0:mtim) 
  TYPE(fset)        :: file  (0:mgrd,0:mtim) 

  ! Overload assignment operator with deep copy for mset:
  interface assignment (=)
     module procedure copy_mset
     module procedure copy_mset_arr
  end interface

contains

  subroutine copy_mset(mset_out, mset_in)

    ! Does a deep copy of mset derived type
    ! Needed because of the allocatable array mset_in/_out%wrfvcoords
    
    implicit none

    type(mset), intent(in) :: mset_in
    type(mset), intent(inout) :: mset_out

    mset_out%type=mset_in%type
    mset_out%z_flag=mset_in%z_flag
    mset_out%zmdlt=mset_in%zmdlt
    mset_out%height(:)=mset_in%height(:)
    mset_out%num_varb(:)=mset_in%num_varb(:)
    mset_out%varb_id(:,:)=mset_in%varb_id(:,:)
    mset_out%chk_sum(:,:)=mset_in%chk_sum(:,:)
    mset_out%rec_per=mset_in%rec_per
    mset_out%delta=mset_in%delta
    mset_out%offset=mset_in%offset
    mset_out%acycle=mset_in%acycle
    mset_out%kdef=mset_in%kdef
    mset_out%kzmix=mset_in%kzmix
    mset_out%kbls=mset_in%kbls
    mset_out%kblt=mset_in%kblt
    mset_out%wrf=mset_in%wrf
    mset_out%tvmix=mset_in%tvmix
    mset_out%dzdt=mset_in%dzdt
    mset_out%uflg=mset_in%uflg
    mset_out%tflg=mset_in%tflg
    mset_out%mixd=mset_in%mixd
    mset_out%dswf=mset_in%dswf
    mset_out%eflx=mset_in%eflx
    mset_out%hflx=mset_in%hflx
    mset_out%uflx=mset_in%uflx
    mset_out%ustr=mset_in%ustr
    mset_out%tstr=mset_in%tstr
    mset_out%velv=mset_in%velv
    mset_out%tken=mset_in%tken
    mset_out%shgt=mset_in%shgt
    mset_out%prss=mset_in%prss
    mset_out%zflg=mset_in%zflg
    mset_out%qflg=mset_in%qflg
    mset_out%tkef=mset_in%tkef
    mset_out%rflg(:)=mset_in%rflg(:)
    mset_out%wflg(:)=mset_in%wflg(:)
    mset_out%cflg=mset_in%cflg
    mset_out%tclf=mset_in%tclf
    mset_out%lclf=mset_in%lclf
    mset_out%radf=mset_in%radf
    mset_out%slmf=mset_in%slmf
    mset_out%wrfptop=mset_in%wrfptop
    if (allocated(mset_out%wrfvcoords)) deallocate(mset_out%wrfvcoords)
    if (allocated(mset_in%wrfvcoords)) then
       allocate(mset_out%wrfvcoords(size(mset_in%wrfvcoords,1),size(mset_in%wrfvcoords,2)))
       mset_out%wrfvcoords(:,:)=mset_in%wrfvcoords(:,:)
    end if
  end subroutine copy_mset

  subroutine copy_mset_arr(mset_arr_out, mset_arr_in)

    ! Does a deep copy of an array of mset derived type
    ! Needed because of the allocatable array mset_in/_out%wrfvcoords
    
    implicit none

    type(mset), intent(in) :: mset_arr_in(:)
    type(mset), intent(inout) :: mset_arr_out(:)

    integer :: i, n

    n=size(mset_arr_in)
    do i=1,n
       call copy_mset(mset_arr_out(i),mset_arr_in(i))
    end do

  end subroutine copy_mset_arr
