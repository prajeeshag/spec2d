module phys_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain
use mpp_domains_mod, only : mpp_get_global_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use constants_mod, only : PI

use time_manager_mod, only : time_type, set_time

use diag_manager_mod, only : diag_axis_init

use radiation_mod, only : init_radiation, con_solr, radiation, NF_ALBD

use sfc_mod, only : init_sfc, get_land_frac, get_albedo, get_tskin

use astronomy_mod, only : astronomy_init, diurnal_solar

implicit none
private

public :: init_phys, phys

type(domain2D), pointer :: domain
integer :: ntrac, nlev
integer :: is, ie, ilen, js, je, jlen
integer :: nlon, nlat

integer :: dt_rad
real :: deltim, deltimr
type(time_type) :: rad_dt

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg

logical :: initialized=.false.

namelist/phys_nml/dt_rad

contains

!--------------------------------------------------------------------------------   
subroutine init_phys(Time, deltim_in, domain_in, ntrac_in, nlev_in, &
                     lat_deg_in, lon_deg_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2d), target :: domain_in
    real, intent(in) :: deltim_in 
    integer, intent(in) :: ntrac_in, nlev_in
    real, intent(in) :: lat_deg_in(:), lon_deg_in(:)
    
    real, allocatable :: fland(:,:)
    integer :: axes(2), id_lev, unit, i, j
    integer :: jsg,jeg,isg,ieg

    unit = open_namelist_file()

    read(unit,nml=phys_nml)

    call close_file(unit)

    deltim = deltim_in

    deltimr = real(dt_rad)

    if (mod(deltimr,deltim)/=0.) call mpp_error('phys_mod', 'deltim for radiation ' &
                                        //'should be a multiple of model timestep', FATAL)

    rad_dt = set_time(dt_rad)

    domain => domain_in

    ntrac = ntrac_in
    
    nlev = nlev_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    call mpp_get_global_domain(domain,jsg,jeg,isg,ieg)

    jlen = je - js + 1
    ilen = ie - is + 1

    !axes(1) = diag_axis_init('lat',lat_deg_in(js:je),'degrees_N','Y',long_name='latitude',domain2=domain)
    axes(1) = diag_axis_init('lat',lat_deg_in(js:je),'degrees_N','Y', & 
              long_name='latitude',domain_decomp=[jsg,jeg,js,je])
    !axes(2) = diag_axis_init('lon',lon_deg_in(is:ie),'degrees_E','X',long_name='longitude',domain2=domain)
    axes(2) = diag_axis_init('lon',lon_deg_in(is:ie),'degrees_E','X', &
              long_name='longitude',domain_decomp=[isg,ieg,is,ie])

    call init_sfc(Time,deltim,domain,axes(1:2))

    allocate(fland(js:je,is:ie))

    call get_land_frac(fland)

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))

    do i = is, ie
        lat_deg(js:je,i) = lat_deg_in(js:je)
        lat_rad(js:je,i) = lat_deg_in(js:je)*PI/180.
    enddo

    do j = js, je
        lon_deg(j,is:ie) = lon_deg_in(is:ie)
        lon_rad(j,is:ie) = lon_deg_in(is:ie)*PI/180.
    enddo

    call astronomy_init()
    call init_radiation(Time, deltim, domain, ntrac, nlev, lat_deg_in, &
                        lon_deg_in, fland, axes, ind_q_in=1, ind_clw_in=2) 

    initialized = .true.

    deallocate(fland)

end subroutine init_phys


!--------------------------------------------------------------------------------   
subroutine phys(Time,tlyr,tr,p)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in) :: tlyr(1:nlev,js:je,is:ie)
    real, intent(in) :: tr(1:nlev,js:je,is:ie,1:ntrac)
    real, intent(in) :: p(js:je,is:ie) 

    real, dimension(js:je,is:ie) :: tskin, coszen, coszdg
    real, dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real :: solcon, rrsun

    coszen(:,:) = 0.0
    coszdg(:,:) = 0.0
    call diurnal_solar(lat_rad, lon_rad, Time, coszen, coszdg, rrsun, rad_dt)
    solcon = con_solr * rrsun
    coszdg = coszen * coszdg

    call get_tskin(Time,tskin)

    call get_albedo(Time,coszen,sfcalb)

    call radiation(Time, tlyr, tr, p, tskin, coszen, coszdg, sfcalb, solcon)

end subroutine phys

end module phys_mod
