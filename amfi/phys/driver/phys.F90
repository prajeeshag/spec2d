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

use time_manager_mod, only : time_type, set_time, operator(==), operator(+), assignment(=)

use diag_manager_mod, only : diag_axis_init, reg_df=>register_diag_field, send_data

use radiation_mod, only : init_radiation, con_solr, radiation, NF_ALBD

use sfc_mod, only : init_sfc, get_land_frac, set_surface

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
type(time_type) :: time_step_rad, rad_time, time_step

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg
real, allocatable, dimension(:,:) :: rsds, rsus, rsns
real, allocatable, dimension(:,:,:) :: htsw
real, allocatable, dimension(:,:) :: rlds, rlus
real, allocatable, dimension(:,:,:) :: htlw

integer :: id_rsds, id_rsus, id_rsns, id_htsw

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
    
    real, allocatable :: fland(:,:), ilevs(:), ilevsp(:)
    integer :: axes(4), id_lev, unit, i, j
    integer :: jsg,jeg,isg,ieg
    character (len=8) :: rou='am_phys'

    unit = open_namelist_file()

    read(unit,nml=phys_nml)

    call close_file(unit)

    deltim = deltim_in

    deltimr = real(dt_rad)

    if (mod(deltimr,deltim)/=0.) call mpp_error('phys_mod', 'deltim for radiation ' &
                                        //'should be a multiple of model timestep', FATAL)

    time_step_rad = set_time(dt_rad)
    time_step = set_time(int(dt_rad))
    rad_time = Time

    domain => domain_in

    ntrac = ntrac_in
    
    nlev = nlev_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    call mpp_get_global_domain(domain,jsg,jeg,isg,ieg)

    jlen = je - js + 1
    ilen = ie - is + 1

    allocate(ilevs(nlev))
    allocate(ilevsp(nlev+1))

    forall(i=1:nlev) ilevs(i) = i
    forall(i=1:nlev+1) ilevsp(i) = i

    axes(1) = diag_axis_init('lat',lat_deg_in(js:je),'degrees_N','Y', & 
              long_name='latitude',domain_decomp=[jsg,jeg,js,je])
    axes(2) = diag_axis_init('lon',lon_deg_in(is:ie),'degrees_E','X', &
              long_name='longitude',domain_decomp=[isg,ieg,is,ie])
    axes(3) = diag_axis_init('lev',ilevs,'','Z',long_name='')
    axes(4) = diag_axis_init('levp',ilevsp,'','Z',long_name='')

    deallocate(ilevs,ilevsp)

    call init_sfc(Time,deltim,domain,axes(1:2))

    allocate(fland(js:je,is:ie))

    call get_land_frac(fland)

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))
    allocate(htsw(nlev,js:je,is:ie))
    allocate(rsds(js:je,is:ie))
    allocate(rsus(js:je,is:ie))
    allocate(rsns(js:je,is:ie))
    allocate(htlw(nlev,js:je,is:ie))
    allocate(rlds(js:je,is:ie))
    allocate(rlus(js:je,is:ie))

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


    deallocate(fland)

    ! Diag out
    ! -------------------------------------------------------------------------------- 
    
    id_rsds = reg_df(rou, 'rsds', axes(1:2), Time, 'Surface Downwelling Shortwave',  'W m-2')

    id_rsus = reg_df(rou, 'rsus', axes(1:2), Time, 'Surface Upwelling Shortwave',  'W m-2')

    id_rsns = reg_df(rou, 'rsns', axes(1:2), Time, 'Surface Net Shortwave',  'W m-2')

    id_htsw = reg_df(rou, 'htsw', [axes(3),axes(1),axes(2)], Time, 'Shortwave Heating Rate',  '?')

    !--------------------------------------------------------------------------------    


    initialized = .true.


end subroutine init_phys


!--------------------------------------------------------------------------------   
subroutine phys(Time,tlyr,tr,p)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in) :: tlyr(1:nlev,js:je,is:ie)
    real, intent(in) :: tr(1:nlev,js:je,is:ie,1:ntrac)
    real, intent(in) :: p(js:je,is:ie) 

    real, dimension(js:je,is:ie) :: tskin, fracday, coszen, rcoszen
    real, dimension(js:je,is:ie) :: rsdsz, rsusz, rsnsz
    real, dimension(nlev,js:je,is:ie) :: htswz
    real, dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real, dimension(js:je,is:ie) :: sfcemis
    real :: solcon, rrsun
    integer :: k
    logical :: used


    if (rad_time==Time) then
        coszen(:,:) = 0.0
        rcoszen(:,:) = 0.0
        call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step_rad)
        call set_surface(Time,tskin,coszen,sfcalb,sfcemis)
        solcon = con_solr * rrsun
        where(coszen>0.) rcoszen = 1./coszen 
        call radiation(Time, tlyr, tr, p, tskin, coszen, fracday, sfcalb, sfcemis, &
                       solcon, htsw, rsds, rsus, htlw, rlds, rlus)
        do k = 1, size(htsw,1)
            htsw(k,:,:) = htsw(k,:,:) * rcoszen
        enddo
        rsds = rsds * rcoszen
        rsus = rsus * rcoszen
        rad_time = rad_time + time_step_rad
    else
        call set_surface(Time,tskin)
    endif

    call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step)
    coszen = coszen * fracday
    do k = 1, size(htsw,1)
        htswz(k,:,:) = htsw(k,:,:) * coszen
    enddo
    rsdsz = rsds * coszen
    rsusz = rsus * coszen
    rsnsz = rsdsz - rsusz


    used = send_data(id_rsds, rsdsz, Time) 
    used = send_data(id_rsus, rsusz, Time) 
    used = send_data(id_rsns, rsnsz, Time) 
    used = send_data(id_htsw, htswz, Time) 

end subroutine phys

end module phys_mod
