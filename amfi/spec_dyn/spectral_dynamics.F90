
module spectral_dynamics_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

use constants_mod, only : RVGAS, RDGAS

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init

use fms_io_mod, only : fms_io_exit 

use transforms_mod, only : get_spherical_wave
use transforms_mod, only : compute_ucos_vcos, compute_vor_div
use transforms_mod, only : spherical_to_grid, grid_to_spherical
use transforms_mod, only : init_transforms, get_lats, get_lons
use transforms_mod, only : register_spec_restart, restore_spec_state, save_spec_restart

use vertical_levels_mod, only: init_vertical_levels, get_ak_bk

use implicit_mod, only : init_implicit, do_implicit

use gfidi_mod, only : gfidi_drv

use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

implicit none
private

public :: init_spectral_dynamics, spectral_dynamics
public :: get_lats, get_lons

type satm_type
    complex, dimension(:,:,:),   allocatable :: vor
    complex, dimension(:,:,:),   allocatable :: div
    complex, dimension(:,:,:),   allocatable :: tem
    complex, dimension(:,:,:,:), allocatable :: tr
    complex, dimension(:,:,:),   allocatable :: prs
    integer :: ntrac = 0
end type satm_type

type(satm_type) :: satm(3)

type gatm_type
    real, dimension(:,:,:),   allocatable :: u
    real, dimension(:,:,:),   allocatable :: v
    real, dimension(:,:,:),   allocatable :: tem
    real, dimension(:,:,:,:), allocatable :: tr
    real, dimension(:,:,:),   allocatable :: prs
    integer :: ntrac = 0
end type gatm_type

type(gatm_type) :: gatm(2), dphi, dlam, dt

complex, dimension(:,:,:),   allocatable :: sucos, svcos
complex, dimension(:,:,:),   allocatable :: stopo

real, allocatable, dimension(:,:,:) :: div, vor
real, allocatable, dimension(:,:) :: gtopo

real, allocatable, dimension(:) :: spdmax

integer :: nlon, nlat, nlev
integer :: trunc
integer :: ntrac=0
  
integer :: nwaves_oe=0
    
integer :: isc, iec, ilen
integer :: jsc, jec, jlen
real :: deltim=600.
real :: filta = 0.85 

type(domain2d), pointer :: domain_g => NULL()

real, allocatable :: sin_lat(:), cosm2_lat(:), deg_lat(:), cosm_lat(:)

real, parameter :: fv = RVGAS/RDGAS-1.

contains

!--------------------------------------------------------------------------------
subroutine init_spectral_dynamics(nlon_in,nlat_in,nlev_in,trunc_in, &
                                  ntrac_in,domain,deltim_in)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: nlon_in, nlat_in, nlev_in, trunc_in, ntrac_in
    real, intent(in) :: deltim_in
    type(domain2d), target :: domain

    integer :: i, j, k, ntr
    character(len=8) :: fldnm
    real :: ref_temp
    integer, allocatable :: sph_wave(:,:)
    real, allocatable :: ak(:), bk(:), sl(:)
    character(len=2) :: nm
    integer :: idx, tr
    
    call mpp_init()
    call fms_init()

    nlon = nlon_in 
    nlat = nlat_in 
    nlev = nlev_in 
    trunc = trunc_in 
    ntrac = ntrac_in 
    deltim = deltim_in 

    domain_g => domain

    call mpp_get_compute_domain(domain_g, jsc, jec, isc, iec)
    
    ilen = iec-isc+1
    jlen = jec-jsc+1
    
    call init_transforms(domain_g,trunc,nwaves_oe)
    
    call init_vertical_levels(nlev)
    
    allocate(ak(nlev+1),bk(nlev+1),sl(nlev))
    
    allocate(sph_wave(nwaves_oe,2))
    
    allocate(sin_lat(jsc:jec),cosm2_lat(jsc:jec))
    allocate(deg_lat(jsc:jec), cosm_lat(jsc:jec))
    
    call get_lats(sinlat=sin_lat,cosm2lat=cosm2_lat, &
                  deglat=deg_lat,cosmlat=cosm_lat)

    if(mpp_pe()==mpp_root_pe()) print *, 'deg_lat=', deg_lat
    
    call get_ak_bk(ak_out=ak,bk_out=bk,sl_out=sl)
    call get_spherical_wave(sph_wave)
    
    ref_temp = 300.
    if (nlev>100) ref_temp=1500.
    
    call init_implicit(ak,bk,ref_temp,deltim,trunc)
    
    call init_horiz_diffusion(trunc,deltim,sl,sph_wave,bk)
    
    allocate(sucos(nlev,nwaves_oe,2))
    allocate(svcos(nlev,nwaves_oe,2))
    
    do i = 1, 3
        allocate(satm(i)%vor(nlev,nwaves_oe,2))
        allocate(satm(i)%div(nlev,nwaves_oe,2))
        allocate(satm(i)%tem(nlev,nwaves_oe,2))
        allocate(satm(i)%tr(nlev,nwaves_oe,2,ntrac))
        allocate(satm(i)%prs(1,nwaves_oe,2))
    enddo
    
    allocate(stopo(1,nwaves_oe,2))
    
    do i = 1, 2
        allocate(gatm(i)%u(nlev,jsc:jec,isc:iec))
        allocate(gatm(i)%v(nlev,jsc:jec,isc:iec))
        allocate(gatm(i)%tem(nlev,jsc:jec,isc:iec))
        allocate(gatm(i)%tr(nlev,jsc:jec,isc:iec,ntrac))
        allocate(gatm(i)%prs(1,jsc:jec,isc:iec))
    enddo
    
    allocate(dphi%u(nlev,jsc:jec,isc:iec))
    allocate(dphi%v(nlev,jsc:jec,isc:iec))
    allocate(dphi%tem(nlev,jsc:jec,isc:iec))
    allocate(dphi%tr(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dphi%prs(1,jsc:jec,isc:iec))
    
    allocate(dlam%u(nlev,jsc:jec,isc:iec))
    allocate(dlam%v(nlev,jsc:jec,isc:iec))
    allocate(dlam%tem(nlev,jsc:jec,isc:iec))
    allocate(dlam%tr(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dlam%prs(1,jsc:jec,isc:iec))
    
    allocate(dt%u(nlev,jsc:jec,isc:iec))
    allocate(dt%v(nlev,jsc:jec,isc:iec))
    allocate(dt%tem(nlev,jsc:jec,isc:iec))
    allocate(dt%tr(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dt%prs(1,jsc:jec,isc:iec))
    
    allocate(div(nlev,jsc:jec,isc:iec))
    allocate(vor(nlev,jsc:jec,isc:iec))
    allocate(spdmax(nlev))
    spdmax = 0.
    
    deallocate(ak,bk,sl)
    deallocate(sph_wave)

    idx=register_spec_restart('spec_res','topo',stopo,.false.,0.)
    do i = 1, 2
        nm='_m'
        if (i==2) nm='_n'
        idx=register_spec_restart('spec_res','vor'//nm,satm(i)%vor,.false.,0.)
        idx=register_spec_restart('spec_res','div'//nm,satm(i)%div,.false.,0.)
        idx=register_spec_restart('spec_res','tem'//nm,satm(i)%tem,.true.,273.15)
        idx=register_spec_restart('spec_res','prs'//nm,satm(i)%prs,.false.,0.)
        do tr = 1, ntrac
            write(fldnm,'(A,I3.3,A)') 'tr',tr,nm
            idx=register_spec_restart('spec_res',fldnm,satm(i)%tr(:,:,:,tr),.false.,0.)
        enddo
    enddo
            
    !call read_specdata('specdata','topo',stopo)
    !
    !call read_specdata('specdata','lnp_1',satm(1)%prs)
    !call read_specdata('specdata','lnp_2',satm(2)%prs)
    !
    !call read_specdata('specdata','vor_1',satm(1)%vor)
    !call read_specdata('specdata','vor_2',satm(2)%vor)
    !
    !call read_specdata('specdata','div_1',satm(1)%div)
    !call read_specdata('specdata','div_2',satm(2)%div)
    !
    !call read_specdata('specdata','tem_1',satm(1)%tem)
    !call read_specdata('specdata','tem_2',satm(2)%tem)
    !
    !do ntr = 1, ntrac
    !    write(fldnm,'(A,I1)') 'tr',ntr
    !    call read_specdata('specdata',trim(fldnm)//'_1',satm(1)%tr(:,:,:,ntr))
    !    call read_specdata('specdata',trim(fldnm)//'_2',satm(2)%tr(:,:,:,ntr))
    !enddo

    call restore_spec_state()
    call save_spec_restart('init')

end subroutine init_spectral_dynamics


!--------------------------------------------------------------------------------
subroutine spectral_dynamics(u,v,tem,tr,p,u1,v1,tem1,tr1,p1)
!--------------------------------------------------------------------------------
    real, intent(out), dimension(nlev,jsc:jec,isc:iec)       :: u, v, tem
    real, intent(out), dimension(nlev,jsc:jec,isc:iec)       :: u1, v1, tem1
    real, intent(out), dimension(nlev,jsc:jec,isc:iec,ntrac) :: tr, tr1
    real, intent(out), dimension(jsc:jec,isc:iec)            :: p, p1

    integer :: i, j, k, ntr
    character(len=8) :: fldnm

    call compute_ucos_vcos(satm(2)%vor,satm(2)%div,sucos,svcos,do_trunc=.false.)
    
    call spherical_to_grid(satm(2)%div,grid=div)
    
    call spherical_to_grid(satm(2)%vor,grid=vor)
    
    call spherical_to_grid(sucos,grid=gatm(1)%u,lon_deriv=dlam%u)
    
    call spherical_to_grid(svcos,grid=gatm(1)%v,lon_deriv=dlam%v)
    
    call spherical_to_grid(satm(2)%prs,grid=gatm(1)%prs,lat_deriv=dphi%prs,lon_deriv=dlam%prs)
    
    call spherical_to_grid(satm(2)%tem,grid=gatm(1)%tem,lat_deriv=dphi%tem,lon_deriv=dlam%tem)
    
    do ntr = 1, ntrac
        call spherical_to_grid(satm(2)%tr(:,:,:,ntr),grid=gatm(1)%tr(:,:,:,ntr), &
            lat_deriv=dphi%tr(:,:,:,ntr),lon_deriv=dlam%tr(:,:,:,ntr))
    enddo
    
    do j = jsc, jec
        dphi%tem(:,j,:) = dphi%tem(:,j,:) * cosm2_lat(j)
        dphi%tr(:,j,:,:) = dphi%tr(:,j,:,:) * cosm2_lat(j)
    
        dlam%tem(:,j,:) = dlam%tem(:,j,:) * cosm2_lat(j)
        dlam%tr(:,j,:,:) = dlam%tr(:,j,:,:) * cosm2_lat(j)
    
        dlam%u(:,j,:) = dlam%u(:,j,:) * cosm2_lat(j)
        dlam%v(:,j,:) = dlam%v(:,j,:) * cosm2_lat(j)
    enddo
    
    dphi%u = dlam%v - vor
    dphi%v = div - dlam%u
    
    call gfidi_drv(nlev, ntrac, ilen, jlen, deltim, sin_lat(jsc:jec), cosm2_lat(jsc:jec), &
            div, gatm(1)%tem, gatm(1)%u, gatm(1)%v, gatm(1)%tr, dphi%prs, dlam%prs, gatm(1)%prs, &
            dphi%tem, dlam%tem, dphi%tr, dlam%tr, dlam%u, dlam%v, dphi%u, dphi%v, &
            dt%prs, dt%tem, dt%tr, dt%u, dt%v, spdmax)
    
    call grid_to_spherical(dt%prs, satm(3)%prs, do_trunc=.true.)
    call grid_to_spherical(dt%tem, satm(3)%tem, do_trunc=.true.)
    call grid_to_spherical(dt%u, sucos, do_trunc=.false.)
    call grid_to_spherical(dt%v, svcos, do_trunc=.false.)
    
    do ntr = 1, ntrac
        call grid_to_spherical(dt%tr(:,:,:,ntr),satm(3)%tr(:,:,:,ntr),do_trunc=.true.)
    enddo
    
    call compute_vor_div(sucos,svcos,satm(3)%vor,satm(3)%div)
    
    satm(3)%vor = satm(1)%vor + 2.*deltim*satm(3)%vor
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr
    do k = 1, size(satm(3)%div,1)
        satm(3)%div(k,:,:) = satm(3)%div(k,:,:) + stopo(1,:,:)
    enddo
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr
    
    call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)
    
    call horiz_diffusion(satm(3)%tr,satm(3)%vor,satm(3)%div,satm(3)%tem,satm(1)%prs(1,:,:))
    
    call time_filter1()
    
    call compute_ucos_vcos(satm(3)%vor,satm(3)%div,sucos,svcos,do_trunc=.false.)
    call spherical_to_grid(sucos,grid=gatm(2)%u)
    call spherical_to_grid(svcos,grid=gatm(2)%v)
    call spherical_to_grid(satm(3)%tem,grid=gatm(2)%tem)
    call spherical_to_grid(satm(3)%prs,grid=gatm(2)%prs)
    do ntr = 1, ntrac
        call spherical_to_grid(satm(3)%tr(:,:,:,ntr),grid=gatm(2)%tr(:,:,:,ntr))
    enddo

    p(jsc:jec,isc:iec) = exp(gatm(1)%prs(1,jsc:jec,isc:iec))
    do j = jsc, jec
        u(1:nlev,j,isc:iec) = gatm(1)%u(1:nlev,j,isc:iec) * cosm_lat(j)
        v(1:nlev,j,isc:iec) = gatm(1)%v(1:nlev,j,isc:iec) * cosm_lat(j)
    enddo

    tr(1:nlev,jsc:jec,isc:iec,1:ntrac) = gatm(1)%tr(1:nlev,jsc:jec,isc:iec,1:ntrac)
    where(tr<0.) tr=0.

    tem(1:nlev,jsc:jec,isc:iec) = gatm(1)%tem(1:nlev,jsc:jec,isc:iec)/(1.0+fv*tr(:,:,:,1))
    
    p1(jsc:jec,isc:iec) = exp(gatm(2)%prs(1,jsc:jec,isc:iec))
    do j = jsc, jec
        u1(1:nlev,j,isc:iec) = gatm(2)%u(1:nlev,j,isc:iec) * cosm_lat(j)
        v1(1:nlev,j,isc:iec) = gatm(2)%v(1:nlev,j,isc:iec) * cosm_lat(j)
    enddo
    tr1(1:nlev,jsc:jec,isc:iec,1:ntrac) = gatm(2)%tr(1:nlev,jsc:jec,isc:iec,1:ntrac)
    where(tr1<0.) tr1=0.
    
    tem1(1:nlev,jsc:jec,isc:iec) = gatm(2)%tem(1:nlev,jsc:jec,isc:iec)/(1.0+fv*tr1(:,:,:,1))
    
end subroutine spectral_dynamics

!--------------------------------------------------------------------------------   
subroutine time_filter1()
!--------------------------------------------------------------------------------   
    real :: filtb

    filtb = (1.-filta)*0.5

    satm(1)%tem = filtb * satm(1)%tem + filta * satm(2)%tem
    satm(1)%div = filtb * satm(1)%div + filta * satm(2)%div
    satm(1)%vor = filtb * satm(1)%vor + filta * satm(2)%vor
    satm(1)%tr = filtb * satm(1)%tr + filta * satm(2)%tr

    satm(1)%prs = satm(2)%prs
    satm(2)%prs = satm(3)%prs

end subroutine time_filter1


end module spectral_dynamics_mod
