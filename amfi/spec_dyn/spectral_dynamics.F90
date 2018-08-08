
module spectral_dynamics_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_max
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain, &
                            mpp_global_field

use constants_mod, only : RVGAS, RDGAS, GRAV, RADIUS

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, &
                    stdlog, stdout, stderr
use fms_io_mod, only : fms_io_exit, restart_file_type, register_restart_field

use transforms_mod, only : get_spherical_wave
use transforms_mod, only : compute_ucos_vcos, compute_vor_div
use transforms_mod, only : spherical_to_grid, grid_to_spherical
use transforms_mod, only : init_transforms, get_lats, get_lons
use transforms_mod, only : register_spec_restart, save_spec_restart, restore_spec_restart 

use vertical_levels_mod, only: init_vertical_levels, get_ak_bk, get_vertical_vel
use vertical_levels_mod, only: get_pressure_at_levels

use implicit_mod, only : init_implicit, do_implicit, do_implicit_adj

use gfidi_mod, only : gfidi_drv

use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

implicit none
private

public :: init_spectral_dynamics, spectral_dynamics
public :: get_lats, get_lons, finish_spectral_dynamics
public :: save_spec_restart, restore_spec_restart

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

real :: pdryini

complex, dimension(:,:,:),   allocatable :: sucos, svcos
complex, dimension(:,:,:),   allocatable :: stopo

real, allocatable, dimension(:,:,:) :: div, vor
real, allocatable, dimension(:,:) :: gtopo

real, allocatable, dimension(:) :: spdmax

real, allocatable :: bfilt(:,:,:)

integer :: nlon, nlat, nlev
integer :: trunc
integer :: ntrac=0
integer, allocatable :: moist_ind(:)
  
integer :: nwaves_oe=0
    
integer :: isc, iec, ilen
integer :: jsc, jec, jlen
real :: deltim=600.
real :: filta = 0.85 
logical :: mass_corr=.true.

integer :: i0 = -1

type(domain2d), pointer :: domain_g => NULL()

real, allocatable :: sin_lat(:), cosm2_lat(:), deg_lat(:), cosm_lat(:), wts_lat(:), &
                     cos_lat(:)
real, allocatable :: typdel(:)

integer, allocatable :: sph_wave(:,:)

real, parameter :: fv = RVGAS/RDGAS-1.

contains

!--------------------------------------------------------------------------------
subroutine init_spectral_dynamics(nlon_in,nlat_in,nlev_in,trunc_in, &
                                  ntrac_in,domain,deltim_in,rstrt,mi)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: nlon_in, nlat_in, nlev_in, trunc_in, ntrac_in
    real, intent(in) :: deltim_in
    type(domain2d), target :: domain
    type(restart_file_type), intent(inout) :: rstrt
    integer, intent(in) :: mi(:) ! indices of moisture fields in tracer array

    integer :: i, j, k, ntr
    character(len=8) :: fldnm
    real :: ref_temp
    real, allocatable :: ak(:), bk(:), sl(:), si(:)
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

    allocate(moist_ind(size(mi)))
    moist_ind = mi

    call mpp_get_compute_domain(domain_g, jsc, jec, isc, iec)
    
    ilen = iec-isc+1
    jlen = jec-jsc+1
    
    call init_transforms(domain_g,trunc,nwaves_oe)
 
    call init_vertical_levels(nlev)
    
    allocate(ak(nlev+1),bk(nlev+1),sl(nlev),si(nlev+1),typdel(nlev))
    
    allocate(sph_wave(nwaves_oe,2))
    
    allocate(sin_lat(jsc:jec), cosm2_lat(jsc:jec))
    allocate(deg_lat(jsc:jec), cosm_lat(jsc:jec))
    allocate(wts_lat(jsc:jec), cos_lat(jsc:jec))
    
    call get_lats(sinlat=sin_lat,cosm2lat=cosm2_lat, &
                  deglat=deg_lat,cosmlat=cosm_lat,wtslat=wts_lat, &
                  coslat=cos_lat)

    call get_ak_bk(ak_out=ak,bk_out=bk,sl_out=sl,si_out=si)
    typdel(:) = si(1:nlev) - si(2:nlev+1)
     
    call get_spherical_wave(sph_wave)

    do i = 1, size(sph_wave,1)
        if (sph_wave(i,1) == 0) then
            i0 = .true.
            exit
        endif
    enddo

    call init_bfiltr(sph_wave,trunc)

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
        satm(i)%vor = 0.
        satm(i)%div = 0.
        satm(i)%tem = 0.
        satm(i)%tr = 0.
        satm(i)%prs = 0.
    enddo
    
    allocate(stopo(1,nwaves_oe,2))
    stopo = 0.
    
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

    idx=register_spec_restart('topo',stopo,.false.,0.)
    do i = 1, 2
        nm='_m'
        if (i==2) nm='_n'
        idx=register_spec_restart('vor'//nm,satm(i)%vor,.false.,0.)
        idx=register_spec_restart('div'//nm,satm(i)%div,.false.,0.)
        idx=register_spec_restart('tem'//nm,satm(i)%tem,.true.,0.)
        idx=register_spec_restart('prs'//nm,satm(i)%prs,.false.,0.)
        do tr = 1, ntrac
            write(fldnm,'(A,I3.3,A)') 'tr',tr,nm
            idx=register_spec_restart(fldnm,satm(i)%tr(:,:,:,tr),.false.,0.)
        enddo
    enddo
    
    pdryini = 0.
    idx = register_restart_field(rstrt,'','pdryini',pdryini,domain_g,mandatory=.false.)

end subroutine init_spectral_dynamics


!--------------------------------------------------------------------------------
subroutine spectral_dynamics(u,v,tem,tr,p,u1,v1,tem1,tr1,p1,vvel1)
!--------------------------------------------------------------------------------
    real, intent(out), dimension(nlev,jsc:jec,isc:iec)       :: u, v, tem
    real, intent(out), dimension(nlev,jsc:jec,isc:iec)       :: u1, v1, tem1, vvel1
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
    spdmax = 0. 
    call gfidi_drv(nlev, ntrac, ilen, jlen, deltim, sin_lat(jsc:jec), cosm2_lat(jsc:jec), &
            div, gatm(1)%tem, gatm(1)%u, gatm(1)%v, gatm(1)%tr, dphi%prs, dlam%prs, gatm(1)%prs, &
            dphi%tem, dlam%tem, dphi%tr, dlam%tr, dlam%u, dlam%v, dphi%u, dphi%v, &
            dt%prs, dt%tem, dt%tr, dt%u, dt%v, spdmax)
  
    do k = 1, size(spdmax) 
        call mpp_max(spdmax(k)) 
    enddo
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
  
    if (nwaves_oe>0) &  
        call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)
    
    call horiz_diffusion(satm(3)%tr,satm(3)%vor,satm(3)%div,satm(3)%tem,satm(1)%prs(1,:,:))
    
    call time_filter1()
    
    call compute_ucos_vcos(satm(3)%vor,satm(3)%div,sucos,svcos,do_trunc=.false.)
    call spherical_to_grid(sucos,grid=gatm(2)%u)
    call spherical_to_grid(svcos,grid=gatm(2)%v)
    call spherical_to_grid(satm(3)%tem,grid=gatm(2)%tem)
    call spherical_to_grid(satm(3)%prs,grid=gatm(2)%prs,lat_deriv=dphi%prs,lon_deriv=dlam%prs)

    
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

    call spherical_to_grid(satm(3)%div,grid=div)

    call get_vertical_vel(p1,dphi%prs(1,:,:),dlam%prs(1,:,:),div,u1,v1,vvel1)

    tr1(1:nlev,jsc:jec,isc:iec,1:ntrac) = gatm(2)%tr(1:nlev,jsc:jec,isc:iec,1:ntrac)
    where(tr1<0.) tr1=0.
    
    tem1(1:nlev,jsc:jec,isc:iec) = gatm(2)%tem(1:nlev,jsc:jec,isc:iec)/(1.0+fv*tr1(:,:,:,1))
    
end subroutine spectral_dynamics


!--------------------------------------------------------------------------------   
subroutine finish_spectral_dynamics(tem, tr, u, v)
!--------------------------------------------------------------------------------   
    real, intent(in), dimension(nlev,jsc:jec,isc:iec)       :: u, v, tem
    real, intent(in), dimension(nlev,jsc:jec,isc:iec,ntrac) :: tr

    complex, dimension(nlev,nwaves_oe,2) :: rqt

    real :: pcorr, plvl1(nlev+1)
    integer :: k, j, ntr

    call calc_mass_corr(exp(gatm(2)%prs(1,:,:)), tr, moist_ind, pcorr)

    do j = jsc, jec
        gatm(2)%u(1:nlev,j,isc:iec) = u(1:nlev,j,isc:iec) * cos_lat(j) !-> to ucos
        gatm(2)%v(1:nlev,j,isc:iec) = v(1:nlev,j,isc:iec) * cos_lat(j) !-> to vcos
    enddo

    gatm(2)%tem(1:nlev,jsc:jec,isc:iec) = tem(1:nlev,jsc:jec,isc:iec)*(1.0+fv*tr(:,:,:,1)) !-> virtual temp

    call grid_to_spherical(gatm(2)%tem, satm(2)%tem, do_trunc=.true.)
    call grid_to_spherical(gatm(2)%u, sucos, do_trunc=.false.)
    call grid_to_spherical(gatm(2)%v, svcos, do_trunc=.false.)
    do ntr = 1, ntrac
        call grid_to_spherical(gatm(2)%tr(:,:,:,ntr),satm(2)%tr(:,:,:,ntr),do_trunc=.true.)
    enddo
 
    call compute_vor_div(sucos,svcos,satm(2)%vor,satm(2)%div)

    satm(3)%vor = satm(3)%vor + (satm(2)%vor-satm(3)%vor) * bfilt

    rqt = (satm(2)%tr(:,:,:,1)-satm(3)%tr(:,:,:,1)) * bfilt 
    
    do ntr = 1, ntrac
        satm(3)%tr(:,:,:,ntr) = satm(3)%tr(:,:,:,ntr) + &
            (satm(2)%tr(:,:,:,ntr)-satm(3)%tr(:,:,:,ntr)) * bfilt
    enddo

    if (mass_corr) then
        satm(3)%prs = cmplx(0.,0.)
        if (i0>0) satm(3)%prs(1,0,1) = cmplx(pcorr,0.)
        
        do k = 1, nlev
            satm(3)%prs(1,:,:) = satm(3)%prs(1,:,:) + typdel(k) * rqt(k,:,:)
        enddo
  
        satm(2)%div = (satm(2)%div-satm(3)%div) * bfilt
        satm(2)%tem = (satm(2)%tem-satm(3)%tem) * bfilt

        if (nwaves_oe>0) &
            call do_implicit_adj(satm(3)%div, satm(3)%tem, satm(2)%prs, &
                             satm(2)%div, satm(2)%tem, satm(3)%prs, deltim)
    else
        call mpp_error(FATAL,'mass correction should be true.')
    endif

    call damp_speed(satm(3)%div, satm(3)%vor, satm(3)%tem, satm(3)%tr, &
                    sph_wave, spdmax, trunc, deltim)

    call time_filter2()

    return
end subroutine finish_spectral_dynamics


!--------------------------------------------------------------------------------   
subroutine calc_mass_corr(ps, trc, mi, pcorr)
!--------------------------------------------------------------------------------
    real, intent(in), dimension(:,:) :: ps
    real, intent(in), dimension(:,:,:,:) :: trc
    integer, intent(in), dimension(:) :: mi
    real, intent(out) :: pcorr

    real, dimension(size(ps,1),size(ps,2)) :: pwat
    real, dimension(size(trc,1),size(trc,2),size(trc,3)) :: delp
    real, dimension(size(trc,1)+1,size(trc,2),size(trc,3)) :: plvl
    real, dimension(nlat,nlon) :: pwatg, psg
    real, dimension(nlat) :: pwatl, psl
    real :: pwattot, pstot, pdryg

    integer :: levs

    levs = size(trc,1)

    call get_pressure_at_levels(ps,plvl)
    delp = plvl(1:levs,:,:) - plvl(2:levs+1,:,:)

    call calc_integral_moisture(delp, trc, pwat, mi)

    call mpp_global_field(domain_g, pwat, pwatg)

    call mpp_global_field(domain_g, ps, psg)

    pwatl = sum(pwatg,2) * GRAV * 0.5 * 1.e-3 / nlon
    psl = sum(psg,2)        * 0.5         / nlon

    pwattot = sum(pwatl*wts_lat)
    pstot = sum(psl*wts_lat)

    pdryg = pstot - pwattot

    if (pdryini<=0.) pdryini = pdryg

    pcorr = (pdryini - pdryg)/pstot*sqrt(2.)

end subroutine calc_mass_corr


!--------------------------------------------------------------------------------   
subroutine calc_integral_moisture(delp, trc, pwat, mi)
!--------------------------------------------------------------------------------   
    real, intent(in), dimension(:,:,:) :: delp
    real, intent(in), dimension(:,:,:,:) :: trc
    real, intent(out), dimension(:,:) :: pwat
    integer, intent(in), dimension(:) :: mi

    integer :: k, n, ni

    pwat = 0.

    do ni = 1, size(mi)
        n = mi(ni)
        do k = 1, size(trc,1)
            pwat = pwat + delp(k,:,:) * trc(k,:,:,n)
        enddo
    enddo

    return
end subroutine calc_integral_moisture

!--------------------------------------------------------------------------------   
subroutine time_filter2()
!--------------------------------------------------------------------------------   
    real :: filtb

    filtb = (1.-filta)*0.5

    satm(2)%tem = satm(3)%tem 
    satm(2)%div = satm(3)%div 
    satm(2)%vor = satm(3)%vor 
    satm(2)%tr  = satm(3)%tr  

    satm(1)%tem = satm(1)%tem + filtb * satm(2)%tem
    satm(1)%div = satm(1)%div + filtb * satm(2)%div
    satm(1)%vor = satm(1)%vor + filtb * satm(2)%vor
    satm(1)%tr  = satm(1)%tr  + filtb * satm(2)%tr

end subroutine time_filter2

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

!--------------------------------------------------------------------------------   
subroutine init_bfiltr(sph_wave, jcap)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: sph_wave(:,:)
    integer, intent(in) :: jcap

    real :: fd2
    integer :: nf0, nf1, i, n, nw, nwaves

    nwaves = size(sph_wave,1)

    allocate(bfilt(nlev,nwaves,2))

    nf0=(jcap+1)*2/3
    nf1=(jcap+1)
    fd2=1./(nf1-nf0)**2

    bfilt = 0.

    do i = 1, 2
        do n = 1, nwaves
            nw = sph_wave(n,i)
            bfilt(:,n,i) = max(1.-fd2*max(nw-nf0,0)**2,0.)
        enddo
    enddo

end subroutine init_bfiltr

subroutine damp_speed(dive,vore,teme,rte,ndexev,spdmax,jcap,deltim)
    implicit none
    complex, intent(inout), dimension(:,:,:) :: dive, vore, teme
    complex, intent(inout), dimension(:,:,:,:) :: rte
    integer, intent(in), dimension(:,:) :: ndexev
    real, intent(in), dimension(:) :: spdmax
    real, intent(in) :: deltim
    integer, intent(in) :: jcap

    integer :: it, k, nw, i
    integer :: nwaves, nlev, ntrac

    real :: alfa,alfadt,beta,coef,factor,rk,rncrit,sf,tk
    real :: cons0,cons1,cons1p009
    real :: cons2,cons2p5  

    nwaves = size(dive,2) 
    nlev = size(dive,1)
    ntrac = size(rte,4)

    cons0     = 0.d0        !constant
    cons1     = 1.d0        !constant
    cons1p009 = 1.009d0     !constant
    cons2     = 2.d0        !constant
    cons2p5   = 2.5d0       !constant

    alfa=cons2p5                    !constant
    beta=RADIUS*cons1p009/deltim     !constant
    alfadt=alfa*deltim/RADIUS

    do k = 1, nlev 
        rncrit=beta/spdmax(k)
        if (rncrit.lt.jcap) then
            coef=alfadt*spdmax(k)
            do nw = 1, 2
                do i = 1, nwaves
                    if (ndexev(i,nw).gt.rncrit) then
                        factor=1./(1.+((ndexev(i,nw)-rncrit)*coef))
                        dive(k,i,nw)=dive(k,i,nw)*factor
                        vore(k,i,nw)=vore(k,i,nw)*factor
                        teme(k,i,nw)=teme(k,i,nw)*factor
                        do it = 1, ntrac
                            rte(k,i,nw,it)=rte(k,i,nw,it)*factor
                        enddo
                    endif
                enddo
            enddo
        endif
    enddo

    return
end subroutine damp_speed

end module spectral_dynamics_mod
