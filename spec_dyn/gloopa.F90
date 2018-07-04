
!module spectral_dynamics_mod

program main

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init

use fms_io_mod, only : fms_io_exit 

use transforms_mod, only : read_specdata, get_spherical_wave
use transforms_mod, only : compute_ucos_vcos, compute_vor_div
use transforms_mod, only : spherical_to_grid, grid_to_spherical
use transforms_mod, only : init_transforms, get_lats

use vertical_levels_mod, only: init_vertical_levels, get_ak_bk

use implicit_mod, only : init_implicit, do_implicit

use gfidi_mod, only : gfidi_drv

use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

implicit none
!private

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
  
integer :: nwaves_oe=0
    
integer :: isc, iec, ilen, i, k
integer :: jsc, jec, jlen, j
integer :: unit, trunc
integer :: pe, ntrac=3, ntr
real :: deltim=600., ref_temp
real :: filta = 0.85 
character(len=16) :: rfile='gloopa', wfile='rgloopa'
character(len=8) :: fldnm
integer :: ishuff

type(domain2d) :: domain_g

real, allocatable :: ak(:), bk(:), sl(:)
real, allocatable :: sin_lat(:), cosm2_lat(:)

integer, allocatable :: sph_wave(:,:)

namelist/gloopa_nml/ trunc, nlon, nlat, nlev, deltim

call mpp_init()
call fms_init()

unit = open_namelist_file()
read(unit,nml=gloopa_nml)
call close_file(unit)

ishuff=1
if(mpp_npes()==1) ishuff=0

call mpp_define_domains( [1,nlat,1,nlon], [1,mpp_npes()], domain_g, kxy=1, ishuff=ishuff)
call mpp_get_compute_domain(domain_g, jsc, jec, isc, iec)
ilen = iec-isc+1
jlen = jec-jsc+1

call init_transforms(domain_g,trunc,nwaves_oe)
jlen = jec-jsc+1
ilen = iec-isc+1

call init_vertical_levels(nlev)

allocate(ak(nlev+1),bk(nlev+1),sl(nlev))

allocate(sph_wave(nwaves_oe,2))

allocate(sin_lat(nlat),cosm2_lat(nlat))

call get_lats(sinlat=sin_lat,cosm2lat=cosm2_lat)

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

call read_specdata('specdata','topo',stopo)

call read_specdata('specdata','lnp_1',satm(1)%prs)
call read_specdata('specdata','lnp_2',satm(2)%prs)

call read_specdata('specdata','vor_1',satm(1)%vor)
call read_specdata('specdata','vor_2',satm(2)%vor)

call read_specdata('specdata','div_1',satm(1)%div)
call read_specdata('specdata','div_2',satm(2)%div)

call read_specdata('specdata','tem_1',satm(1)%tem)
call read_specdata('specdata','tem_2',satm(2)%tem)

do ntr = 1, ntrac
    write(fldnm,'(A,I1)') 'tr',ntr
    call read_specdata('specdata',trim(fldnm)//'_1',satm(1)%tr(:,:,:,ntr))
    call read_specdata('specdata',trim(fldnm)//'_2',satm(2)%tr(:,:,:,ntr))
enddo

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

call write_data('rgloopa','div',div,domain_g)
call write_data('rgloopa','vor',vor,domain_g)

call write_data('rgloopa','dudt',dt%u,domain_g)
call write_data('rgloopa','dudphi',dphi%u,domain_g)
call write_data('rgloopa','dudlam',dlam%u,domain_g)
call write_data('rgloopa','u',gatm(1)%u,domain_g)

call write_data('rgloopa','dvdt',dt%v,domain_g)
call write_data('rgloopa','dvdphi',dphi%v,domain_g)
call write_data('rgloopa','dvdlam',dlam%v,domain_g)
call write_data('rgloopa','v',gatm(1)%v,domain_g)

call write_data('rgloopa','dpdt',dt%prs,domain_g)
call write_data('rgloopa','dpdphi',dphi%prs,domain_g)
call write_data('rgloopa','dpdlam',dlam%prs,domain_g)
call write_data('rgloopa','p',gatm(1)%prs,domain_g)

call write_data('rgloopa','dtemdt',dt%tem,domain_g)
call write_data('rgloopa','dtemdphi',dphi%tem,domain_g)
call write_data('rgloopa','dtemdlam',dlam%tem,domain_g)
call write_data('rgloopa','tem',gatm(1)%tem,domain_g)

do ntr = 1, ntrac
    write(fldnm,'(A,I2.2)') 'tr',ntr
    call write_data('rgloopa','d'//trim(fldnm)//'dt',dt%tr(:,:,:,ntr),domain_g)
    call write_data('rgloopa','d'//trim(fldnm)//'dphi',dphi%tr(:,:,:,ntr),domain_g)
    call write_data('rgloopa','d'//trim(fldnm)//'dlam',dlam%tr(:,:,:,ntr),domain_g)
    call write_data('rgloopa',trim(fldnm),gatm(1)%tr(:,:,:,ntr),domain_g)
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
satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr

call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                 satm(2)%div, satm(2)%tem, satm(2)%prs, &
                 satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)

call horiz_diffusion(satm(3)%tr,satm(3)%vor,satm(3)%div,satm(3)%tem,satm(1)%prs(1,:,:))

call time_filter1()

call spherical_to_grid(satm(3)%div,grid=div)
call spherical_to_grid(satm(3)%tem,grid=gatm(2)%tem)
call spherical_to_grid(satm(3)%prs,grid=gatm(2)%prs)

call write_data('rgloopa','divdt',div,domain_g)
call write_data('rgloopa','temdt',gatm(2)%tem,domain_g)
call write_data('rgloopa','prsdt',gatm(2)%prs,domain_g)

call fms_io_exit()
call mpp_exit()

contains

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


end program main
