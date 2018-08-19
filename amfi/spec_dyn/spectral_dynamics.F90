
module spectral_dynamics_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error, &
                    mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_max, &
                    mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                    mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather, &
                    mpp_declare_pelist, mpp_set_current_pelist, &
                    mpp_get_current_pelist, mpp_sum

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain, &
                            mpp_global_field, mpp_get_global_domain

use time_manager_mod, only : time_type

use diag_manager_mod, only : reg_df=>register_diag_field, send_data, diag_axis_init

use constants_mod, only : RVGAS, RDGAS, GRAV, RADIUS

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, &
                    stdlog, stdout, stderr
use fms_io_mod, only : fms_io_exit, restart_file_type, register_restart_field

use tracer_manager_mod, only : get_tracer_index, get_tracer_name, get_number_tracers

use field_manager_mod, only : MODEL_ATMOS

use transforms_mod, only : get_spherical_wave, get_lonsP, compute_ucos_vcos, compute_vor_div, &
                           spherical_to_grid, grid_to_spherical, init_transforms, get_latsF, &
                           register_spec_restart, save_spec_restart, get_latsP, &
                           restore_spec_restart 

use vertical_levels_mod, only: init_vertical_levels, get_ak_bk, get_vertical_vel, &
                               get_pressure_at_levels

use implicit_mod, only : init_implicit, do_implicit, do_implicit_adj

use gfidi_mod, only : gfidi_drv

use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

use ocpack_mod, only : oc_ny, oc_nx, oc_nfour, ocpack_typeP, npack=>oc_npack, get_ocpackP

implicit none
private

public :: init_spectral_dynamics, spectral_dynamics, get_latsP, get_lonsP, &
          finish_spectral_dynamics, save_spec_restart, restore_spec_restart

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

complex, dimension(:,:,:), allocatable :: sucos, svcos, stopo, stmp3d
complex, dimension(:,:,:,:), allocatable :: stmp3d1

real, allocatable, dimension(:,:,:) :: div, vor
real, allocatable, dimension(:,:) :: gtopo

real, allocatable, dimension(:) :: spdmax

real, allocatable, dimension(:,:) :: wtsbynlon

real, allocatable :: bfilt(:,:,:)

integer :: nlev
integer :: trunc
integer :: ntrac=0
integer :: nwaves_oe=0
    
integer :: isp, iep, ilenp, ocnx
integer :: jsp, jep, jlenp, ocny
type(ocpack_typeP) , allocatable :: ocpkP(:,:)
integer, allocatable :: pelist(:)

integer :: commID
real :: deltim=600.
real :: filta = 0.85 
logical :: mass_corr=.true.

integer :: i0 = -1

type(domain2d), pointer :: domain_g => NULL()

real, allocatable :: sin_latP(:,:), cosm2_latP(:,:), cosm_latP(:,:), wts_latP(:,:), &
                     cos_latP(:,:)
real, allocatable :: typdel(:)

integer, allocatable :: sph_wave(:,:)

real, parameter :: fv = RVGAS/RDGAS-1.

character(len=8) :: moist_tracer_names(10)
integer, allocatable :: moist_ind(:)
integer :: nmoist_tracers = 0

integer :: id_dtadv, id_dttot, id_dtdyn
integer, allocatable, dimension(:) :: id_dtradv, id_dtrtot, id_dtrdyn

character(len=8), parameter :: rou='sdyn'

contains

!--------------------------------------------------------------------------------
subroutine init_spectral_dynamics(Time, nlev_in, trunc_in, domain, deltim_in, rstrt, gaxis)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer, intent(in) :: nlev_in, trunc_in
    real, intent(in) :: deltim_in
    type(domain2d), target :: domain
    type(restart_file_type), intent(inout) :: rstrt
    integer, intent(out) :: gaxis(4)

    integer :: i, j, k, ntr
    real :: ref_temp
    real, dimension(nlev_in+1) :: ak, bk, si, plevp
    real, dimension(nlev_in) :: sl, plev
    real, allocatable :: tmpg(:,:)
    integer :: idx, tr, n, is, ie, nlon
    integer :: num_prog, num_diag
    
    call mpp_init()
    call fms_init()

    nlev = nlev_in 
    trunc = trunc_in 
    deltim = deltim_in 

    domain_g => domain

    moist_tracer_names(:) = ''
    moist_tracer_names(1:2) = ['sphum','clw']

    nmoist_tracers=count(moist_tracer_names/='')

    call get_number_tracers(MODEL_ATMOS,ntrac,num_prog,num_diag)

    if(mpp_pe()==mpp_root_pe()) then
        print *, ' num_tracers,num_prog,num_diag=', ntrac,num_prog,num_diag
        print *, 'nmoist_tracers =', nmoist_tracers
    endif
    
    allocate(moist_ind(nmoist_tracers))
    do n = 1, nmoist_tracers
        idx = get_tracer_index(MODEL_ATMOS, moist_tracer_names(n))
        if (idx<=0) call mpp_error(FATAL,'init_atmos: tracer not found: '// &
                          trim(moist_tracer_names(n)))
        moist_ind(n) = idx
    enddo

    ocnx = oc_nx()
    ocny = oc_ny()
    allocate(ocpkP(npack(),ocny))
    call get_ocpackP(ocpkP)

    call mpp_get_compute_domain(domain_g, jsp, jep, isp, iep)
    ilenp = iep - isp + 1
    jlenp = jep - jsp + 1
   
    allocate(pelist(mpp_npes())) 
    call mpp_get_current_pelist(pelist, commid=commID)

    call init_transforms(domain_g,trunc,nwaves_oe)
 
    call init_vertical_levels(nlev)
    
    allocate(sph_wave(nwaves_oe,2))
   
    allocate(tmpg(ocny,ocnx)) 
    allocate(sin_latP(jsp:jep,isp:iep), cosm2_latP(jsp:jep,isp:iep))
    allocate(cosm_latP(jsp:jep,isp:iep))
    allocate(wts_latP(jsp:jep,isp:iep), cos_latP(jsp:jep,isp:iep))
    allocate(wtsbynlon(jsp:jep,isp:iep))
    
    call get_latsP(sinlat = tmpg)
    sin_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)

    call get_latsP(cosm2lat = tmpg)
    cosm2_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)

    call get_latsP(cosmlat = tmpg)
    cosm_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)

    call get_latsP(wtslat = tmpg)
    wts_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)

    call get_latsP(coslat = tmpg)
    cos_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)

    deallocate(tmpg)

    do j = jsp, jep
        do i = isp, iep
            nlon = ocpkP(1,j)%ilen
            if (i > ocpkP(1,j)%ie) nlon = ocpkP(2,j)%ilen 
            wtsbynlon(j,i) = wts_latP(j,i)/nlon
        end do
    end do

    call get_ak_bk(ak_out=ak,bk_out=bk,sl_out=sl,si_out=si)
    allocate(typdel(nlev))
    typdel(:) = si(1:nlev) - si(2:nlev+1)
    
    !call get_pressure_at_levels(100.,prsi=plevp, prsl=plev)
    forall(k=1:nlev) plev(k) = k
    forall(k=1:nlev+1) plevp(k) = k
    
    gaxis(1) = diag_axis_init('ocy',[(float(i),i=jsp,jep)],'degrees_N','Y', &
                long_name='y', domain_decomp=[1,ocny,jsp,jep])

    gaxis(2) = diag_axis_init('ocx',[(float(i),i=isp,iep)],'degrees_E','X', &
                long_name='x', domain_decomp=[1,ocnx,isp,iep])

    gaxis(3) = diag_axis_init('lev',plev,'','Z',long_name='')

    gaxis(4) = diag_axis_init('levp',plevp,'','Z',long_name='')

    call init_diag_out(gaxis,Time) 

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
   
    call init_data()

    pdryini = 0.
    idx = register_restart_field(rstrt,'','pdryini',pdryini,domain_g,mandatory=.false.)

end subroutine init_spectral_dynamics


!--------------------------------------------------------------------------------   
subroutine init_data() 
!--------------------------------------------------------------------------------   
    integer :: i, idx, tr
    character(len=2) :: nm
    character(len=8) :: fldnm

    allocate(sucos(nlev,nwaves_oe,2))
    allocate(svcos(nlev,nwaves_oe,2))
    allocate(stmp3d(nlev,nwaves_oe,2))
    allocate(stmp3d1(nlev,nwaves_oe,2,ntrac))
    
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
        allocate(gatm(i)%u(nlev,jsp:jep,isp:iep))
        allocate(gatm(i)%v(nlev,jsp:jep,isp:iep))
        allocate(gatm(i)%tem(nlev,jsp:jep,isp:iep))
        allocate(gatm(i)%tr(nlev,jsp:jep,isp:iep,ntrac))
        allocate(gatm(i)%prs(1,jsp:jep,isp:iep))
    enddo
    
    allocate(dphi%u(nlev,jsp:jep,isp:iep))
    allocate(dphi%v(nlev,jsp:jep,isp:iep))
    allocate(dphi%tem(nlev,jsp:jep,isp:iep))
    allocate(dphi%tr(nlev,jsp:jep,isp:iep,ntrac))
    allocate(dphi%prs(1,jsp:jep,isp:iep))
    
    allocate(dlam%u(nlev,jsp:jep,isp:iep))
    allocate(dlam%v(nlev,jsp:jep,isp:iep))
    allocate(dlam%tem(nlev,jsp:jep,isp:iep))
    allocate(dlam%tr(nlev,jsp:jep,isp:iep,ntrac))
    allocate(dlam%prs(1,jsp:jep,isp:iep))
    
    allocate(dt%u(nlev,jsp:jep,isp:iep))
    allocate(dt%v(nlev,jsp:jep,isp:iep))
    allocate(dt%tem(nlev,jsp:jep,isp:iep))
    allocate(dt%tr(nlev,jsp:jep,isp:iep,ntrac))
    allocate(dt%prs(1,jsp:jep,isp:iep))
    
    allocate(div(nlev,jsp:jep,isp:iep))
    allocate(vor(nlev,jsp:jep,isp:iep))
    allocate(spdmax(nlev))

    spdmax = 0.

    idx=register_spec_restart('topo',stopo,.true.,0.)

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
    
    return
end subroutine init_data


!--------------------------------------------------------------------------------   
subroutine init_diag_out(axis, Time)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: axis(4)
    type(time_type), intent(in) :: Time
    integer :: n, axis1(3)
    character (len=256) :: fldnm, name, longname, units

    axis1 = [axis(3),axis(1),axis(2)]

    allocate(id_dtradv(ntrac))
    allocate(id_dtrtot(ntrac))
    allocate(id_dtrdyn(ntrac))

    id_dtadv = reg_df(rou, 'dtadv', axis1, Time, 'Temperature tendency (adv)', 'Ks-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'spectral_dynamics_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'adv'
        longname = trim(longname)//' tendency (advection)'
        units = trim(units)//' s-1'
        id_dtradv(n) = reg_df(rou, fldnm, axis1, Time, longname, units)
    enddo
    
    id_dttot = reg_df(rou, 'dttot', axis1, Time, 'Temperature tendency (total)', 'Ks-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'spectral_dynamics_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'tot'
        longname = trim(longname)//' tendency (total)'
        units = trim(units)//' s-1'
        id_dtrtot(n) = reg_df(rou, fldnm, axis1, Time, longname, units)
    enddo

    id_dtdyn = reg_df(rou, 'dtdyn', axis1, Time, 'Temperature tendency (dynamics)', 'Ks-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'spectral_dynamics_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'dyn'
        longname = trim(longname)//' tendency (dynamics)'
        units = trim(units)//' s-1'
        id_dtrdyn(n) = reg_df(rou, fldnm, axis1, Time, longname, units)
    enddo

    return
end subroutine init_diag_out


!--------------------------------------------------------------------------------
subroutine spectral_dynamics(Time,u,v,tem,tr,p,u1,v1,tem1,tr1,p1,vvel1)
!--------------------------------------------------------------------------------
    type(time_type), intent(in) :: Time
    real, intent(out), dimension(nlev,jsp:jep,isp:iep)       :: u, v, tem
    real, intent(out), dimension(nlev,jsp:jep,isp:iep)       :: u1, v1, tem1, vvel1
    real, intent(out), dimension(nlev,jsp:jep,isp:iep,ntrac) :: tr, tr1
    real, intent(out), dimension(jsp:jep,isp:iep)            :: p, p1

    real, dimension(nlev,jsp:jep,isp:iep) :: gtmp1
    complex, dimension(nlev,nwaves_oe,2) :: stmp1
    complex, dimension(nlev,nwaves_oe,2,ntrac) :: stmp2
    real :: val
    integer :: i, j, k, ntr
    logical :: used

    stmp3d = satm(2)%tem
    stmp3d1 = satm(2)%tr

    call compute_ucos_vcos(satm(2)%vor,satm(2)%div,sucos,svcos,do_trunc=.false.)
    
    call spherical_to_grid(satm(2)%prs,grid=gatm(1)%prs,lat_deriv=dphi%prs,lon_deriv=dlam%prs)

    call spherical_to_grid(satm(2)%tem,grid=gatm(1)%tem,lat_deriv=dphi%tem,lon_deriv=dlam%tem)

    call spherical_to_grid(satm(2)%div,grid=div)
    
    call spherical_to_grid(satm(2)%vor,grid=vor)
    
    call spherical_to_grid(sucos,grid=gatm(1)%u,lon_deriv=dlam%u)
    
    call spherical_to_grid(svcos,grid=gatm(1)%v,lon_deriv=dlam%v)
    
    do ntr = 1, ntrac
        call spherical_to_grid(satm(2)%tr(:,:,:,ntr),grid=gatm(1)%tr(:,:,:,ntr), &
            lat_deriv=dphi%tr(:,:,:,ntr),lon_deriv=dlam%tr(:,:,:,ntr))
    enddo
    
    do k = 1, nlev
        dphi%tem(k, jsp:jep, isp:iep) = dphi%tem(k, jsp:jep, isp:iep) * cosm2_latP(jsp:jep, isp:iep)
        dlam%tem(k, jsp:jep, isp:iep) = dlam%tem(k, jsp:jep, isp:iep) * cosm2_latP(jsp:jep, isp:iep)
        dlam%u(k, jsp:jep, isp:iep) = dlam%u(k, jsp:jep, isp:iep) * cosm2_latP(jsp:jep, isp:iep)
        dlam%v(k, jsp:jep, isp:iep) = dlam%v(k, jsp:jep, isp:iep) * cosm2_latP(jsp:jep, isp:iep)
    enddo

    do ntr = 1, ntrac
        do k = 1, nlev
            dphi%tr(k, jsp:jep, isp:iep, ntr) = dphi%tr(k, jsp:jep, isp:iep, ntr) * cosm2_latP(jsp:jep, isp:iep)
            dlam%tr(k, jsp:jep, isp:iep, ntr) = dlam%tr(k, jsp:jep, isp:iep, ntr) * cosm2_latP(jsp:jep, isp:iep)
        end do
    end do
    
    dphi%u = dlam%v - vor
    dphi%v = div - dlam%u
    spdmax = 0. 
    call gfidi_drv(nlev, ntrac, ilenp, jlenp, deltim, sin_latP(jsp:jep,isp:iep), cosm2_latP(jsp:jep,isp:iep), &
            div, gatm(1)%tem, gatm(1)%u, gatm(1)%v, gatm(1)%tr, dphi%prs, dlam%prs, gatm(1)%prs, &
            dphi%tem, dlam%tem, dphi%tr, dlam%tr, dlam%u, dlam%v, dphi%u, dphi%v, &
            dt%prs, dt%tem, dt%tr, dt%u, dt%v, spdmax)
  
    call mpi_max_arr(spdmax,size(spdmax))
    spdmax = sqrt(spdmax)

    call grid_to_spherical(dt%prs, satm(3)%prs, do_trunc=.true.)
    call grid_to_spherical(dt%tem, satm(3)%tem, do_trunc=.true.)
    call grid_to_spherical(dt%u, sucos, do_trunc=.false.)
    call grid_to_spherical(dt%v, svcos, do_trunc=.false.)
    
    do ntr = 1, ntrac
        call grid_to_spherical(dt%tr(:,:,:,ntr),satm(3)%tr(:,:,:,ntr),do_trunc=.true.)
    enddo
    
    call compute_vor_div(sucos,svcos,satm(3)%vor,satm(3)%div)
    
    satm(3)%vor = satm(1)%vor + 2.*deltim*satm(3)%vor

    do ntr = 1, ntrac
        if (id_dtradv(ntr)>0) then
            call spherical_to_grid(satm(3)%tr(:,:,:,ntr),grid=gtmp1)
            used = send_data(id_dtradv(ntr),gtmp1,Time)
        endif
    enddo

    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr

    do k = 1, size(satm(3)%div,1)
        satm(3)%div(k,:,:) = satm(3)%div(k,:,:) + stopo(1,:,:)
    enddo
  
    if (nwaves_oe>0) &  
        call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)
    
    if (id_dtadv>0) then
        call spherical_to_grid((satm(3)%tem-satm(1)%tem)/(2.*deltim),grid=gtmp1)
        used = send_data(id_dtadv,gtmp1,Time)
    endif

    call horiz_diffusion(satm(3)%tr,satm(3)%vor,satm(3)%div,satm(3)%tem,satm(1)%prs(1,:,:))
    
    if (id_dtdyn>0) then
        call spherical_to_grid((satm(3)%tem-satm(1)%tem)/(2.*deltim),grid=gtmp1)
        used = send_data(id_dtdyn,gtmp1,Time)
    endif

    do ntr = 1, ntrac
        if (id_dtrdyn(ntr)>0) then
            call spherical_to_grid((satm(3)%tr(:,:,:,ntr) &
                - satm(1)%tr(:,:,:,ntr))/(2.*deltim),grid=gtmp1)
            used = send_data(id_dtrdyn(ntr),gtmp1,Time)
        endif
    enddo

    call time_filter1()
    
    call compute_ucos_vcos(satm(3)%vor,satm(3)%div,sucos,svcos,do_trunc=.false.)
    call spherical_to_grid(sucos,grid=gatm(2)%u)
    call spherical_to_grid(svcos,grid=gatm(2)%v)
    call spherical_to_grid(satm(3)%tem,grid=gatm(2)%tem)
    call spherical_to_grid(satm(3)%prs,grid=gatm(2)%prs,lat_deriv=dphi%prs,lon_deriv=dlam%prs)

    do ntr = 1, ntrac
        call spherical_to_grid(satm(3)%tr(:,:,:,ntr),grid=gatm(2)%tr(:,:,:,ntr))
    enddo

    p(jsp:jep,isp:iep) = exp(gatm(1)%prs(1,jsp:jep,isp:iep))

    do k = 1, nlev 
        u(k,jsp:jep,isp:iep) = gatm(1)%u(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
        v(k,jsp:jep,isp:iep) = gatm(1)%v(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
    enddo

    tr(1:nlev,jsp:jep,isp:iep,1:ntrac) = gatm(1)%tr(1:nlev,jsp:jep,isp:iep,1:ntrac)
    where(tr<0.) tr=0.

    tem(1:nlev,jsp:jep,isp:iep) = gatm(1)%tem(1:nlev,jsp:jep,isp:iep) &
                                / (1.0+fv*tr(1:nlev,jsp:jep,isp:iep,1))
    
    p1(jsp:jep,isp:iep) = exp(gatm(2)%prs(1,jsp:jep,isp:iep))

    do k = 1, nlev 
        u1(k,jsp:jep,isp:iep) = gatm(2)%u(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
        v1(k,jsp:jep,isp:iep) = gatm(2)%v(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
    enddo

    call spherical_to_grid(satm(3)%div,grid=div)

    call get_vertical_vel(p1,dphi%prs(1,:,:),dlam%prs(1,:,:),div,u1,v1,vvel1)

    tr1(1:nlev,jsp:jep,isp:iep,1:ntrac) = gatm(2)%tr(1:nlev,jsp:jep,isp:iep,1:ntrac)
    where(tr1<0.) tr1=0.
    
    tem1(1:nlev,jsp:jep,isp:iep) = gatm(2)%tem(1:nlev,jsp:jep,isp:iep) &
                                 / (1.0+fv*tr1(1:nlev,jsp:jep,isp:iep,1))
    
end subroutine spectral_dynamics


!--------------------------------------------------------------------------------   
subroutine finish_spectral_dynamics(Time, tem, tr, u, v)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(nlev,jsp:jep,isp:iep)       :: u, v, tem
    real, intent(in), dimension(nlev,jsp:jep,isp:iep,ntrac) :: tr

    complex, dimension(nlev,nwaves_oe,2) :: rqt
    real, dimension(nlev,jsp:jep,isp:iep) :: gtmp1

    real :: pcorr, plvl1(nlev+1)
    integer :: k, j, ntr
    integer :: used

    call calc_mass_corr(exp(gatm(2)%prs(1,:,:)), tr, moist_ind, pcorr)

    do k = 1, size(u,1)
        gatm(2)%u(k,jsp:jep,isp:iep) = u(k,jsp:jep,isp:iep) * cos_latP(jsp:jep,isp:iep) !-> to ucos
        gatm(2)%v(k,jsp:jep,isp:iep) = v(k,jsp:jep,isp:iep) * cos_latP(jsp:jep,isp:iep) !-> to vcos
    enddo

    gatm(2)%tem(1:nlev,jsp:jep,isp:iep) = tem(1:nlev,jsp:jep,isp:iep) &
                                        * (1.0+fv*tr(1:nlev,jsp:jep,isp:iep,1)) !-> virtual temp

    call grid_to_spherical(gatm(2)%tem, satm(2)%tem, do_trunc=.true.)
    call grid_to_spherical(gatm(2)%u, sucos, do_trunc=.false.)
    call grid_to_spherical(gatm(2)%v, svcos, do_trunc=.false.)

    do ntr = 1, ntrac
        call grid_to_spherical(tr(1:nlev,jsp:jep,isp:iep,ntr), satm(2)%tr(:,:,:,ntr), do_trunc=.true.)
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

    if (id_dttot>0) then
        call spherical_to_grid((satm(2)%tem-stmp3d)/deltim,grid=gtmp1)
        used = send_data(id_dttot,gtmp1,Time)
    endif

    do ntr = 1, ntrac
        if (id_dtrtot(ntr)>0) then
            call spherical_to_grid((satm(2)%tr(:,:,:,ntr) &
                    - stmp3d1(:,:,:,ntr))/deltim,grid=gtmp1)
            used = send_data(id_dtrtot(ntr),gtmp1,Time)
        endif
    enddo

    return
end subroutine finish_spectral_dynamics


!--------------------------------------------------------------------------------   
subroutine calc_mass_corr(ps, trc, mi, pcorr)
!--------------------------------------------------------------------------------
    real, intent(in), dimension(jsp:jep,isp:iep) :: ps
    real, intent(in), dimension(nlev,jsp:jep,isp:iep,ntrac) :: trc
    integer, intent(in), dimension(:) :: mi
    real, intent(out) :: pcorr

    real, dimension(jsp:jep,isp:iep) :: pwat
    real, dimension(jsp:jep,isp:iep) :: pstmp
    real, dimension(nlev,jsp:jep,isp:iep) :: delp
    real, dimension(nlev+1,jsp:jep,isp:iep) :: plvl
    !real, dimension(ocny,ocnx) :: pwatg, psg
    !real, dimension(ocny) :: pwatl, psl
    real :: pwattot, pstot, pdryg

    integer :: levs

    levs = size(trc,1)

    call get_pressure_at_levels(ps,plvl)

    delp = plvl(1:levs,:,:) - plvl(2:levs+1,:,:)

    call calc_integral_moisture(delp, trc, pwat, mi)

    pwat = pwat * GRAV * 0.5 * 1.e-3 * wtsbynlon
    pwattot = sum(pwat)
    pstot = sum(ps * 0.5 * wtsbynlon)
    call mpp_sum(pwattot)
    call mpp_sum(pstot)

    !call mpp_global_field(domain_g, pwat, pwatg)
    !call mpp_global_field(domain_g, ps, psg)
    !pwatl = sum(pwatg,2) * GRAV * 0.5 * 1.e-3 / nlon
    !psl = sum(psg,2)        * 0.5         / nlon
    !pwattot = sum(pwatl*wts_lat)
    !pstot = sum(psl*wts_lat)

    pdryg = pstot - pwattot

    if (pdryini<=0.) pdryini = pdryg

    pcorr = (pdryini - pdryg)/pstot*sqrt(2.)

    return
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

    return
end subroutine init_bfiltr


subroutine damp_speed(dive,vore,teme,rte,ndexev,spdmax,jcap,delt)
    implicit none
    complex, intent(inout), dimension(:,:,:) :: dive, vore, teme
    complex, intent(inout), dimension(:,:,:,:) :: rte
    integer, intent(in), dimension(:,:) :: ndexev
    real, intent(in), dimension(:) :: spdmax
    real, intent(in) :: delt
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
    beta=RADIUS*cons1p009/delt     !constant
    alfadt=alfa*delt/RADIUS

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

subroutine mpi_max_arr(arr,n)
    include 'mpif.h'
    integer, intent(in) :: n
    real(kind=8), intent(inout) :: arr(n) 
    real(kind=8) :: buff(n) 
    integer :: ierr
 
    call MPI_ALLREDUCE(arr, buff, n, MPI_REAL8, MPI_MAX, commID, ierr)

    arr = buff

    return
end subroutine mpi_max_arr

end module spectral_dynamics_mod
