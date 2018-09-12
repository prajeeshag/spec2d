
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

use diag_manager_mod, only : reg_df=>register_diag_field, send_data, diag_axis_init, register_static_field

use constants_mod, only : RVGAS, RDGAS, GRAV, RADIUS

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, &
                    stdlog, stdout, stderr
use fms_io_mod, only : fms_io_exit, restart_file_type, register_restart_field, &
        file_exist, field_exist, field_size

use tracer_manager_mod, only : get_tracer_index, get_tracer_name, get_number_tracers

use field_manager_mod, only : MODEL_ATMOS

use transforms_mod, only : get_spherical_wave, get_lonsP, compute_ucos_vcos, compute_vor_div, &
                           spherical_to_grid, grid_to_spherical, init_transforms, get_latsF, &
                           register_spec_restart, save_spec_restart, get_latsP, &
                           restore_spec_restart, end_transforms, save_wisdom

use vertical_levels_mod, only: init_vertical_levels, get_ak_bk, get_vertical_vel, &
                               get_pressure_at_levels

use implicit_mod, only : init_implicit, do_implicit, do_implicit_adj

use gfidi_mod, only : gfidi_drv

use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

use ocpack_mod, only : oc_ny, oc_nx, oc_nfour, ocpack_typeP, npack=>oc_npack, get_ocpackP

use spec_comm_mod, only : spec_comm_max, spec_comm_min

#ifdef AQUAPLANET
use aqua_planet_mod, only : aquape_init_temp
#endif

implicit none
private

public :: init_spectral_dynamics, spectral_dynamics, get_latsP, get_lonsP, &
          finish_spectral_dynamics, save_spec_restart, end_spectral_dynamics

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

real :: pdryini = 98.2633

complex, dimension(:,:,:), allocatable :: sucos, svcos, stopo, soro, stmp3d
complex, dimension(:,:,:,:), allocatable :: stmp3d1

real, allocatable, dimension(:,:,:) :: div, vor
real, allocatable, dimension(:,:,:) :: gtopo

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
logical :: zero_topo=.false.
integer :: niter_topo = 10

integer :: i0 = -1

type(domain2d), pointer :: domain_g => NULL()

real, allocatable :: sin_latP(:,:), cosm2_latP(:,:), cosm_latP(:,:), wts_latP(:,:), &
                     cos_latP(:,:), deg_latP(:,:)
real, allocatable :: typdel(:)

integer, allocatable :: sph_wave(:,:), nnp1(:,:)

character(len=8) :: moist_tracer_names(10)
integer, allocatable :: moist_ind(:)
integer :: nmoist_tracers = 0, twotimes=0

integer :: id_dtadv, id_dttot, id_dtdyn, id_topo, id_prs_ini, id_tem_ini
integer, allocatable, dimension(:) :: id_dtradv, id_dtrtot, id_dtrdyn

real, parameter :: fv = RVGAS/RDGAS-1., GA2=GRAV/(RADIUS*RADIUS)

character(len=8), parameter :: rou='sdyn'
character(len=256) :: topo_file='INPUT/topography.nc', topo_field='topo'
character(len=256) :: temp_file='INPUT/temp_pres.nc'

contains

!--------------------------------------------------------------------------------
subroutine init_spectral_dynamics(Time, nlev_in, trunc_in, domain, deltim_in, gaxis)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer, intent(in) :: nlev_in, trunc_in
    real, intent(in) :: deltim_in
    type(domain2d), target :: domain
    integer, intent(out) :: gaxis(4)

    integer :: i, j, k, ntr
    real :: ref_temp, tmpmx(1)
    real, dimension(nlev_in+1) :: ak, bk, si, plevp
    real, dimension(nlev_in) :: sl, plev
    real, allocatable :: tmpg(:,:)
    integer :: idx, tr, n, is, ie, nlon
    integer :: num_prog, num_diag, unit, stat
    logical :: used

    namelist/spectral_dynamics_nml/zero_topo, niter_topo, topo_file, topo_field, pdryini
    
    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=spectral_dynamics_nml)
    call close_file(unit) 

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
    allocate(nnp1(nwaves_oe,2))
   
    allocate(tmpg(ocny,ocnx)) 
    allocate(sin_latP(jsp:jep,isp:iep), cosm2_latP(jsp:jep,isp:iep))
    allocate(cosm_latP(jsp:jep,isp:iep))
    allocate(deg_latP(jsp:jep,isp:iep))
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

    call get_latsP(deglat = tmpg)
    deg_latP(jsp:jep,isp:iep) = tmpg(jsp:jep,isp:iep)
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
    
    forall(k=1:nlev) plev(k) = k
    forall(k=1:nlev+1) plevp(k) = k
    
    gaxis(1) = diag_axis_init('ocy',[(float(i),i=jsp,jep)],'degrees_N','Y', &
                long_name='y', domain_decomp=[1,ocny,jsp,jep])

    gaxis(2) = diag_axis_init('ocx',[(float(i),i=isp,iep)],'degrees_E','X', &
                long_name='x', domain_decomp=[1,ocnx,isp,iep])

    gaxis(3) = diag_axis_init('lev',plev,'','Z',long_name='')

    gaxis(4) = diag_axis_init('levp',plevp,'','Z',long_name='')

    call init_diag_out(gaxis,Time) 

    call get_spherical_wave(sph_wave,nnp1)
    
    do i = 1, size(sph_wave,1)
        if (sph_wave(i,1) == 0) then
            i0 = i
            exit
        endif
    enddo

    call init_bfiltr(sph_wave,trunc)

    ref_temp = 300.
    if (nlev>100) ref_temp=1500.
    
    call init_implicit(ak,bk,ref_temp,deltim,trunc)
    
    call init_horiz_diffusion(trunc,deltim,sl,sph_wave,bk)
   
    call init_data()

    call restore_spec_restart()

    gtopo = 0.; stopo = 0.; soro = 0.
#ifndef AQUAPLANET
    if (.not.zero_topo) then
        call read_topo()
    endif
#endif

    used = send_data(id_topo, gtopo(1,:,:))
    
    tmpmx = 0.
    if (nwaves_oe>0) tmpmx(1) = maxval(abs(satm(2)%prs))
    call spec_comm_max(tmpmx,1,commID)

    if (tmpmx(1)<=0.) then
        call mpp_error(WARNING,'Initial conditions not given or not proper, assuming a cold start')
        call mpp_error(WARNING,'-------------------COLD START-------------------------')
        call mpp_error(NOTE,'All tracers set to Zero')
        satm(1)%tr = 0.; satm(2)%tr = 0.
        call set_prs_temp_prof()
    end if
        
end subroutine init_spectral_dynamics

!--------------------------------------------------------------------------------   
subroutine set_prs_temp_prof()
!--------------------------------------------------------------------------------   
    character(len=100) :: cval
    real, allocatable :: tmp2d(:,:), tmp3d(:,:,:), axin(:)
    real, allocatable :: prsl(:,:,:)
    integer :: i, siz(4), j
    logical :: used

#ifndef AQUAPLANET
    if (.not.file_exist(temp_file)) call mpp_error(FATAL,trim(temp_file)//' does not exist')
    if (.not.field_exist(temp_file,'temp')) call mpp_error(FATAL,'field temp does not exist in file '&
                                          //trim(temp_file))
    if (.not.field_exist(temp_file,'pres')) call mpp_error(FATAL,'field pres does not exist in file '&
                                          //trim(temp_file))

    allocate(tmp2d(oc_ny(),oc_nx()))
    call read_data(temp_file,'pres',tmp2d)
    gatm(1)%prs(1,:,:) = tmp2d(jsp:jep,isp:iep) * 0.001 !->pascals to centibar
    deallocate(tmp2d)
    allocate(prsl(nlev,jsp:jep,isp:iep))
    call get_pressure_at_levels(gatm(1)%prs(1,:,:), prsl=prsl)
    call mpp_error(NOTE,'Surface pressure is set from file '//trim(temp_file))


    call field_size(temp_file,'level',siz)
    allocate(axin(siz(1)), tmp3d(siz(1),oc_ny(),oc_nx()))
    allocate(prsl(nlev,jsp:jep,isp:iep))
    call read_data(temp_file,'level',axin)
    call read_data(temp_file,'temp',tmp3d)
    axin = axin * 0.1 !-> mb to cb
    do i = isp, iep
        do j = jsp, jep
            call interp_vert(tmp3d(:,j,i), gatm(1)%tem(:,j,i), axin, prsl(:,j,i))
        end do
    end do
    call mpp_error(NOTE,'Temperature profiles are set from file '//trim(temp_file))
    deallocate(axin, tmp3d, prsl)
#else
    call mpp_error(NOTE,'-----AQUAPLANET RUN--------')

    call mpp_error(NOTE,'Setting up a constant initial Surface Pressure')
    gatm(1)%prs(1,:,:) = pdryini 
    allocate(prsl(nlev,jsp:jep,isp:iep))
    call get_pressure_at_levels(gatm(1)%prs(1,:,:), prsl=prsl)

    call mpp_error(NOTE,'Setting up AQUAPLANET Temperature profiles')
    do i = isp, iep
        do j = jsp, jep
            call aquape_init_temp(deg_latP(j,i), prsl(:,j,i)*10., gatm(1)%tem(:,j,i))
        end do
    end do
    deallocate(prsl)
#endif

    do i = 1, size(satm)
        call grid_to_spherical(gatm(1)%tem, satm(i)%tem, do_trunc=.true.)
    end do
    call spherical_to_grid(satm(2)%tem,grid=gatm(1)%tem)
    used =  send_data(id_tem_ini,gatm(1)%tem)

    gatm(1)%prs(1,:,:) = log(gatm(1)%prs(1,:,:))
    do i = 1, size(satm)
        call grid_to_spherical(gatm(1)%prs, satm(i)%prs, do_trunc=.true.)
    end do
    call spherical_to_grid(satm(2)%prs,grid=gatm(1)%prs)
    gatm(1)%prs = exp(gatm(1)%prs)
    used =  send_data(id_prs_ini,gatm(1)%prs(1,:,:))

    return
end subroutine set_prs_temp_prof


!--------------------------------------------------------------------------------   
subroutine interp_vert (fldin,fldout,axin,axout)
!--------------------------------------------------------------------------------   
    real, dimension(:), intent(in) :: fldin, axin, axout
    real, dimension(:), intent(out) :: fldout

    integer :: i1, i2, ni, j
    real :: w1, w2, tmp(size(axin)), w12
   
    ni = size(axin) 

    do j = 1, size(axout)
        tmp = axin-axout(j)
        i1 = 0; i2 = 0
        if (any(tmp>=0)) i1 = minloc(tmp, 1, mask=tmp>=0.)
        if (any(tmp<=0)) i2 = maxloc(tmp, 1, mask=tmp<=0.)

        if(i1<=0.and.i2<=0) call mpp_error(FATAL,'interp_vert: both i1 and i2 <= 0') 

        if (i1<=0) i1=i2
        if (i2<=0) i2=i1

        w1 = abs(axout(j)-axin(i2))
        w2 = abs(axout(j)-axin(i1))

        if (w1==0.or.w2==0.) then
            w1 = 1.; w2 = 1.
        endif

        w12 = w1 + w2

        w1 = w1/w12; w2=w2/w12

        fldout(j) = fldin(i1)*w1+fldin(i2)*w2
    end do

    return

end subroutine interp_vert

subroutine read_topo()
    real, allocatable :: topog(:,:)
    real :: topomin(1)
    character(len=100) :: cval
    integer :: n

    if (.not.file_exist(topo_file)) call mpp_error(FATAL,'init_spectral_dynamics: zero_topo is False but '// &
                                            trim(topo_file)//' not present')

    if (.not.field_exist(topo_file,trim(topo_field))) call mpp_error(FATAL,'init_spectral_dynamics: '//&
                'cannot find field '//trim(topo_field)//' in file '//trim(topo_file))

    call mpp_error(NOTE,'init_spectral_dynamics: reading new topography from '//trim(topo_file))

    allocate(topog(oc_ny(),oc_nx())) 

    call read_data(topo_file,topo_field,topog)

    gtopo(1,:,:) = topog(jsp:jep,isp:iep)

    call grid_to_spherical(gtopo,soro,do_trunc=.true.)

    do n = 1, niter_topo
        call spherical_to_grid(soro,grid=gtopo)
        topomin = minval(gtopo)
        call spec_comm_min(topomin, 1, commID)
        if (topomin(1) > 0.) exit
        where(gtopo<0.) gtopo = 0.
        call grid_to_spherical(gtopo,soro,do_trunc=.true.)
    end do

    call spherical_to_grid(soro,grid=gtopo)
    topomin = minval(gtopo)
    call spec_comm_min(topomin, 1, commID)
    write(cval,*) topomin(1)
    call mpp_error(NOTE,'Minimum value of Topography is '//trim(adjustl(cval)))

    stopo(1,:,:) = soro(1,:,:)*(GA2*real(nnp1))
    deallocate(topog)

    return
end subroutine read_topo


subroutine end_spectral_dynamics()
    integer :: i

    call save_spec_restart() 
 
    deallocate(sucos)
    deallocate(svcos)
    deallocate(stmp3d)
    deallocate(stmp3d1)
    
    do i = 1, 3
        deallocate(satm(i)%vor)
        deallocate(satm(i)%div)
        deallocate(satm(i)%tem)
        deallocate(satm(i)%tr)
        deallocate(satm(i)%prs)
    enddo
    deallocate(stopo)
    deallocate(soro)
    do i = 1, 2
        deallocate(gatm(i)%u)
        deallocate(gatm(i)%v)
        deallocate(gatm(i)%tem)
        deallocate(gatm(i)%tr)
        deallocate(gatm(i)%prs)
    enddo
    
    deallocate(dphi%u)
    deallocate(dphi%v)
    deallocate(dphi%tem)
    deallocate(dphi%tr)
    deallocate(dphi%prs)
    
    deallocate(dlam%u)
    deallocate(dlam%v)
    deallocate(dlam%tem)
    deallocate(dlam%tr)
    deallocate(dlam%prs)
    
    deallocate(dt%u)
    deallocate(dt%v)
    deallocate(dt%tem)
    deallocate(dt%tr)
    deallocate(dt%prs)
    
    deallocate(div)
    deallocate(vor)
    deallocate(spdmax)

    call end_transforms()
    
end subroutine end_spectral_dynamics 

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
    allocate(soro(1,nwaves_oe,2))
    stopo = 0.
    soro = 0.
    
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
    allocate(gtopo(1,jsp:jep,isp:iep))
    allocate(spdmax(nlev))

    spdmax = 0.

    do i = 1, 2
        nm='_m'
        if (i==2) nm='_n'
        satm(i)%vor = 0.
        satm(i)%div = 0.
        satm(i)%tem = 0.
        satm(i)%prs = 0.
        satm(i)%tr = 0.
        idx=register_spec_restart('vor'//nm,satm(i)%vor,.false.,0.)
        idx=register_spec_restart('div'//nm,satm(i)%div,.false.,0.)
        idx=register_spec_restart('tem'//nm,satm(i)%tem,.false.,0.)
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

    id_tem_ini = register_static_field(rou, 'tem_ini', axis1, 'Temperature', 'K')
    id_prs_ini = register_static_field(rou, 'prs_ini', axis1(2:3), 'Surface Pressure', 'cbar')
    id_topo = register_static_field(rou, 'topo', axis1(2:3), 'Topography', 'm')

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
    real :: val
    integer :: i, j, k, ntr
    logical :: used

    twotimes = twotimes + 1
    if (twotimes==2) then
        call save_wisdom()
    endif

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
  
    call spec_comm_max(spdmax, size(spdmax), commID)
    
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
    logical :: used

    call calc_mass_corr(exp(gatm(2)%prs(1,:,:)), tr, moist_ind, pcorr)

    do k = 1, size(u,1)
        gatm(2)%u(k,jsp:jep,isp:iep) = u(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
        gatm(2)%v(k,jsp:jep,isp:iep) = v(k,jsp:jep,isp:iep) * cosm_latP(jsp:jep,isp:iep)
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
        if (i0>0) satm(3)%prs(1,i0,1) = cmplx(pcorr,0.)
        
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
        if (spdmax(k)==0.) cycle
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
