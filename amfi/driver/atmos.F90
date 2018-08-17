
module atmos_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_io_mod, only : mpp_open, MPP_RDONLY

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use fms_io_mod, only : restart_file_type, reg_rf => register_restart_field
use fms_io_mod, only : restore_state, save_restart, fms_io_exit

use data_override_mod, only : data_override_init

use time_manager_mod, only : time_type

use tracer_manager_mod, only : tracer_manager_init, get_number_tracers, get_tracer_index

use field_manager_mod, only : MODEL_ATMOS

use spectral_dynamics_mod, only : init_spectral_dynamics, spectral_dynamics
use spectral_dynamics_mod, only : get_latsP, get_lonsP, finish_spectral_dynamics, &
                                  save_spec_restart, restore_spec_restart
  
use phys_mod, only : init_phys, phys

use ocpack_mod, only : oc_nx, oc_ny, npack=>oc_npack, init_ocpack 

implicit none
private

public :: init_atmos, update_atmos, end_atmos

integer :: nlev=64
integer :: isp, iep, ilenp, ocnx
integer :: jsp, jep, jlenp, ocny
integer :: ishuff, ntrac=2
integer :: unit, trunc=62
integer :: num_lat=94, maxlon=0
logical :: reduced=.true.
logical :: packed=.true.

real :: deltim=600.

type(domain2d) :: domain_g

real, allocatable :: u(:,:,:)
real, allocatable :: u1(:,:,:)
real, allocatable :: u2(:,:,:)
real, allocatable :: v(:,:,:)
real, allocatable :: v1(:,:,:)
real, allocatable :: v2(:,:,:)
real, allocatable :: vvel1(:,:,:)
real, allocatable :: p(:,:)
real, allocatable :: p1(:,:)
real, allocatable :: tem(:,:,:)
real, allocatable :: tem1(:,:,:)
real, allocatable :: tem2(:,:,:)
real, allocatable :: tr(:,:,:,:)
real, allocatable :: tr1(:,:,:,:)
real, allocatable :: tr2(:,:,:,:)

real, allocatable :: lat_deg(:,:), lon_deg(:,:)

type(restart_file_type) :: rstrt

character(len=8) :: fldnm

character(len=8) :: moist_tracer_names(10)
integer :: moist_tracer_ind(10) = 0
integer :: nmoist_tracers = 0

namelist/atmos_nml/ trunc, maxlon, num_lat, reduced, packed, nlev

contains


!--------------------------------------------------------------------------------   
subroutine init_atmos(Time,deltim_in)
!--------------------------------------------------------------------------------
    type(time_type), intent(in) :: Time
    real, intent(in) :: deltim_in
    integer :: layout(2), tmp, idx

    integer :: num_prog, num_diag, n
    integer :: axis(4)

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=atmos_nml)
    call close_file(unit)

    if (maxlon>0) then
        call init_ocpack(num_lat, trunc, 1, max_lon=maxlon, isreduced=reduced, ispacked=packed)
    else
        call init_ocpack(num_lat, trunc, 1, isreduced=reduced, ispacked=packed)
    end if

    ocnx = oc_nx()
    ocny = oc_ny()

    layout = [1,mpp_npes()]

    deltim = deltim_in

    call tracer_manager_init()
    call get_number_tracers(MODEL_ATMOS,ntrac,num_prog,num_diag)

    call mpp_define_domains( [1,ocny,1,ocnx], layout, domain_g, kxy=1)
    call mpp_get_compute_domain(domain_g, jsp, jep, isp, iep)
    ilenp = iep - isp + 1
    jlenp = jep - jsp + 1

    call data_override_init(Atm_domain_in=domain_g)
 
    allocate(u(nlev,jsp:jep,isp:iep))
    allocate(v(nlev,jsp:jep,isp:iep))
    allocate(tem(nlev,jsp:jep,isp:iep))
    allocate(tr(nlev,jsp:jep,isp:iep,ntrac))
    allocate(p(jsp:jep,isp:iep))
    
    allocate(u1(nlev,jsp:jep,isp:iep))
    allocate(v1(nlev,jsp:jep,isp:iep))
    allocate(vvel1(nlev,jsp:jep,isp:iep))
    allocate(tem1(nlev,jsp:jep,isp:iep))
    allocate(tr1(nlev,jsp:jep,isp:iep,ntrac))
    allocate(p1(jsp:jep,isp:iep))

    allocate(u2(nlev,jsp:jep,isp:iep))
    allocate(v2(nlev,jsp:jep,isp:iep))
    allocate(tem2(nlev,jsp:jep,isp:iep))
    allocate(tr2(nlev,jsp:jep,isp:iep,ntrac))

    allocate(lat_deg(ocny, ocnx), lon_deg(ocny,ocnx))

    idx = reg_rf(rstrt, 'amfi_res', 'tmp', tmp, domain_g, mandatory=.false.) !-> Just for registering the restart filename

    call init_spectral_dynamics(Time, nlev, trunc, domain_g, deltim, rstrt, axis)

    call get_lonsP(deglon=lon_deg)
    call get_latsP(deglat=lat_deg)
    
    !call init_phys(Time,deltim*2,deltim,domain_g,nlev,lat_deg,lon_deg,rstrt,axis)

    !call restore_state(rstrt)
    call restore_spec_restart()
    !call save_restart(rstrt,'Ini')
    call save_spec_restart('Ini')

end subroutine init_atmos


!--------------------------------------------------------------------------------   
subroutine update_atmos(Time)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer :: ntr

    call spectral_dynamics(Time,u,v,tem,tr,p,u1,v1,tem1,tr1,p1,vvel1)

    !call write_data('rgloopa','u',u,domain_g)
    !call write_data('rgloopa','v',v,domain_g)
    !call write_data('rgloopa','p',p,domain_g)
    !call write_data('rgloopa','tem',tem,domain_g)
    !
    !do ntr = 1, ntrac
    !    write(fldnm,'(A,I1)') 'tr',ntr
    !    call write_data('rgloopa',fldnm,tr(:,:,:,ntr),domain_g) 
    !enddo
    !
    !call write_data('rgloopa','u1',u1,domain_g)
    !call write_data('rgloopa','v1',v1,domain_g)
    !call write_data('rgloopa','p1',p1,domain_g)
    !call write_data('rgloopa','tem1',tem1,domain_g)
    !
    !do ntr = 1, ntrac
    !    write(fldnm,'(A,I1)') 'tr',ntr
    !    call write_data('rgloopa',trim(fldnm)//'_1',tr1(:,:,:,ntr),domain_g) 
    !enddo

    !call phys(Time,tem,tr,p,tem1,tr1,p1,u1,v1,vvel1,tem2,tr2,u2,v2)

    call finish_spectral_dynamics(Time,tem2,tr2,u2,v2)

end subroutine update_atmos


!--------------------------------------------------------------------------------   
subroutine end_atmos(Time)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time

    call save_restart(rstrt)
    call save_spec_restart()

    deallocate(u)
    deallocate(v)
    deallocate(tem)
    deallocate(tr)
    deallocate(p)
    deallocate(u1)
    deallocate(v1)
    deallocate(tem1)
    deallocate(tr1)
    deallocate(p1)

end subroutine end_atmos


end module atmos_mod
