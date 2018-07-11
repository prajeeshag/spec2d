
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
use fms_io_mod, only : fms_io_exit 

use time_manager_mod, only : time_type

use spectral_dynamics_mod, only : init_spectral_dynamics, spectral_dynamics

use radiation_mod, only : init_radiation

implicit none
private

public :: init_atmos, update_atmos, end_atmos

integer :: nlon=192, nlat=94, nlev=64
integer :: isc, iec, ilen
integer :: jsc, jec, jlen
integer :: ishuff, ntrac=3
integer :: unit, trunc=62
real :: deltim=600.

type(domain2d) :: domain_g

real, allocatable :: u(:,:,:)
real, allocatable :: u1(:,:,:)
real, allocatable :: v(:,:,:)
real, allocatable :: v1(:,:,:)
real, allocatable :: p(:,:)
real, allocatable :: p1(:,:)
real, allocatable :: tem(:,:,:)
real, allocatable :: tem1(:,:,:)
real, allocatable :: tr(:,:,:,:)
real, allocatable :: tr1(:,:,:,:)

integer :: ntr
character(len=8) :: fldnm


namelist/atmos_nml/ trunc, nlon, nlat, nlev, deltim


contains


!--------------------------------------------------------------------------------   
subroutine init_atmos(Time,deltim_in)
!--------------------------------------------------------------------------------
    type(time_type), intent(in) :: Time
    real, intent(in) :: deltim_in
    call mpp_init()
    call fms_init()
    
    unit = open_namelist_file()
    read(unit,nml=atmos_nml)
    call close_file(unit)

    deltim = deltim_in
    
    ntrac = 3
    ishuff = 2
    if(mpp_npes()==1) ishuff=0
    
    call mpp_define_domains( [1,nlat,1,nlon], [1,mpp_npes()], domain_g, kxy=1, ishuff=ishuff)
    call mpp_get_compute_domain(domain_g, jsc, jec, isc, iec)
    ilen = iec-isc+1
    jlen = jec-jsc+1
    
    allocate(u(nlev,jsc:jec,isc:iec))
    allocate(v(nlev,jsc:jec,isc:iec))
    allocate(tem(nlev,jsc:jec,isc:iec))
    allocate(tr(nlev,jsc:jec,isc:iec,ntrac))
    allocate(p(jsc:jec,isc:iec))
    
    allocate(u1(nlev,jsc:jec,isc:iec))
    allocate(v1(nlev,jsc:jec,isc:iec))
    allocate(tem1(nlev,jsc:jec,isc:iec))
    allocate(tr1(nlev,jsc:jec,isc:iec,ntrac))
    allocate(p1(jsc:jec,isc:iec))
    
    call init_spectral_dynamics(nlon,nlat,nlev,trunc,ntrac,domain_g,deltim)
    
    call init_radiation(deltim,domain_g,ntrac,nlev)

end subroutine init_atmos


!--------------------------------------------------------------------------------   
subroutine update_atmos(Time)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time

    call spectral_dynamics(u,v,tem,tr,p,u1,v1,tem1,tr1,p1)
    
    call write_data('rgloopa','u',u,domain_g)
    call write_data('rgloopa','v',v,domain_g)
    call write_data('rgloopa','p',p,domain_g)
    call write_data('rgloopa','tem',tem,domain_g)
    
    do ntr = 1, ntrac
        write(fldnm,'(A,I1)') 'tr',ntr
        call write_data('rgloopa',fldnm,tr(:,:,:,ntr),domain_g) 
    enddo
    
    call write_data('rgloopa','u1',u1,domain_g)
    call write_data('rgloopa','v1',v1,domain_g)
    call write_data('rgloopa','p1',p1,domain_g)
    call write_data('rgloopa','tem1',tem1,domain_g)
    
    do ntr = 1, ntrac
        write(fldnm,'(A,I1)') 'tr',ntr
        call write_data('rgloopa',trim(fldnm)//'_1',tr1(:,:,:,ntr),domain_g) 
    enddo

end subroutine update_atmos



!--------------------------------------------------------------------------------   
subroutine end_atmos(Time)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time

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
