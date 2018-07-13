module sfc_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain
use mpp_domains_mod, only : mpp_get_global_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use time_manager_mod, only : time_type

use fms_mod, only : open_namelist_file, close_file, file_exist

use fms_io_mod, only : read_data, restart_file_type, register_restart_field
use fms_io_mod, only : restore_state, save_restart

use albedo_mod, only : init_albedo, setalb_lnd, setalb_sice, setalb_ocean

implicit none
private

public :: init_sfc, get_albedo

real, allocatable, dimension(:,:)   :: fland, focn, fice, fsea, cellarea
real, allocatable, dimension(:,:)   :: tslnd, sst, tsice
real, allocatable, dimension(:,:,:) :: smc, stc, slc
real, allocatable, dimension(:,:)   :: vegfrac, tg3, sheleg, snwdph
real, allocatable, dimension(:,:)   :: canopy, trans, sncovr, zorl
real, allocatable, dimension(:,:)   :: alvsf, alvwf, alnsf, alnwf, facsf, facwf
real, allocatable, dimension(:,:)   :: hprif

integer, allocatable, dimension(:,:) :: soiltype, vegtype, slopetype

logical, allocatable, dimension(:,:)   :: lland, locn, lice, lsea

type(domain2D), pointer :: domain => NULL()

real :: deltim=0.

integer :: is, ie, ilen
integer :: js, je, jlen
integer, parameter :: nllnd=4, nlice=2

type(restart_file_type) :: sfcres 

contains

!--------------------------------------------------------------------------------   
subroutine init_sfc(Time,deltim_in,domain_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2D), target :: domain_in
    real, intent(in) :: deltim_in

    deltim = deltim_in 
    
    domain => domain_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    ilen = ie-is+1
    jlen = je-js+1
  
    allocate( cellarea(js:je,is:ie) )
    allocate( fland(js:je,is:ie) )
    allocate( focn(js:je,is:ie) )
    allocate( fice(js:je,is:ie) )
    allocate( fsea(js:je,is:ie) )
    allocate( tslnd(js:je,is:ie) )
    allocate( sst(js:je,is:ie) )
    allocate( tsice(js:je,is:ie) )
    allocate( smc(nllnd,js:je,is:ie) )
    allocate( stc(nllnd,js:je,is:ie) )
    allocate( slc(nllnd,js:je,is:ie) )
    allocate( vegfrac(js:je,is:ie) )
    allocate( tg3(js:je,is:ie) )
    allocate( sheleg(js:je,is:ie) )
    allocate( snwdph(js:je,is:ie) )
    allocate( canopy(js:je,is:ie) )
    allocate( trans(js:je,is:ie) )
    allocate( sncovr(js:je,is:ie) )
    allocate( zorl(js:je,is:ie) )
    allocate( alvsf(js:je,is:ie) )
    allocate( alvwf(js:je,is:ie) )
    allocate( alnsf(js:je,is:ie) )
    allocate( alnwf(js:je,is:ie) )
    allocate( facsf(js:je,is:ie) )
    allocate( facwf(js:je,is:ie) )
    allocate( soiltype(js:je,is:ie) )
    allocate( vegtype(js:je,is:ie) )
    allocate( slopetype(js:je,is:ie) )
    allocate( lland(js:je,is:ie) )
    allocate( locn(js:je,is:ie) )
    allocate( lice(js:je,is:ie) )
    allocate( lsea(js:je,is:ie) )
    allocate( hprif(js:je,is:ie) )
   
    call init_land()
    call init_albedo() 

end subroutine init_sfc

!--------------------------------------------------------------------------------   
subroutine init_land()
!-------------------------------------------------------------------------------- 
    character(len=32) :: gridfile='INPUT/grid_spec.nc'
    real, allocatable :: tmp(:,:), tmpt(:,:) 
    integer :: isg, ieg, jsg, jeg, iglen
    integer :: indx

    if (.not.file_exist(trim(gridfile))) &
        call mpp_error(FATAL,trim(gridfile)//' does not exist!!')

    call mpp_get_global_domain(domain,jsg,jeg,isg,ieg)

    allocate(tmp(isg:ieg,jsg:jeg))
    allocate(tmpt(jsg:jeg,isg:ieg))

    call read_data(gridfile,'AREA_LND',tmp)
    tmpt = transpose(tmp)
    fland(js:je,is:ie) = tmpt(js:je,is:ie)

    call read_data(gridfile,'AREA_LND_CELL',tmp)
    tmpt = transpose(tmp)
    cellarea(js:je,is:ie) = tmpt(js:je,is:ie)
    fland = fland/cellarea

    fsea = 1.0 - fland
    focn = fsea
    fice = 0.
    
    lland = fland>0.
    lsea = fsea>0.
    locn = focn>0.
    lice = fice>0.

    call read_data('INPUT/soiltype','soiltype',tmpt)
    soiltype(js:je,is:ie) = int(tmpt(js:je,is:ie)) 

    call read_data('INPUT/vegtype','vegtype',tmpt)
    vegtype(js:je,is:ie) = int(tmpt(js:je,is:ie)) 

    call read_data('INPUT/slopetype','slopetype',tmpt)
    slopetype(js:je,is:ie) = int(tmpt(js:je,is:ie)) 

    call read_data('INPUT/tg3','tg3',tmpt)
    tg3(js:je,is:ie) = tmpt(js:je,is:ie)

    call read_data('INPUT/mtn.nc','mtn0',tmpt)
    hprif(js:je,is:ie) = tmpt(js:je,is:ie)

    !indx = register_restart_field(sfcres, 'sfc_res', 'sheleg', sheleg, &
    !                domain, mandatory=.false., data_default=0.)

    !indx = register_restart_field(sfcres, 'sfc_res', 'snwdph', snwdph, &
    !                domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'tslnd', tslnd, &
                    domain, mandatory=.true., data_default=283.15)

    !indx = register_restart_field(sfcres, 'sfc_res', 'stc', stc, &
    !                domain, mandatory=.true., data_default=283.15)

    !indx = register_restart_field(sfcres, 'sfc_res', 'smc', smc, &
    !                domain, mandatory=.false., data_default=0.)

    !indx = register_restart_field(sfcres, 'sfc_res', 'slc', slc, &
    !                domain, mandatory=.false., data_default=0.)

    !indx = register_restart_field(sfcres, 'sfc_res', 'canopy', canopy, &
    !                domain, mandatory=.false., data_default=0.)

    !indx = register_restart_field(sfcres, 'sfc_res', 'trans', trans, &
    !                domain, mandatory=.false., data_default=0.)

    !indx = register_restart_field(sfcres, 'sfc_res', 'sncovr', sncovr, &
    !                domain, mandatory=.false., data_default=0.)

    call restore_state(sfcres)

    call write_data('check_read','stc',stc,domain=domain)
    call write_data('check_read','tslnd',tslnd,domain=domain)

    call save_restart(sfcres,'init')
    return
end subroutine init_land

subroutine get_albedo(Time,sfcalb)
    type(time_type), intent(in) :: Time
    real, intent(in) :: sfcalb(:,:,:)

end subroutine get_albedo


end module sfc_mod
