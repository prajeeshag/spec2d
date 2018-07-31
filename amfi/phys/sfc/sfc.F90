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

use data_override_mod, only : data_override

use diag_manager_mod, only : register_static_field, register_diag_field
use diag_manager_mod, only : diag_axis_init, send_data

use albedo_mod, only : init_albedo, setalb_lnd, setalb_sice, setalb_ocean

implicit none
private

public :: init_sfc, get_albedo, get_tskin, get_land_frac, get_emis

real, allocatable, dimension(:,:)   :: fland, focn, fice, cellarea
real, allocatable, dimension(:,:)   :: tslnd, sst, tsice
real, allocatable, dimension(:,:,:) :: smc, stc, slc
real, allocatable, dimension(:,:)   :: vegfrac, tg3, sheleg, snwdph
real, allocatable, dimension(:,:)   :: canopy, trans, sncovr, zorl
real, allocatable, dimension(:,:)   :: alvsf, alvwf, alnsf, alnwf, facsf, facwf
real, allocatable, dimension(:,:)   :: hprif, emis_ref
real, allocatable, dimension(:,:)   :: hsnow_sice, hice

integer, allocatable, dimension(:,:) :: soiltype, vegtype, slopetype

logical, allocatable, dimension(:,:)   :: lland, locn, lice, lsea

type(domain2D), pointer :: domain => NULL()

real :: deltim=0.

integer :: is, ie, ilen
integer :: js, je, jlen
integer, parameter :: nllnd=4, nlice=2

type(restart_file_type) :: sfcres 

integer :: id_fland, id_focn, id_fice, id_cellarea
integer :: id_tslnd, id_sst, id_tsice
integer :: id_smc, id_stc, id_slc
integer :: id_vegfrac, id_tg3, id_sheleg, id_snwdph
integer :: id_canopy, id_trans, id_sncovr, id_zorl
integer :: id_alvsf, id_alvwf, id_alnsf, id_alnwf, id_facsf, id_facwf
integer :: id_hprif
integer :: id_hsnow_sice, id_hice

logical :: initialized=.false.

contains

!--------------------------------------------------------------------------------   
subroutine init_sfc(Time,deltim_in,domain_in,axes_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2D), target :: domain_in
    real, intent(in) :: deltim_in
    integer, intent(in) :: axes_in(:)

    deltim = deltim_in 
    
    domain => domain_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    ilen = ie-is+1
    jlen = je-js+1

    call sfc_diag_init(axes_in,Time)
  
    allocate( cellarea(js:je,is:ie) )
    allocate( fland(js:je,is:ie) )
    allocate( focn(js:je,is:ie) )
    allocate( fice(js:je,is:ie) )
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
    allocate( hsnow_sice(js:je,is:ie) )
    allocate( hice(js:je,is:ie) )
    allocate( emis_ref(js:je,is:ie) )
   
    call init_land(Time)
    call init_albedo() 

    initialized = .true.

end subroutine init_sfc

!--------------------------------------------------------------------------------   
subroutine sfc_diag_init(axes,Time)
    integer, intent(in) :: axes(2)
    type(time_type), intent(in) :: Time

    integer :: id_lice, id_lland

    id_facsf = register_static_field('amfi_sfc','facsf',axes,'facsf','1')
    id_facwf = register_static_field('amfi_sfc','facwf',axes,'facwf','1')

    return
end subroutine sfc_diag_init


!--------------------------------------------------------------------------------   
subroutine init_land(Time)
!-------------------------------------------------------------------------------- 
    type(time_type), intent(in) :: Time
    character(len=32) :: gridfile='INPUT/grid_spec.nc'
    real, allocatable :: tmp(:,:), tmpt(:,:) 
    integer :: isg, ieg, jsg, jeg, iglen
    integer :: indx, k
    logical :: ov, used

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

    lland = .false.
    locn = .false.
    lice = .false.

    focn = 1.0 - fland
    fice = 0.
    
    lland = fland>0.
    locn = focn>0.
    lice = fice>0.

    if (any(lland)) then
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

        call data_override('ATM','facsf',facsf,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'init_land: data_override failed for facsf !')
        used = send_data(id_facsf,facsf)

        call data_override('ATM','facwf',facwf,Time,ov)
        if(.not.ov) facwf = 1. - facsf
        used = send_data(id_facwf,facwf)
    endif

    if (file_exist('INPUT/emis_ref.nc')) then
        call read_data('INPUT/emis_ref.nc','emis',tmpt)
        emis_ref(js:je,is:ie) = tmpt(js:je,is:ie)
        call mpp_error(NOTE,'Using background land emissivity from emis_ref.nc')
    else
        emis_ref = 0.95
        call mpp_error(NOTE,'Using constant land emissivity = 0.95')
    endif

    indx = register_restart_field(sfcres, 'sfc_res', 'sheleg', sheleg, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'snwdph', snwdph, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'tslnd', tslnd, &
                    domain, mandatory=.true., data_default=283.15)

    indx = register_restart_field(sfcres, 'sfc_res', 'stc', stc, &
                    domain, mandatory=.true., data_default=283.15)

    indx = register_restart_field(sfcres, 'sfc_res', 'smc', smc, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'slc', slc, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'canopy', canopy, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'trans', trans, &
                    domain, mandatory=.false., data_default=0.)

    indx = register_restart_field(sfcres, 'sfc_res', 'sncovr', sncovr, &
                    domain, mandatory=.false., data_default=0.)

    call restore_state(sfcres)

    call save_restart(sfcres,'init')
    return

end subroutine init_land


!--------------------------------------------------------------------------------   
subroutine get_land_frac(frac)
!--------------------------------------------------------------------------------   
    real, intent(out) :: frac(:,:)

    if (.not.initialized) call mpp_error(FATAL,'sfc_mod: not initialized')
   
    frac = fland

end subroutine get_land_frac


!--------------------------------------------------------------------------------   
subroutine get_tskin(Time,tskin)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(out) :: tskin(:,:)
    logical :: ov

    if (.not.initialized) call mpp_error(FATAL,'sfc_mod: not initialized')

    call data_override('ATM','sst',sst,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'get_tskin: data_override failed for sst !')
    call data_override('ATM','tsice',tsice,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'get_tskin: data_override failed for tsice !')
    call data_override('ATM','fice',fice,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'get_tskin: data_override failed for fice !')

    tskin = tslnd * fland &
          + sst * focn * (1.-fice) &
          + tsice * focn * fice
            
end subroutine get_tskin


!--------------------------------------------------------------------------------   
subroutine get_albedo(Time,coszen,sfcalb)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in)  :: coszen(:,:)
    real, intent(out) :: sfcalb(:,:,:)
    logical :: ov
    real, dimension(size(sfcalb,2),size(sfcalb,3),size(sfcalb,1)) :: alblnd, albsice, albocn
 
    integer :: imax, i

    if (.not.initialized) call mpp_error(FATAL,'sfc_mod: not initialized')

    imax = size(lland,1) * size(lland,2)

    alblnd = 0.
    albsice = 0.
    albocn = 0.

    if (any(lland)) then 
        call data_override('ATM','alvsf',alvsf,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'get_albedo: data_override failed for alvsf !')
        call data_override('ATM','alnsf',alnsf,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'get_albedo: data_override failed for alnsf !')
        call data_override('ATM','alvwf',alvwf,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'get_albedo: data_override failed for alvwf !')
        call data_override('ATM','alnwf',alnwf,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'get_albedo: data_override failed for alnwf !')

        call setalb_lnd(imax,lland,sheleg,zorl,coszen,hprif, &
                        alvsf,alnsf,alvwf,alnwf,facsf,facwf,  &
                        alblnd)
    endif

    call data_override('ATM','tsice',tsice,Time)
    call data_override('ATM','hsnow_sice',hsnow_sice,Time)
    call data_override('ATM','hice',hice,Time)
    lice = (fice>0.)
    call setalb_sice(imax,lice,hsnow_sice,hice,tsice,albsice)

    locn = (focn>0.).and.((1.-fice)>0.)
    call setalb_ocean(imax,locn,coszen,albocn)

    do i = 1, size(sfcalb,1)
    sfcalb(i,:,:) = fland * alblnd(:,:,i) &
                  + focn * fice * albsice(:,:,i) &
                  + focn * (1.-fice) * albocn(:,:,i)
    enddo

    return

end subroutine get_albedo

!--------------------------------------------------------------------------------   
subroutine get_emis(sfcemis)
!--------------------------------------------------------------------------------   
    implicit none
    real, dimension(:,:), intent(out) :: sfcemis
    real :: fsno1(size(fland,1),size(fland,2))
    real, parameter ::  emsref(8) = [0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99]

    fsno1 = 1.0 - sncovr
    sfcemis = emis_ref*fsno1 + emsref(8)*sncovr ! land

    fsno1 = 1.0 - fland

    sfcemis = sfcemis * fland  &
            + fsno1 * (1.-fice) * emsref(1) & ! sea point
            + fsno1 * fice * emsref(7) ! sea-ice

    return
end subroutine get_emis

end module sfc_mod

