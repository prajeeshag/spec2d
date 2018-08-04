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

use fms_io_mod, only : read_data, restart_file_type, reg_rf => register_restart_field
use fms_io_mod, only : restore_state, save_restart

use data_override_mod, only : data_override

use diag_manager_mod, only : register_static_field, reg_df => register_diag_field
use diag_manager_mod, only : diag_axis_init, send_data

use constants_mod, only : KELVIN, RDGAS, GRAV, RVGAS, HLV, CP_AIR

use sat_vapor_pres_mod, only : compute_qs

use albedo_mod, only : init_albedo, setalb_lnd, setalb_sice, setalb_ocean

use seaice_mod, only : sfc_sice_drv, kmi, tgice, himin

use sfc_ocean_mod, only : sfc_ocean

use sfc_diff_mod, only : sfc_diff

implicit none
private

public :: init_sfc, get_land_frac, set_surface, do_surface

real, allocatable, dimension(:,:)   :: fland, cellarea
real, allocatable, dimension(:,:)   :: tslnd
real, allocatable, dimension(:,:,:) :: smc, stc, slc
real, allocatable, dimension(:,:)   :: tg3, sheleg, snwdph
real, allocatable, dimension(:,:)   :: canopy, trans, sncovr, zorl, zorlocn
real, allocatable, dimension(:,:)   :: alvsf, alvwf, alnsf, alnwf, facsf, facwf
real, allocatable, dimension(:,:)   :: hprif, emis_ref

real, allocatable, dimension(:,:)   :: hsnow_sice, hice, tsice, fice
real, allocatable, dimension(:,:,:)   :: tice

real, allocatable, dimension(:,:)   :: sst, focn

integer, allocatable, dimension(:,:) :: soiltype, vegtype, slopetype

logical, allocatable, dimension(:,:)   :: lland, locn
logical, allocatable, dimension(:,:,:)   :: lland3d

real, parameter :: con_fvirt = RVGAS/RDGAS-1.
real, parameter :: ELOCP = HLV/CP_AIR
real, parameter :: HVAPI = 1./HLV
real, parameter :: CPINV = 1./CP_AIR
real, parameter :: zorl_min = 0.001
real, parameter :: fice_min = 0.05

type(domain2D), pointer :: domain => NULL()

real :: deltim=0.

integer :: is, ie, ilen
integer :: js, je, jlen
integer, parameter :: lsoil=4

real, parameter :: zsoil(lsoil) = [-0.1,-0.4,-1.0,-2.0]
real :: sldpth(lsoil)

integer :: id_cd, id_cdq, id_rb, id_stress, id_ffmm, id_ffhh, id_uustar, &
           id_wind, id_fm10, id_fh2, id_qss, id_cmm, id_chh, id_evap, id_hflx

integer :: id_vegtype, id_soiltype, id_slopetype, id_canopy, id_tslnd, id_sheleg, id_stc,  &
           id_smc, id_slc, id_vegfrac, id_snwdph, id_nroot, id_eta, id_sheat, id_ec, id_edir,  &
           id_et, id_trans, id_esnow, id_drip, id_dew, id_beta, id_etp, id_ssoil, id_flx1,  &
           id_flx2, id_flx3, id_runoff1, id_runoff2, id_runoff3, id_snomlt, id_sncovr, id_rc,  &
           id_pc, id_rsmin, id_xlai, id_rcs, id_rct, id_rcq, id_rcsoil, id_soilw, id_soilm,  &
           id_smcwlt, id_smcdry, id_smcref, id_smcmax

integer :: id_facsf, id_facwf

logical :: initialized=.false.

contains

!--------------------------------------------------------------------------------   
subroutine init_sfc(Time,deltim_in,domain_in,axes_in,rstrt)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2D), target :: domain_in
    real, intent(in) :: deltim_in
    integer, intent(in) :: axes_in(:)
    type(restart_file_type), intent(inout) :: rstrt
    
    integer :: indx

    deltim = deltim_in 
    
    domain => domain_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    ilen = ie-is+1
    jlen = je-js+1

    allocate( cellarea(js:je,is:ie) )
    allocate( fland(js:je,is:ie) )
    allocate( focn(js:je,is:ie) )
    allocate( fice(js:je,is:ie) )
    allocate( tslnd(js:je,is:ie) )
    allocate( sst(js:je,is:ie) )
    allocate( tsice(js:je,is:ie) )
    allocate( tice(kmi,js:je,is:ie) )
    allocate( smc(lsoil,js:je,is:ie) )
    allocate( stc(lsoil,js:je,is:ie) )
    allocate( slc(lsoil,js:je,is:ie) )
    allocate( tg3(js:je,is:ie) )
    allocate( sheleg(js:je,is:ie) )
    allocate( snwdph(js:je,is:ie) )
    allocate( canopy(js:je,is:ie) )
    allocate( trans(js:je,is:ie) )
    allocate( sncovr(js:je,is:ie) )
    allocate( zorl(js:je,is:ie) )
    allocate( zorlocn(js:je,is:ie) )
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
    allocate( lland3d(lsoil,js:je,is:ie) )
    allocate( hprif(js:je,is:ie) )
    allocate( hsnow_sice(js:je,is:ie) )
    allocate( hice(js:je,is:ie) )
    allocate( emis_ref(js:je,is:ie) )
   
    sheleg = 0.
    indx = reg_rf(rstrt, '', 'sheleg', sheleg, domain, mandatory=.false.)

    snwdph = 0.
    indx = reg_rf(rstrt, '', 'snwdph', snwdph, domain, mandatory=.false.)

    tslnd = KELVIN
    indx = reg_rf(rstrt, '', 'tslnd', tslnd, domain, mandatory=.true.)

    stc = KELVIN
    indx = reg_rf(rstrt, '', 'stc', stc, domain, mandatory=.true.)

    smc = 0.
    indx = reg_rf(rstrt, '', 'smc', smc, domain, mandatory=.false.)

    stc = 0.
    indx = reg_rf(rstrt, '', 'slc', slc, domain, mandatory=.false.)

    canopy = 0.
    indx = reg_rf(rstrt, '', 'canopy', canopy, domain, mandatory=.false.)

    trans = 0.
    indx = reg_rf(rstrt, '', 'trans', trans, domain, mandatory=.false.)

    sncovr = 0.
    indx = reg_rf(rstrt, '', 'sncovr', sncovr, domain, mandatory=.false.)

    zorlocn = zorl_min
    indx = reg_rf(rstrt, '', 'zorlocn', zorlocn, domain, mandatory=.false.)

    tsice = tgice
    indx = reg_rf(rstrt, '', 'tsice', tsice, domain, mandatory=.false.)

    tice = tgice
    indx = reg_rf(rstrt, '', 'tice', tice, domain, mandatory=.false.)

    hice = himin
    indx = reg_rf(rstrt, '', 'hice', hice, domain, mandatory=.false.)

    hsnow_sice = 0.
    indx = reg_rf(rstrt, '', 'hsnow_sice', hsnow_sice, domain, mandatory=.false.)

    call sfc_diag_init(axes_in,Time)
  
    call init_land(Time)

    call init_albedo() 

    initialized = .true.

end subroutine init_sfc


!--------------------------------------------------------------------------------   
subroutine sfc_diag_init(axes,Time)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: axes(2)
    type(time_type), intent(in) :: Time

    integer :: id_kmi, id_lsoil, k, axesl(3), axesi(3)
    real :: tmp(kmi)

    id_lsoil = diag_axis_init('lsoil', zsoil, 'm', 'Z', long_name='Soil Levels')

    forall(k=1:kmi) tmp(k) = k 
    id_kmi = diag_axis_init('kmi', tmp, '', 'Z', long_name='Sea-Ice Levels')

    axesl(:) = [id_lsoil,axes(1),axes(2)]
    axesi(:) = [id_kmi,axes(1),axes(2)]

    id_facsf  = register_static_field('amfi_sfc','facsf',axes,'facsf','1')
    id_facwf  = register_static_field('amfi_sfc','facwf',axes,'facwf','1')

    id_cd      =  reg_df('am_sfc',  'cd',      axes,  Time,  'noname',  'nounit')
    id_cdq     =  reg_df('am_sfc',  'cdq',     axes,  Time,  'noname',  'nounit')
    id_rb      =  reg_df('am_sfc',  'rb',      axes,  Time,  'noname',  'nounit')
    id_stress  =  reg_df('am_sfc',  'stress',  axes,  Time,  'noname',  'nounit')
    id_ffmm    =  reg_df('am_sfc',  'ffmm',    axes,  Time,  'noname',  'nounit')
    id_ffhh    =  reg_df('am_sfc',  'ffhh',    axes,  Time,  'noname',  'nounit')
    id_uustar  =  reg_df('am_sfc',  'uustar',  axes,  Time,  'noname',  'nounit')
    id_wind    =  reg_df('am_sfc',  'wind',    axes,  Time,  'noname',  'nounit')
    id_fm10    =  reg_df('am_sfc',  'fm10',    axes,  Time,  'noname',  'nounit')
    id_fh2     =  reg_df('am_sfc',  'fh2',     axes,  Time,  'noname',  'nounit')
    id_qss     =  reg_df('am_sfc',  'qss',     axes,  Time,  'noname',  'nounit')
    id_cmm     =  reg_df('am_sfc',  'cmm',     axes,  Time,  'noname',  'nounit')
    id_chh     =  reg_df('am_sfc',  'chh',     axes,  Time,  'noname',  'nounit')
    id_evap    =  reg_df('am_sfc',  'evap',    axes,  Time,  'noname',  'nounit')
    id_hflx    =  reg_df('am_sfc',  'hflx',    axes,  Time,  'noname',  'nounit')

    id_vegtype   = reg_df('am_sfc', 'vegtype',   axes, Time, & 
                      'vegtype',    '1', missing_value = -99999.0)
    id_soiltype  = reg_df('am_sfc', 'soiltype',  axes, Time, & 
                      'soiltype',   '1', missing_value = -99999.0)
    id_slopetype = reg_df('am_sfc', 'slopetype', axes, Time, & 
                      'slopetype',  '1', missing_value = -99999.0)
    id_canopy    = reg_df('am_sfc', 'canopy',    axes, Time, & 
                      'canopy',     'mm', missing_value = -99999.0)
    id_tslnd     = reg_df('am_sfc', 'tslnd',     axes, Time, & 
                      'Land Surface Temperature', 'K', missing_value = -99999.0)
    id_sheleg    = reg_df('am_sfc', 'sheleg',    axes, Time, & 
                      'sheleg',     'mm', missing_value = -99999.0)
    id_stc       = reg_df('am_sfc', 'stc',       axesl, Time, & 
                      'Soil Temperature', 'K', missing_value = -99999.0)
    id_smc       = reg_df('am_sfc', 'smc',       axesl, Time, & 
                      'Soil Moisture', 'nounit', missing_value = -99999.0)
    id_slc       = reg_df('am_sfc', 'slc',       axesl, Time, & 
                      'Soil Liquid Water', 'nounit', missing_value = -99999.0)
    id_vegfrac   = reg_df('am_sfc', 'vegfrac',   axes, Time, & 
                      'Green Vegetation Fraction', '1', missing_value = -99999.0)
    id_snwdph    = reg_df('am_sfc', 'snwdph',    axes, Time, & 
                      'Snow Height', 'mm', missing_value = -99999.0)
    id_nroot     = reg_df('am_sfc', 'nroot',     axes, Time, & 
                      'No: of root layers', '1', missing_value = -99999.0)
    id_eta       = reg_df('am_sfc', 'eta',       axes, Time, & 
                      'Latent Heat Flux (Land', 'nounit', missing_value = -99999.0)
    id_sheat     = reg_df('am_sfc', 'sheat',     axes, Time, & 
                      'Sensible Heat Flux (Land)', 'nounit', missing_value = -99999.0)
    id_ec        = reg_df('am_sfc', 'ec',        axes, Time, & 
                      'ec',         'nounit', missing_value = -99999.0)
    id_edir      = reg_df('am_sfc', 'edir',      axes, Time, & 
                      'edir',       'nounit', missing_value = -99999.0)
    id_et        = reg_df('am_sfc', 'et',        axes, Time, & 
                      'et',         'nounit', missing_value = -99999.0)
    id_trans     = reg_df('am_sfc', 'trans',     axes, Time, & 
                      'trans',      'nounit', missing_value = -99999.0)
    id_esnow     = reg_df('am_sfc', 'esnow',     axes, Time, & 
                      'esnow',      'nounit', missing_value = -99999.0)
    id_drip      = reg_df('am_sfc', 'drip',      axes, Time, & 
                      'drip',       'nounit', missing_value = -99999.0)
    id_dew       = reg_df('am_sfc', 'dew',       axes, Time, & 
                      'dew',        'nounit', missing_value = -99999.0)
    id_beta      = reg_df('am_sfc', 'beta',      axes, Time, & 
                      'beta',       'nounit', missing_value = -99999.0)
    id_etp       = reg_df('am_sfc', 'etp',       axes, Time, & 
                      'etp',        'nounit', missing_value = -99999.0)
    id_ssoil     = reg_df('am_sfc', 'ssoil',     axes, Time, & 
                      'ssoil',      'nounit', missing_value = -99999.0)
    id_flx1      = reg_df('am_sfc', 'flx1',      axes, Time, & 
                      'flx1',       'nounit', missing_value = -99999.0)
    id_flx2      = reg_df('am_sfc', 'flx2',      axes, Time, & 
                      'flx2',       'nounit', missing_value = -99999.0)
    id_flx3      = reg_df('am_sfc', 'flx3',      axes, Time, & 
                      'flx3',       'nounit', missing_value = -99999.0)
    id_runoff1   = reg_df('am_sfc', 'runoff1',   axes, Time, & 
                      'runoff1',    'nounit', missing_value = -99999.0)
    id_runoff2   = reg_df('am_sfc', 'runoff2',   axes, Time, & 
                      'runoff2',    'nounit', missing_value = -99999.0)
    id_runoff3   = reg_df('am_sfc', 'runoff3',   axes, Time, & 
                      'runoff3',    'nounit', missing_value = -99999.0)
    id_snomlt    = reg_df('am_sfc', 'snomlt',    axes, Time, & 
                      'snomlt',     'nounit', missing_value = -99999.0)
    id_sncovr    = reg_df('am_sfc', 'sncovr',    axes, Time, & 
                      'sncovr',     'nounit', missing_value = -99999.0)
    id_rc        = reg_df('am_sfc', 'rc',        axes, Time, & 
                      'rc',         'nounit', missing_value = -99999.0)
    id_pc        = reg_df('am_sfc', 'pc',        axes, Time, & 
                      'pc',         'nounit', missing_value = -99999.0)
    id_rsmin     = reg_df('am_sfc', 'rsmin',     axes, Time, & 
                      'rsmin',      'nounit', missing_value = -99999.0)
    id_xlai      = reg_df('am_sfc', 'xlai',      axes, Time, & 
                      'xlai',       'nounit', missing_value = -99999.0)
    id_rcs       = reg_df('am_sfc', 'rcs',       axes, Time, & 
                      'rcs',        'nounit', missing_value = -99999.0)
    id_rct       = reg_df('am_sfc', 'rct',       axes, Time, & 
                      'rct',        'nounit', missing_value = -99999.0)
    id_rcq       = reg_df('am_sfc', 'rcq',       axes, Time, & 
                      'rcq',        'nounit', missing_value = -99999.0)
    id_rcsoil    = reg_df('am_sfc', 'rcsoil',    axes, Time, & 
                      'rcsoil',     'nounit', missing_value = -99999.0)
    id_soilw     = reg_df('am_sfc', 'soilw',     axes, Time, & 
                      'soilw',      'nounit', missing_value = -99999.0)
    id_soilm     = reg_df('am_sfc', 'soilm',     axes, Time, & 
                      'soilm',      'nounit', missing_value = -99999.0)
    id_smcwlt    = reg_df('am_sfc', 'smcwlt',    axes, Time, & 
                      'smcwlt',     'nounit', missing_value = -99999.0)
    id_smcdry    = reg_df('am_sfc', 'smcdry',    axes, Time, & 
                      'smcdry',     'nounit', missing_value = -99999.0)
    id_smcref    = reg_df('am_sfc', 'smcref',    axes, Time, & 
                      'smcref',     'nounit', missing_value = -99999.0)
    id_smcmax    = reg_df('am_sfc', 'smcmax',    axes, Time, & 
                      'smcmax',     'nounit', missing_value = -99999.0)

    return
end subroutine sfc_diag_init


!--------------------------------------------------------------------------------   
subroutine init_land(Time)
!-------------------------------------------------------------------------------- 
    type(time_type), intent(in) :: Time
    character(len=32) :: gridfile='INPUT/grid_spec.nc'
    real, allocatable :: tmp(:,:), tmpt(:,:) 
    integer :: isg, ieg, jsg, jeg, iglen
    integer :: indx, k, me
    logical :: ov, used

    me = 1
    if (mpp_pe()==mpp_root_pe()) me = 0

    call set_soilveg(me)

    sldpth(1) = - zsoil(1)
    sldpth(2:lsoil) = zsoil(1:lsoil-1) - zsoil(2:lsoil)

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

    focn = 1.0 - fland
    fice = 0.
    
    lland = fland>0.
    forall(k=1:lsoil) lland3d(k,:,:) = lland
    locn = focn>0.

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


    return

end subroutine init_land


!--------------------------------------------------------------------------------   
subroutine do_surface(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, &
                      sfcdsw, sfcnsw, sfcdlw, fprcp, lprcp, rb, ffmm, ffhh, &
                      qss, hflx, evap, stress, wind)
!--------------------------------------------------------------------------------   
    implicit none
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(js:je,is:ie) :: pgr, ugrs, vgrs, tgrs, qgrs, prsl
    real, intent(in), dimension(js:je,is:ie) :: prslki, sfcdsw, sfcnsw, sfcdlw
    real, intent(in), dimension(js:je,is:ie) :: fprcp, lprcp
    real, intent(out), dimension(js:je,is:ie) :: rb, stress, ffmm, ffhh, wind, &
                                                 qss, evap, hflx

    real, dimension(js:je,is:ie) :: cd, cdq, uustar, fm10, fh2, cmm, chh, fwat

    real, dimension(js:je,is:ie) :: cda, cdqa, rba, stressa, ffmma, ffhha, uustara, &
                                    winda, fm10a, fh2a, qssa, cmma, chha, evapa, hflxa

    logical :: used

    cd = 0.; cdq = 0.; rb = 0.; stress = 0.; ffmm = 0.; ffhh = 0.; uustar = 0.; wind = 0.;
    fm10 = 0.; fh2 = 0.; qss = 0.; cmm = 0.; chh = 0.; evap = 0.; hflx = 0.

    call do_ocean(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, &
                  cda, cdqa, rba, stressa, ffmma, ffhha, uustara, &
                  winda, fm10a, fh2a, qssa, cmma, chha, evapa, hflxa)

    fwat = 1. - fice
    
    cd      =  cd      +  cda      *  fwat  *  focn
    cdq     =  cdq     +  cdqa     *  fwat  *  focn
    rb      =  rb      +  rba      *  fwat  *  focn
    stress  =  stress  +  stressa  *  fwat  *  focn
    ffmm    =  ffmm    +  ffmma    *  fwat  *  focn
    ffhh    =  ffhh    +  ffhha    *  fwat  *  focn
    uustar  =  uustar  +  uustara  *  fwat  *  focn
    wind    =  wind    +  winda    *  fwat  *  focn
    fm10    =  fm10    +  fm10a    *  fwat  *  focn
    fh2     =  fh2     +  fh2a     *  fwat  *  focn
    qss     =  qss     +  qssa     *  fwat  *  focn
    cmm     =  cmm     +  cmma     *  fwat  *  focn
    chh     =  chh     +  chha     *  fwat  *  focn
    evap    =  evap    +  evapa    *  fwat  *  focn
    hflx    =  hflx    +  hflxa    *  fwat  *  focn

    call do_seaice(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, sfcdlw, &
                   sfcnsw, sfcdsw, fprcp, cda, cdqa, rba, stressa, ffmma, &
                   ffhha, uustara, winda, fm10a, fh2a, qssa, cmma, chha, evapa, &
                   hflxa)

    cd      =  cd      +  cda      *  fice  *  focn
    cdq     =  cdq     +  cdqa     *  fice  *  focn
    rb      =  rb      +  rba      *  fice  *  focn
    stress  =  stress  +  stressa  *  fice  *  focn
    ffmm    =  ffmm    +  ffmma    *  fice  *  focn
    ffhh    =  ffhh    +  ffhha    *  fice  *  focn
    uustar  =  uustar  +  uustara  *  fice  *  focn
    wind    =  wind    +  winda    *  fice  *  focn
    fm10    =  fm10    +  fm10a    *  fice  *  focn
    fh2     =  fh2     +  fh2a     *  fice  *  focn
    qss     =  qss     +  qssa     *  fice  *  focn
    cmm     =  cmm     +  cmma     *  fice  *  focn
    chh     =  chh     +  chha     *  fice  *  focn
    evap    =  evap    +  evapa    *  fice  *  focn
    hflx    =  hflx    +  hflxa    *  fice  *  focn

    call do_land(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, sfcdlw, &
                 sfcnsw, sfcdsw, fprcp, lprcp, cda, cdqa, rba, stressa, ffmma, &
                 ffhha, uustara, winda, fm10a, fh2a, qssa, cmma, chha, evapa, hflxa)

    cd      =  cd      +  cda      *  fland
    cdq     =  cdq     +  cdqa     *  fland
    rb      =  rb      +  rba      *  fland
    stress  =  stress  +  stressa  *  fland
    ffmm    =  ffmm    +  ffmma    *  fland
    ffhh    =  ffhh    +  ffhha    *  fland
    uustar  =  uustar  +  uustara  *  fland
    wind    =  wind    +  winda    *  fland
    fm10    =  fm10    +  fm10a    *  fland
    fh2     =  fh2     +  fh2a     *  fland
    qss     =  qss     +  qssa     *  fland
    cmm     =  cmm     +  cmma     *  fland
    chh     =  chh     +  chha     *  fland
    evap    =  evap    +  evapa    *  fland
    hflx    =  hflx    +  hflxa    *  fland

    !diagnostics
    used = send_data(id_cd    , cd    , Time) 
    used = send_data(id_cdq   , cdq   , Time)  
    used = send_data(id_rb    , rb    , Time)  
    used = send_data(id_stress, stress, Time)  
    used = send_data(id_ffmm  , ffmm  , Time)  
    used = send_data(id_ffhh  , ffhh  , Time)  
    used = send_data(id_uustar, uustar, Time)  
    used = send_data(id_wind  , wind  , Time)  
    used = send_data(id_fm10  , fm10  , Time)  
    used = send_data(id_fh2   , fh2   , Time)  
    used = send_data(id_qss   , qss   , Time)  
    used = send_data(id_cmm   , cmm   , Time)  
    used = send_data(id_chh   , chh   , Time)  
    used = send_data(id_evap  , evap  , Time)  
    used = send_data(id_hflx  , hflx  , Time)  

end subroutine do_surface


!--------------------------------------------------------------------------------   
subroutine do_ocean(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, &
                    cd, cdq, rb, stress, ffmm, ffhh, uustar, wind, fm10, fh2, &
                    qss, cmm, chh, evap, hflx)
!--------------------------------------------------------------------------------   
    implicit none
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(js:je,is:ie) :: pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki 
    real, dimension(js:je,is:ie), intent(out) :: cd, cdq, rb, stress, ffmm, ffhh, &
                                                 uustar, wind, fm10, fh2, qss, cmm, &
                                                 chh, evap, hflx

    logical, dimension(js:je,is:ie) :: mask
    real, dimension(js:je,is:ie) :: ddvel
    real, dimension(js:je,is:ie) :: slimsk

    integer :: im

    ddvel = 0.
    im = size(pgr,1)*size(pgr,2)

    mask = .false.
    where (focn>0.and.fice<1.) mask = .true.
    slimsk = 0.
    
    cd = 0.; cdq = 0.; rb = 0.; stress = 0.; ffmm = 0.; ffhh = 0.; uustar = 0.
    wind = 0.; fm10 = 0.; fh2 = 0.; qss = 0.; cmm = 0.; chh = 0.; evap = 0.; hflx = 0.

    call sfc_diff(im, pgr, ugrs, vgrs, tgrs, qgrs, sst, zorlocn, cd, cdq, rb, &
                  prsl, prslki, slimsk, stress, ffmm, ffhh, uustar, wind, &
                  ddvel, fm10, fh2, sst, mask)

    call sfc_ocean(im, pgr, wind, tgrs, qgrs, sst, cd, cdq, prsl, &
                           prslki, mask, qss, cmm, chh, evap, hflx)

end subroutine do_ocean


!--------------------------------------------------------------------------------   
subroutine do_seaice(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, sfcdlw, &
                     sfcnsw, sfcdsw, fprecip, cd, cdq, rb, stress, ffmm, ffhh, &
                     uustar, wind, fm10, fh2, qss, cmm, chh, evap, hflx)
!--------------------------------------------------------------------------------   
    implicit none
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(js:je,is:ie) :: pgr, ugrs, vgrs, tgrs
    real, intent(in), dimension(js:je,is:ie) :: qgrs, prsl, prslki
    real, intent(in), dimension(js:je,is:ie) :: sfcdlw, sfcnsw, sfcdsw
    real, intent(in), dimension(js:je,is:ie) :: fprecip
    
    real, dimension(js:je,is:ie), intent(out) :: cd, cdq, rb, stress, ffmm, ffhh, &
                                                 uustar, wind, fm10, fh2, qss, cmm, &
                                                 chh, evap, hflx

    real, dimension(js:je,is:ie) :: ddvel
    real, dimension(js:je,is:ie) :: sice_snwdph, sice_snowmt, gflux, zlvl
    logical, dimension(js:je,is:ie) :: mask
    real, dimension(js:je,is:ie) :: slimsk

    integer :: im

    im = size(pgr,1)*size(pgr,2)

    mask = .false.
    where (focn>0.and.fice>0.) mask = .true.
    slimsk = 1.

    ddvel = 0.

    sice_snwdph = 0.; sice_snowmt = 0.; gflux = 0.; zlvl = 0.
    cd = 0.; cdq = 0.; rb = 0.; stress = 0.; ffmm = 0.; ffhh = 0.; uustar = 0.
    wind = 0.; fm10 = 0.; fh2 = 0.; qss = 0.; cmm = 0.; chh = 0.; evap = 0.; hflx = 0.
    
    call sfc_diff(im, pgr, ugrs, vgrs, tgrs, qgrs, tsice, zorl, cd, cdq, rb, &
                  prsl, prslki, slimsk, stress, ffmm, ffhh, uustar, wind, &
                  ddvel, fm10, fh2, tsice, mask)

    call sfc_sice_drv(im, pgr, wind, tgrs, qgrs, deltim, sfcdlw, sfcnsw, &
                      sfcdsw, cd, cdq, prsl, prslki, mask, hice, tsice, &
                      hsnow_sice, fprecip, tice, sice_snwdph, qss, &
                      sice_snowmt, gflux, cmm, chh, zlvl, evap, hflx)

end subroutine do_seaice


!--------------------------------------------------------------------------------   
subroutine do_land(Time, pgr, ugrs, vgrs, tgrs, qgrs, prsl, prslki, sfcdlw, &
                   sfcnsw, sfcdsw, fprecip, lprecip, cd, cdq, rb, stress, ffmm, ffhh, &
                   uustar, wind, fm10, fh2, qss, cmm, chh, evap, hflx)
!--------------------------------------------------------------------------------   
    implicit none
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(js:je,is:ie) :: pgr, ugrs, vgrs, tgrs
    real, intent(in), dimension(js:je,is:ie) :: qgrs, prsl, prslki
    real, intent(in), dimension(js:je,is:ie) :: sfcdlw, sfcnsw, sfcdsw
    real, intent(in), dimension(js:je,is:ie) :: fprecip, lprecip
    
    real, dimension(js:je,is:ie), intent(out) :: cd, cdq, rb, stress, ffmm, ffhh, &
                                                 uustar, wind, fm10, fh2, qss, cmm, &
                                                 chh, evap, hflx

    real, dimension(js:je,is:ie) :: ddvel, sfcems, vegfrac, slimsk, eta, sheat, &
                ec, edir, et, esnow, drip, dew, beta, etp, ssoil, flx1, flx2, flx3, &
                runoff1, runoff2, runoff3, snomlt, rc, pc, rsmin, xlai, rcs, rcq, &
                rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax, rct, tsurf

    integer, dimension(js:je,is:ie) :: nroot
    logical, dimension(js:je,is:ie) :: mask, guess
    integer, parameter :: niter = 2
    integer :: im, iter, i, j
    logical :: ov, used

    im = size(pgr,1)*size(pgr,2)

    slimsk = 1.

    ddvel = 0.

    call get_emis(sfcems)

    call data_override('ATM','vegfrac',vegfrac,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'do_land: data_override failed for vegfrac!')
    where (vegfrac<0.01) vegfrac = 0.01

    tsurf = tslnd

    mask = .false.
    where (fland>0) mask = .true.
    guess = .false.

    do iter = 1, niter
    
        call sfc_diff(im, pgr, ugrs, vgrs, tgrs, qgrs, tslnd, zorl, cd, cdq, rb, &
                      prsl, prslki, slimsk, stress, ffmm, ffhh, uustar, wind, &
                      ddvel, fm10, fh2, tsurf, mask)

        if (iter == 1.and.niter>1) then
            where(wind<2..and.mask) guess = .true.
        endif

        call sfc_land_drv(im, deltim, fprecip, sfcdsw, sfcnsw, sfcdlw, sfcems, pgr, prsl, tgrs, &
                          prslki, wind, lprecip, qgrs, tg3, cd, cdq, vegtype, soiltype,  &
                          slopetype, canopy, tslnd, sheleg, stc, smc, slc, vegfrac, snwdph, &
                          nroot, eta, sheat, ec, edir, et, trans, esnow, drip, dew, &
                          beta, etp, ssoil, flx1, flx2, flx3, runoff1, runoff2, &
                          runoff3, snomlt, sncovr, rc, pc, rsmin, xlai, rcs, &
                          rct, rcq, rcsoil, soilw, soilm, smcwlt, smcdry, smcref, &
                          smcmax, cmm, chh, qss, evap, hflx, mask, guess)
        mask = .false.
        where (guess) mask = .true.
        guess = .false.
    enddo

    used = send_data(id_vegtype,   real(vegtype),   Time, mask=lland)
    used = send_data(id_soiltype,  real(soiltype),  Time, mask=lland)
    used = send_data(id_slopetype, real(slopetype), Time, mask=lland)
    used = send_data(id_nroot,     real(nroot),     Time, mask=lland)
    used = send_data(id_canopy,    canopy,    Time, mask=lland)
    used = send_data(id_tslnd,     tslnd,     Time, mask=lland)
    used = send_data(id_sheleg,    sheleg,    Time, mask=lland)
    used = send_data(id_stc,       stc,       Time, mask=lland3d)
    used = send_data(id_smc,       smc,       Time, mask=lland3d)
    used = send_data(id_slc,       slc,       Time, mask=lland3d)
    used = send_data(id_vegfrac,   vegfrac,   Time, mask=lland)
    used = send_data(id_snwdph,    snwdph,    Time, mask=lland)
    used = send_data(id_eta,       eta,       Time, mask=lland)
    used = send_data(id_sheat,     sheat,     Time, mask=lland)
    used = send_data(id_ec,        ec,        Time, mask=lland)
    used = send_data(id_edir,      edir,      Time, mask=lland)
    used = send_data(id_et,        et,        Time, mask=lland)
    used = send_data(id_trans,     trans,     Time, mask=lland)
    used = send_data(id_esnow,     esnow,     Time, mask=lland)
    used = send_data(id_drip,      drip,      Time, mask=lland)
    used = send_data(id_dew,       dew,       Time, mask=lland)
    used = send_data(id_beta,      beta,      Time, mask=lland)
    used = send_data(id_etp,       etp,       Time, mask=lland)
    used = send_data(id_ssoil,     ssoil,     Time, mask=lland)
    used = send_data(id_flx1,      flx1,      Time, mask=lland)
    used = send_data(id_flx2,      flx2,      Time, mask=lland)
    used = send_data(id_flx3,      flx3,      Time, mask=lland)
    used = send_data(id_runoff1,   runoff1,   Time, mask=lland)
    used = send_data(id_runoff2,   runoff2,   Time, mask=lland)
    used = send_data(id_runoff3,   runoff3,   Time, mask=lland)
    used = send_data(id_snomlt,    snomlt,    Time, mask=lland)
    used = send_data(id_sncovr,    sncovr,    Time, mask=lland)
    used = send_data(id_rc,        rc,        Time, mask=lland)
    used = send_data(id_pc,        pc,        Time, mask=lland)
    used = send_data(id_rsmin,     rsmin,     Time, mask=lland)
    used = send_data(id_xlai,      xlai,      Time, mask=lland)
    used = send_data(id_rcs,       rcs,       Time, mask=lland)
    used = send_data(id_rct,       rct,       Time, mask=lland)
    used = send_data(id_rcq,       rcq,       Time, mask=lland)
    used = send_data(id_rcsoil,    rcsoil,    Time, mask=lland)
    used = send_data(id_soilw,     soilw,     Time, mask=lland)
    used = send_data(id_soilm,     soilm,     Time, mask=lland)
    used = send_data(id_smcwlt,    smcwlt,    Time, mask=lland)
    used = send_data(id_smcdry,    smcdry,    Time, mask=lland)
    used = send_data(id_smcref,    smcref,    Time, mask=lland)
    used = send_data(id_smcmax,    smcmax,    Time, mask=lland)

end subroutine do_land


!--------------------------------------------------------------------------------   
subroutine sfc_land_drv(imax, dt, fprcp, swdn, swnet, lwdn, sfcems, ps, p1, t1, &
                        prslki, sfcspd, lprcp, q1, tbot, cd, cdq, vegtyp, soiltyp,  &
                        slopetyp, cmc, sfctmp, sneqv, stsoil, smsoil, sh2o, shdfac, snowh, &
                        nroot, eta, sheat, ec, edir, et, ett, esnow, drip, dew, &
                        beta, etp, ssoil, flx1, flx2, flx3, runoff1, runoff2, &
                        runoff3, snomlt, snowc, rc, pc, rsmin, xlai, rcs, &
                        rct, rcq, rcsoil, soilw, soilm, smcwlt, smcdry, smcref, &
                        smcmax, cmm, chh, qss, evap, hflx, mask, guess)
!--------------------------------------------------------------------------------   
    use sflx_mod, only : sflx
    implicit none
    integer, intent(in) :: imax
    real, intent(in) :: dt
    real, dimension(imax), intent(in) :: fprcp, lprcp !precip rate (kg m-2 s-1)
    real, dimension(imax), intent(in) :: swdn, swnet, lwdn ! Radiation Wm-2
    real, dimension(imax), intent(in) :: sfcems !Sfc Emissitivity
    real, dimension(imax), intent(in) :: ps, p1 ! Pressure in cb
    real, dimension(imax), intent(in) :: t1, sfcspd, q1, cd, cdq, prslki

    integer, dimension(imax), intent(in) :: vegtyp, soiltyp, slopetyp
    real, dimension(imax), intent(inout) :: cmc, sfctmp, sneqv, shdfac, snowh, tbot
    real, dimension(lsoil,imax), intent(inout) :: stsoil, smsoil, sh2o
    integer, dimension(imax), intent(out) :: nroot
    real, dimension(imax), intent(out) :: eta, sheat, ec, &
                             edir, et, ett, esnow, drip, dew, &
                             beta, etp, ssoil, flx1, flx2, flx3, &
                             runoff1, runoff2, runoff3, snomlt, &
                             snowc, rc, pc, rsmin, xlai, rcs, &
                             rct, rcq, rcsoil, soilw, soilm, smcwlt, &
                             smcdry, smcref, smcmax, cmm, chh, qss, evap, hflx

    logical, dimension(imax), intent(in) :: mask, guess

    real, dimension(lsoil,imax) :: stsoil_old, smsoil_old, sh2o_old
    real, dimension(imax) :: snowh_old, sneqv_old, cmc_old 
    real :: zlvl, sfcprs, prs1, q1sat, dqsdt1, theta1, ch, cm, tv1, rho, rch

    integer, parameter :: couple = 1, icein = 0
    real, parameter :: cb2pa = 1000., mm2m = 0.001, m2mm = 1000.
    integer :: i

  
    cmc = cmc * mm2m 
    snowh = snowh * mm2m 
    sneqv = sneqv * mm2m 

    do i = 1, imax
        if (.not.guess(i)) cycle
            snowh_old(i)     =  snowh(i)
            sneqv_old(i)     =  sneqv(i)
            cmc_old(i)       =  cmc(i)
            stsoil_old(:,i)  =  stsoil(:,i)
            smsoil_old(:,i)  =  smsoil(:,i)
            sh2o_old(:,i)    =  sh2o(:,i)
    enddo

    do i = 1, imax
        if (.not.mask(i)) cycle

        sfcprs = ps(i) * cb2pa 
        prs1 = p1(i) * cb2pa

        theta1 = t1(i) * prslki(i) 
        tv1 = t1(i) * (1.0 + con_fvirt*q1(i))

        rho = prs1 / (RDGAS * tv1)

        call compute_qs(t1(i),prs1,q1sat,dqsdT=dqsdt1)

        zlvl = -RDGAS * tv1 * log(prs1/sfcprs) / GRAV 

        ch = cdq(i) * sfcspd(i)
        cm = cd(i) * sfcspd(i)
        cmm(i) = cm
        chh(i) = rho * ch

        if (sneqv(i) /= 0.0 .and. snowh(i) == 0.0) then
            snowh(i) = 10.0 * sneqv(i)
        endif

        call sflx(lsoil, couple, icein, fprcp(i), dt, zlvl, sldpth, &
              swdn(i), swnet(i), lwdn(i), sfcems(i), sfcprs, t1(i), &
              sfcspd(i), lprcp(i), q1(i), q1sat, dqsdt1, theta1, &
              vegtyp(i), soiltyp(i), slopetyp(i), tbot(i), &
              !  ---  input/outputs: &
              cmc(i), sfctmp(i), stsoil(:,i), smsoil(:,i), sh2o(:,i), sneqv(i), ch, cm, &
              !  ---  outputs: &
              nroot(i), shdfac(i), snowh(i), eta(i), sheat(i), ec(i), edir(i), &
              et(i), ett(i), esnow(i), drip(i), dew(i), beta(i), etp(i), ssoil(i), &
              flx1(i), flx2(i), flx3(i), runoff1(i), runoff2(i), runoff3(i), &
              snomlt(i), snowc(i), rc(i), pc(i), rsmin(i), xlai(i), rcs(i), rct(i), & 
              rcq(i), rcsoil(i), soilw(i), soilm(i), smcwlt(i), smcdry(i), smcref(i), &
              smcmax(i))
   
        rch = rho * CP_AIR * cdq(i) * sfcspd(i) 
        qss(i) = q1(i) + eta(i)/(ELOCP*rch)
        evap(i) = eta(i) * HVAPI / rho
        hflx(i) = sheat(i) * CPINV / rho
    end do 

    do i = 1, imax
        if (.not.guess(i)) cycle
            snowh(i)     =  snowh_old(i)
            sneqv(i)     =  sneqv_old(i)
            cmc(i)       =  cmc_old(i)
            stsoil(:,i)  =  stsoil_old(:,i)
            smsoil(:,i)  =  smsoil_old(:,i)
            sh2o(:,i)    =  sh2o_old(:,i)
    enddo

    cmc = cmc * m2mm 
    snowh = snowh * m2mm
    sneqv = sneqv * m2mm

    return
end subroutine sfc_land_drv


!--------------------------------------------------------------------------------   
subroutine get_land_frac(frac)
!--------------------------------------------------------------------------------   
    real, intent(out) :: frac(:,:)

    if (.not.initialized) call mpp_error(FATAL,'sfc_mod: not initialized')
   
    frac = fland

end subroutine get_land_frac

!--------------------------------------------------------------------------------   
subroutine set_surface(Time,tskin,coszen,sfcalb,sfcemis)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), optional  :: coszen(:,:)
    real, intent(out) :: tskin(:,:)
    real, intent(out), optional :: sfcalb(:,:,:)
    real, intent(out), optional :: sfcemis(:,:)

    logical :: ov
    integer :: k

    if (any(lland)) then
        call data_override('ATM','zorl',zorl,Time,ov)
        if (.not.ov) call mpp_error(FATAL,'set_surface: data_override failed for zorl !')
    endif
    where(zorl<zorl_min) zorl = zorl_min

    call data_override('ATM','sst',sst,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'set_surface: data_override failed for sst !')


    call data_override('ATM','fice',fice,Time,ov)
    if (.not.ov) call mpp_error(FATAL,'set_surface: data_override failed for fice !')
    if (ov) where(fice<fice_min) fice = 0.

    call get_tskin(tskin)
  
    if (present(sfcalb).and..not.present(coszen)) &
       call mpp_error(FATAL, 'set_surface: sfcalb present but coszen not present')

    if (present(sfcalb)) call get_albedo(Time,coszen,sfcalb)

    if (present(sfcemis)) call get_emis(sfcemis)
     
    return
end subroutine set_surface
    

!--------------------------------------------------------------------------------   
subroutine get_tskin(tskin)
!--------------------------------------------------------------------------------   
    real, intent(out) :: tskin(:,:)
    logical :: ov

    if (.not.initialized) call mpp_error(FATAL,'sfc_mod: not initialized')

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
    logical :: ov, mask(js:je,is:ie)
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

    mask = .false.
    mask = fice > 0.
    call setalb_sice(imax,mask,hsnow_sice,hice,tsice,albsice)

    mask = .false.
    mask = (focn>0.).and.((1.-fice)>0.)
    call setalb_ocean(imax,mask,coszen,albocn)

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

