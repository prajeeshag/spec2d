module phys_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error, mpp_npes, &
        mpp_get_current_pelist, mpp_pe, mpp_exit, mpp_clock_id, &
        mpp_clock_begin, mpp_clock_end, mpp_sync, mpp_root_pe, &
        mpp_broadcast, mpp_gather, mpp_declare_pelist, &
        mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain, &
        mpp_get_global_domain, mpp_global_field

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, file_exist

use fms_io_mod, only : restart_file_type, reg_rf => register_restart_field, restore_state, save_restart, &
    fms_io_exit

use constants_mod, only : PI, CP_AIR, RDGAS, RVGAS, GRAV, RADIUS, HLF

use time_manager_mod, only : time_type, set_time, operator(==), operator(+), assignment(=), &
                             print_date

use diag_manager_mod, only : diag_axis_init, reg_df=>register_diag_field, send_data, diag_manager_end

use tracer_manager_mod, only : get_tracer_index, get_tracer_name, get_number_tracers

use field_manager_mod, only : MODEL_ATMOS

use radiation_mod, only : init_radiation, con_solr, radiation, NF_ALBD

use sfc_mod, only : init_sfc, get_land_frac, set_surface, do_surface, end_sfc, save_restart_sfc

use vertdiff_mod, only : do_vertical_diffusion

use astronomy_mod, only : astronomy_init, diurnal_solar

use vertical_levels_mod, only : get_pressure_at_levels

use gwdrag_mod, only : init_gwdrag, gwdrag

use gwdrag_conv_mod, only : gwdrag_conv

use cu_conv_mod, only : init_cu_conv, cu_conv

use shallow_conv_mod, only : shallow_conv

use micro_phys_mod, only : init_micro_phys, micro_phys, end_micro_phys, save_restart_micro_phys

use ocpack_mod, only : ocpack_typeP, oc_ny, oc_npack, oc_nlat, get_ocpackP

implicit none
private

public :: init_phys, phys, end_phys, save_restart_phys

real, parameter :: rhoh2o = 1000.

type(domain2D), pointer :: domain
integer :: ntrac, nlev
integer :: is, ie, ilen, js, je, jlen
integer :: ind_q = 0, ind_clw = 0, ind_cli = 0

integer :: dt_rad
real :: dt_phys, rdt_phys, dt_atmos, rdt_atmos
type(time_type) :: time_step_rad, rad_time, time_step

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lon_deg, lat_deg, gdlen
real, allocatable, dimension(:,:,:) :: htsw
real, allocatable, dimension(:,:,:) :: htlw
real, allocatable, dimension(:,:) :: rldsz
real, allocatable, dimension(:,:) :: rsdsz, rsusz
real, allocatable, dimension(:,:) :: rsdt, rsut
real, allocatable, dimension(:,:) :: rlut
real, allocatable, dimension(:,:) :: fprcp, lprcp
real, allocatable, dimension(:,:) :: slmsk

character(len=16) :: resfnm = 'amfi_phys_res'
character (len=8) :: rou='am_phys'
type(restart_file_type) :: rstrt

integer :: id_rsds, id_rsus, id_rsns, id_dtlw, id_dtrd, id_shflx, id_lhflx, id_taux, &
           id_tauy, id_dtvd, id_hpbl, id_duvd, id_dvvd, id_rlds, id_rlus, id_tskin, &
           id_sfcemis, id_snlhflx, id_totene, id_rsdt, id_rsut, id_rlut
integer :: id_dugwd, id_dvgwd, id_ducgwd, id_dvcgwd
integer :: id_dtcu, id_ducu, id_dvcu, id_lprcu, id_kcnv, id_fprcu
integer :: id_dtsc, id_dtmp, id_lprmp, id_fprmp, id_inmpclw, id_inmpq
integer :: id_dtphy, id_lpr, id_fpr, id_pr, id_duphy, id_dvphy, id_ua, id_va, id_ta
           

integer, allocatable, dimension(:) :: id_dtrvd, id_dtrmp, id_dtrphy, id_dtrcu, id_dtrsc, id_trin

integer :: clck_phys, clck_rad, clck_sfc, clck_vd, clck_gwd, clck_cu, clck_sc, clck_mp, clck_gwdc

logical :: initialized=.false., debug=.false.

namelist/phys_nml/dt_rad,debug

contains

!--------------------------------------------------------------------------------   
subroutine init_phys(Time, dt_phys_in, dt_atmos_in, domain_in, nlev_in, &
                     lat_deg_in, lon_deg_in, axes)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2d), target :: domain_in
    real, intent(in) :: dt_phys_in, dt_atmos_in 
    integer, intent(in) :: nlev_in
    real, intent(in) :: lat_deg_in(:,:), lon_deg_in(:,:)
    integer, intent(in) :: axes(4)
    
    real, allocatable :: fland(:,:), ilevs(:), ilevsp(:)
    integer :: unit, i, j, nlon
    integer :: indx, num_prog, num_diag
    real :: tem1, tem2
    type(ocpack_typeP), allocatable :: ocpkP(:,:)

    unit = open_namelist_file()

    read(unit,nml=phys_nml)

    call close_file(unit)

    dt_phys = dt_phys_in
    rdt_phys = 1./dt_phys

    dt_atmos = dt_atmos_in
    rdt_atmos = 1./dt_atmos

    domain => domain_in
    nlev = nlev_in

    if (mod(real(dt_rad),dt_atmos)/=0.) &
             call mpp_error('phys_mod', 'dt for radiation ' &
                  //'should be a multiple of model timestep', FATAL)

    time_step_rad = set_time(dt_rad)
    time_step = set_time(int(dt_rad))
    rad_time = Time + set_time(int(dt_atmos))

    call get_number_tracers(MODEL_ATMOS, ntrac, num_prog, num_diag)

    ind_q = get_tracer_index(MODEL_ATMOS, 'sphum')   
    ind_clw = get_tracer_index(MODEL_ATMOS, 'clw')
    ind_cli = get_tracer_index(MODEL_ATMOS, 'cli')
 
    call mpp_get_compute_domain(domain,js,je,is,ie)

    jlen = je - js + 1
    ilen = ie - is + 1

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))
    allocate(gdlen(js:je,is:ie))

    lat_deg(js:je,is:ie) = lat_deg_in(js:je,is:ie)
    lat_rad(js:je,is:ie) = lat_deg_in(js:je,is:ie)*PI/180.

    lon_deg(js:je,is:ie) = lon_deg_in(js:je,is:ie)
    lon_rad(js:je,is:ie) = lon_deg_in(js:je,is:ie)*PI/180.

    allocate(ocpkP(oc_npack(),oc_ny()))
    call get_ocpackP(ocpkP) 

    do i = is, ie
        do j = js, je
            nlon = ocpkP(1,j)%ilen
            if (i>nlon) nlon = ocpkP(2,j)%ilen
            tem1 = 2. * PI * RADIUS * cos(lat_rad(j,i)) / nlon 
            tem2 = PI * RADIUS / oc_nlat()
            gdlen(j,i) = sqrt(tem1**2 + tem2**2)
        enddo
    enddo

    allocate(htsw(nlev,js:je,is:ie))
    allocate(htlw(nlev,js:je,is:ie))
    allocate(rsdsz(js:je,is:ie))
    allocate(rsusz(js:je,is:ie))
    allocate(rldsz(js:je,is:ie))
    allocate(rsdt(js:je,is:ie))
    allocate(rsut(js:je,is:ie))
    allocate(rlut(js:je,is:ie))
    allocate(fprcp(js:je,is:ie))
    allocate(lprcp(js:je,is:ie))

    htsw  = 0. 
    htlw  = 0.
    rsdsz = 0.
    rsusz = 0.
    rldsz = 0.
    rsut = 0.
    rsdt = 0.
    rlut = 0.

    fprcp = 0.
    indx = reg_rf(rstrt, resfnm, 'fprcp', fprcp, domain, mandatory=.false.)

    lprcp = 0.
    indx = reg_rf(rstrt, resfnm, 'lprcp', lprcp, domain, mandatory=.false.)

    call init_sfc(Time,dt_atmos,domain,axes(1:2),lat_deg_in)
    allocate(fland(js:je,is:ie))
    allocate(slmsk(js:je,is:ie))
    call get_land_frac(fland)
    slmsk = 0.
    where(fland>=0.5) slmsk = 1.

    call astronomy_init()
    call init_radiation(Time, domain, ntrac, nlev, lat_deg_in, &
                        lon_deg_in, fland, axes, ind_q_in=ind_q, ind_clw_in=ind_clw) 

    call init_gwdrag(domain)

    call init_cu_conv(is,ie,js,je)

    call init_micro_phys(domain, nlev, lat_deg_in*PI/180.) 

    deallocate(fland)

    clck_phys = mpp_clock_id('Physics')
    clck_rad  = mpp_clock_id(':---> Radiation')
    clck_sfc  = mpp_clock_id(':---> Surface')
    clck_vd   = mpp_clock_id(':---> Vertical Diffusion')
    clck_gwd  = mpp_clock_id(':---> GW Drag')
    clck_gwdc = mpp_clock_id(':---> CGW Drag')
    clck_cu   = mpp_clock_id(':---> Cumulus Convection')
    clck_sc   = mpp_clock_id(':---> Shallow Convection')
    clck_mp   = mpp_clock_id(':---> Micro Physics')

    call init_diag_out(Time, axes)   

    if (file_exist('INPUT/'//trim(resfnm))) call restore_state(rstrt)

    initialized = .true.

end subroutine init_phys



!--------------------------------------------------------------------------------   
subroutine init_diag_out(Time, axes)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer, intent(in) :: axes(4)
    integer :: n
    character(len=256) :: longname, units, fldnm, name
 
    allocate(id_dtrvd(ntrac), id_dtrsc(ntrac), &
             id_dtrmp(ntrac), id_dtrcu(ntrac), &
             id_dtrphy(ntrac), id_trin(ntrac))

    ! Diag out
    ! -------------------------------------------------------------------------------- 
    !->rd 
    id_rsdt = reg_df(rou, 'rsdt', axes(1:2), time, 'TOA downwelling shortwave',  'w m-2')

    id_rsut = reg_df(rou, 'rsut', axes(1:2), time, 'TOA upwelling shortwave',  'w m-2')

    id_rlut = reg_df(rou, 'rlut', axes(1:2), time, 'TOA upwelling longtwave',  'w m-2')

    id_rsds = reg_df(rou, 'rsds', axes(1:2), time, 'surface downwelling shortwave',  'w m-2')

    id_rsus = reg_df(rou, 'rsus', axes(1:2), time, 'surface upwelling shortwave',  'w m-2')

    id_rsns = reg_df(rou, 'rsns', axes(1:2), time, 'surface net shortwave',  'w m-2')

    id_sfcemis = reg_df(rou, 'sfcemis', axes(1:2), time, 'surface emisivity',  'w m-2')

    id_rlds = reg_df(rou, 'rlds', axes(1:2), time, 'surface downwelling longwave',  'w m-2')
    id_rlus = reg_df(rou, 'rlus', axes(1:2), time, 'surface upwelling longwave',  'w m-2')
    id_totene = reg_df(rou, 'totene', axes(1:2), time, 'total energy flux',  'w m-2')
    id_tskin = reg_df(rou, 'tskin', axes(1:2), time, 'Surface Temperature',  'K')

    id_dtrd = reg_df(rou, 'dtrd', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (radiation)',  'k s-1')

    id_dtlw = reg_df(rou, 'dtlw', [axes(3),axes(1),axes(2)], time, &
               'temperature tendency (lw)',  'ks-1')

    !<-rd

    !->vd
    id_shflx = reg_df(rou, 'shflx', [axes(1),axes(2)], time, &
               'sensible heat flux', 'wm-2')
    id_lhflx = reg_df(rou, 'lhflx', [axes(1),axes(2)], time, &
               'latent heat flux', 'wm-2')
    id_snlhflx = reg_df(rou, 'snlhflx', [axes(1),axes(2)], time, &
               'heat released due snow formation)', 'wm-2')
    id_taux = reg_df(rou, 'taux', [axes(1),axes(2)], time, &
               'u-wind stress', 'pa')
    id_tauy = reg_df(rou, 'tauy', [axes(1),axes(2)], time, &
               'v-wind stress', 'pa')
    id_hpbl = reg_df(rou, 'hpbl', [axes(1),axes(2)], time, &
               'height of planetary boundary layer', 'm')
    id_dtvd = reg_df(rou, 'dtvd', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (vertical diffusion)', 'k s-1')

    id_duvd = reg_df(rou, 'duvd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (vertical diffusion)', 'm/s-2')
    id_dvvd = reg_df(rou, 'dvvd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (vertical diffusion)', 'm/s-2')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'vd'
        longname = trim(longname)//' tendency (vertical diffusion)'
        units = trim(units)//' s-1'
        id_dtrvd(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo
    !<-vd
    
    !-> gwd
    id_dugwd = reg_df(rou, 'dugwd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (gravity wave drag)', 'm/s-2')
    id_dvgwd = reg_df(rou, 'dvgwd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (gravity wave drag)', 'm/s-2')
    !<-gwd

    !->gwdc
    id_ducgwd = reg_df(rou, 'ducgwd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (convective gravity wave drag)', 'm/s-2')
    id_dvcgwd = reg_df(rou, 'dvcgwd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (convective gravity wave drag)', 'm/s-2')
    !<-gwdc

    !->cu
    id_dtcu = reg_df(rou, 'dtcu', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (cumulus conv)', 'k s-1')
    id_ducu = reg_df(rou, 'ducu', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (cumulus conv)', 'm/s-2')
    id_dvcu = reg_df(rou, 'dvcu', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (cumulus conv)', 'm/s-2')
    id_lprcu = reg_df(rou, 'lprcu', [axes(1),axes(2)], time, &
                'rainfall rate (cumulus conv)', 'kg m-2 s-1')
    id_fprcu = reg_df(rou, 'fprcu', [axes(1),axes(2)], time, &
                'snowfall rate (cumulus conv)', 'kg m-2 s-1')
    id_kcnv = reg_df(rou, 'kcnv', [axes(1),axes(2)], time, &
                'cumulus conv occurence', '1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'cu'
        longname = trim(longname)//' tendency (cumulus convection)'
        units = trim(units)//' s-1'
        id_dtrcu(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo
    !<-cu

    !->sc
    id_dtsc = reg_df(rou, 'dtsc', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (shallow conv)', 'k s-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'sc'
        longname = trim(longname)//' tendency (shallow convection)'
        units = trim(units)//' s-1'
        id_dtrsc(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo
    !<-sc

    !->mp
    id_inmpclw = reg_df(rou, 'inmpclw', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (micro phys)', 'k s-1')
    id_inmpq = reg_df(rou, 'inmpq', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (micro phys)', 'k s-1')
    id_dtmp = reg_df(rou, 'dtmp', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (micro phys)', 'k s-1')
    id_lprmp = reg_df(rou, 'lprmp', [axes(1),axes(2)], time, &
                'rainfall rate (micro phys)', 'kg m-2 s-1')
    id_fprmp = reg_df(rou, 'fprmp', [axes(1),axes(2)], time, &
                'snowfall rate (micro phys)', 'kg m-2 s-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'mp'
        longname = trim(longname)//' tendency (micro phys)'
        units = trim(units)//' s-1'
        id_dtrmp(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo
    !<-mp

    !->phy
    id_ua = reg_df(rou, 'ua', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity (phys-in)', 'ms-1')
    id_va = reg_df(rou, 'va', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity (phys-in)', 'ms-1')
    id_ta = reg_df(rou, 'ta', [axes(3),axes(1),axes(2)], time,  &
               'Temperature', 'K')
    id_dtphy = reg_df(rou, 'dtphy', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (phys)', 'k s-1')
    id_duphy = reg_df(rou, 'duphy', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (phys)', 'm/s-2')
    id_dvphy = reg_df(rou, 'dvphy', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (phys)', 'm/s-2')
    id_lpr = reg_df(rou, 'lpr', [axes(1),axes(2)], time, &
               'rainfall rate', 'kg m-2 s-1')
    id_fpr = reg_df(rou, 'fpr', [axes(1),axes(2)], time, &
                'snowfall rate', 'kg m-2 s-1')
    id_pr = reg_df(rou, 'pr', [axes(1),axes(2)], time, &
               'precipitation rate', 'kg m-2 s-1')
    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = 'd'//trim(name)//'phy'
        longname = trim(longname)//' tendency (phys)'
        units = trim(units)//' s-1'
        id_dtrphy(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo

    do n = 1, ntrac
        if (.not.get_tracer_name(MODEL_ATMOS,n,name,longname,units)) &
            call mpp_error(FATAL,'phys_mod: get_tracer_name failed')
        fldnm = trim(name)//'_in'
        longname = trim(longname)//' (phys)'
        id_trin(n) = reg_df(rou, fldnm, [axes(3),axes(1),axes(2)], time,  &
                        longname, units)
    enddo
    !<-phy
    
    return
end subroutine init_diag_out
    

!--------------------------------------------------------------------------------   
subroutine phys(Time,tlyr1,tr1,p1,u1,v1,vvel1,dtdt,dqdt,dudt,dvdt,topo,enetot)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(1:nlev,js:je,is:ie) :: tlyr1, u1, v1, vvel1
    real, intent(in), dimension(1:nlev,js:je,is:ie,1:ntrac) :: tr1
    real, intent(in), dimension(js:je,is:ie) :: p1
    real, intent(out), dimension(1:nlev,js:je,is:ie) :: dtdt, dudt, dvdt 
    real, intent(out), dimension(1:nlev,js:je,is:ie,ntrac) :: dqdt
    real, intent(in), optional, dimension(js:je,is:ie) :: topo
    real, intent(out), optional, dimension(js:je,is:ie) :: enetot

    real, dimension(1:nlev,js:je,is:ie) :: plyr, plyrk, prslki, delp, phil, tmp3d
    real, dimension(1:nlev+1,js:je,is:ie) :: plvl, plvlk, phii
    real, dimension(1:nlev,js:je,is:ie) :: dtdt1, dudt1, dvdt1, clw, cli
    real, dimension(1:nlev,js:je,is:ie,ntrac) :: dqdt1

    real, dimension(js:je,is:ie) :: tskin, fracday, coszen, rcoszen, rsds, rsus, rsns, &
                                    rlds, rlus, rldsg, rb, ffmm, ffhh, qss, hflx, evap, &
                                    stress, wind, dusfc1, dvsfc1, dtsfc1, dqsfc1, hpbl, &
                                    gamt, gamq, topo1, sfcemis, cldwrk, rain1, rain, &
                                    snow, snow1, totene, snlhflx

    integer, dimension(js:je,is:ie) :: kpbl, kbot, ktop, kcnv

    real, dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real :: solcon, rrsun
    integer :: k, imax, n
    logical :: used, chck(js:je,is:ie)
    character(len=32) :: fnm


    call mpp_clock_begin(clck_phys)

    imax = size(p1,1)*size(p1,2)

    dtdt = 0.; dudt = 0.; dvdt = 0.; dqdt = 0.
    dtdt1 = 0.; dudt1 = 0.; dvdt1 = 0.; dqdt1 = 0.
    rain = 0.; snow = 0.; rain1 = 0.; snow1 = 0.
    totene = 0.

    topo1 = 0.
    if(present(topo)) topo1 = topo

    if (id_ua>0) used = send_data(id_ua, u1, Time)
    if (id_va>0) used = send_data(id_va, v1, Time)
    if (id_ta>0) used = send_data(id_ta, tlyr1, Time)

    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'------begin phys-----')
        if (any(tlyr1<173.15).or.any(tlyr1>373.15)) then
            fnm=""
            write(fnm,*)mpp_pe()
            fnm="debug_"//trim(adjustl(fnm))
            call write_data(fnm,'tlyr1',tlyr1)
            call fms_io_exit()
            call mpp_error(FATAL,"temperature out of range")
        endif
    endif

    call get_pressure_at_levels(p1,plvl,plyr,plvlk,plyrk)
    prslki = plvlk(1:nlev,:,:)/plyrk
    delp = plvl(1:nlev,:,:) - plvl(2:nlev+1,:,:)

    call get_phi(topo1, tlyr1, tr1(:,:,:,ind_q), plvl, plvlk, plyr, plyrk, phii, phil)

    !Radiation
    !--------------------------------------------------------------------------------   
    call mpp_clock_begin(clck_rad)
    if (rad_time==Time) then
        call print_date(Time,'Calling Radiation at:')
        coszen(:,:) = 0.0
        rcoszen(:,:) = 0.0
        call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step_rad)
        call set_surface(Time,tskin,coszen,sfcalb,sfcemis)
        if (debug) then
            call mpp_sync()
            call mpp_error(NOTE,'after set_surface')
        endif
        solcon = con_solr * rrsun
        where(coszen>0.) rcoszen = 1./coszen 
        call radiation(Time, tlyr1, tr1, plyr, plvl, phil, topo1, tskin, coszen, fracday, sfcalb, sfcemis, &
                       solcon, htsw, rsdsz, rsusz, htlw, rldsz, rlus, rsdt, rsut, rlut)
        do k = 1, size(htsw,1)
            htsw(k,:,:) = htsw(k,:,:) * rcoszen
        enddo
        rsdsz = rsdsz * rcoszen
        rsusz = rsusz * rcoszen
        rldsz = rldsz / (tlyr1(1,:,:)**4.)
        rad_time = rad_time + time_step_rad
    else
        call set_surface(Time,tskin,sfcemis=sfcemis)
        if (debug) then
            call mpp_sync()
            call mpp_error(NOTE,'after set_surface')
        endif
    endif
    call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step)
    coszen = coszen * fracday

    call adjust_rad(coszen, tskin, tlyr1(1,:,:), rsdsz, rsusz, rldsz, htsw,  &
                    htlw, sfcemis, dtdt1, rsds, rsus, rlds, rlus)

    dtdt = dtdt + dtdt1

    rsns = rsds - rsus

    rldsg = rlds * sfcemis !-> rldsg for land,ocn,sea-ice

    totene = totene - rsds + rsus - rlds + rlus + rsdt - rsut - rlut

    used = send_data(id_rsdt, rsdt, Time) 
    used = send_data(id_rsut, rsut, Time) 
    used = send_data(id_rlut, rlut, Time) 
    used = send_data(id_rsds, rsds, Time) 
    used = send_data(id_rsus, rsus, Time) 
    used = send_data(id_rsns, rsns, Time) 
    used = send_data(id_rlds, rlds, Time) 
    used = send_data(id_rlus, rlus, Time) 
    used = send_data(id_tskin, tskin, Time) 
    used = send_data(id_sfcemis, sfcemis, Time) 
    used = send_data(id_dtrd, dtdt1, Time) 
    used = send_data(id_dtlw, htlw, Time) 

    call mpp_clock_end(clck_rad)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after Radiation')
    endif
    !End Radiation

    !Surface Fluxes
    call mpp_clock_begin(clck_sfc)
    call do_surface(Time, p1(:,:), u1(1,:,:), v1(1,:,:), tlyr1(1,:,:), &
            tr1(1,:,:,ind_q), plyr(1,:,:), prslki(1,:,:), rsds, rsns, rldsg, &
            fprcp, lprcp, rb, ffmm, ffhh, qss, hflx, evap, stress, wind)
    call mpp_clock_end(clck_sfc)

    if (debug) then
        call mpp_sync()
        !write(msg,*) tlyr1(1,:,:)
        call mpp_error(NOTE,'after sfc_fluxes')
    endif
    !End surface Fluxes
    !Vertical Diffusion
    call mpp_clock_begin(clck_vd)
    dtdt1 = 0.; dudt1 = 0.; dvdt1 = 0.; dqdt1 = 0.
    dusfc1 = 0.; dvsfc1 = 0.; dtsfc1 = 0.; dqsfc1 = 0.
    call do_vertical_diffusion(imax, nlev, ntrac, dvdt1, dudt1, dtdt1, dqdt1, u1, &
                    v1, tlyr1, tr1, plvlk, rb, ffmm, ffhh, qss, hflx, evap, stress, &
                    wind, kpbl, plvl, delp, plyr, plyrk, phii, phil, dt_phys, dusfc1, &
                    dvsfc1, dtsfc1, dqsfc1, hpbl, gamt, gamq)

    used = send_data(id_shflx, dtsfc1, Time)
    used = send_data(id_lhflx, dqsfc1, Time)
    used = send_data(id_hpbl, hpbl, Time)
    used = send_data(id_taux, dusfc1, Time)
    used = send_data(id_tauy, dvsfc1, Time)
    used = send_data(id_dtvd, dtdt1, Time)
    used = send_data(id_duvd, dudt1, Time)
    used = send_data(id_dvvd, dvdt1, Time)
    do n = 1, ntrac
        used = send_data(id_dtrvd(n), dqdt1(:,:,:,n), Time)
    enddo

    dtdt = dtdt + dtdt1 
    dudt = dudt + dudt1 
    dvdt = dvdt + dvdt1 
    dqdt = dqdt + dqdt1 
    call mpp_clock_end(clck_vd)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after vertical_diffusion')
    endif
    !End Vertical Diffusion
    totene = totene + dtsfc1 + dqsfc1

    !Gravity Wave Drag
    call mpp_clock_begin(clck_gwd)
    dudt1 = 0.; dvdt1 = 0.
    call gwdrag(dvdt1, dudt1, u1, v1, tlyr1, tr1(:,:,:,ind_q), kpbl, &
                plvl, delp, plyr, plyrk, phii, phil, dt_phys, dusfc1, dvsfc1)

    used = send_data(id_dugwd, dudt1, Time)
    used = send_data(id_dvgwd, dvdt1, Time)

    dudt = dudt + dudt1 
    dvdt = dvdt + dvdt1 
    call mpp_clock_end(clck_gwd)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after GW_Drag')
    endif
    !End Gravity Wave Drag

    dtdt = tlyr1 + dtdt*dt_phys 
    dudt = u1 + dudt*dt_phys 
    dvdt = v1 + dvdt*dt_phys 
    dqdt = tr1 + dqdt*dt_phys 
    
    call get_phi(topo1, dtdt, dqdt(:,:,:,ind_q), plvl, plvlk, plyr, plyrk, phii, phil)
    
    !Cumulus convection
    call mpp_clock_begin(clck_cu)
    dtdt1 = dtdt
    dudt1 = dudt
    dvdt1 = dvdt
    dqdt1 = dqdt
    rain1 = 0.
    snow1 = 0.
    kbot = 0
    ktop = 0
    kcnv = 0
    clw = dqdt1(:,:,:,ind_clw)  
    if (ind_cli>0) cli = dqdt1(:,:,:,ind_cli)

    call cu_conv (dt_phys, delp, plyr, p1, phil, clw, cli, dqdt1(:,:,:,ind_q), &
                  dtdt1, dudt1, dvdt1, cldwrk, rain1, kbot, ktop, kcnv, slmsk, vvel1) 

    dqdt1(:,:,:,ind_clw) = clw
    if (ind_cli>0) dqdt1(:,:,:,ind_cli) = cli

    tmp3d = rdt_phys*(dtdt1-dtdt)
    if (id_dtcu>0) used = send_data(id_dtcu, tmp3d, Time)
    do n = 1, ntrac
        if (id_dtrcu(n)>0) used = send_data(id_dtrcu(n), rdt_phys*(dqdt1(:,:,:,n) &
                                            - dqdt(:,:,:,n)), Time)
    enddo
    if (id_ducu>0) used = send_data(id_ducu, rdt_phys*(dudt1-dudt), Time)
    if (id_dvcu>0) used = send_data(id_dvcu, rdt_phys*(dvdt1-dvdt), Time)
    if (id_lprcu>0) used = send_data(id_lprcu, rdt_phys*rain1*rhoh2o, Time)
    if (id_fprcu>0) used = send_data(id_fprcu, rdt_phys*snow1*rhoh2o, Time)
    if (id_kcnv>0) used = send_data(id_kcnv, real(kcnv), Time)

    dtdt = dtdt1 
    dudt = dudt1
    dvdt = dvdt1
    dqdt = dqdt1

    rain = rain + rain1
    snow = snow + snow1
    call mpp_clock_end(clck_cu)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after cumulus_conv')
    endif
    !End Cumulus Convection


    ! Convective Gravity Wave Drag
    call mpp_clock_begin(clck_gwdc)
    dudt1 = 0.
    dvdt1 = 0.
    dusfc1 = 0.
    dvsfc1 = 0.
    call gwdrag_conv(nlev, imax, u1, v1, tlyr1, tr1(:,:,:,ind_q), plyr, delp, &
                     tmp3d, plvl, ktop, kbot, kcnv, dudt1, dvdt1, gdlen, dusfc1, dvsfc1)   
    used = send_data(id_ducgwd, dudt1, Time)
    used = send_data(id_dvcgwd, dvdt1, Time)

    dudt = dudt + dt_phys*dudt1
    dvdt = dvdt + dt_phys*dvdt1
    call mpp_clock_end(clck_gwdc)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after Conv_GW_Drag')
    endif
    !End Convective Gravity Wave Drag


    ! Shallow Convection
    call mpp_clock_begin(clck_sc)
    dtdt1 = dtdt
    dqdt1 = dqdt !Initialize
    call shallow_conv(imax, nlev, dt_phys, delp, plvl, plyr, plyrk, kcnv, &
                      dqdt1(:,:,:,ind_q), dtdt1)
    if (id_dtsc>0) used = send_data(id_dtsc, rdt_phys*(dtdt1-dtdt), Time)
    do n = 1, ntrac
        if (id_dtrsc(n)>0) used = send_data(id_dtrsc(n), rdt_phys*(dqdt1(:,:,:,n) &
                                            - dqdt(:,:,:,n)), Time)
    enddo
    dtdt = dtdt1
    dqdt = dqdt1
    call mpp_clock_end(clck_sc)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after shallow_conv')
    endif
    !End Shallow Convection

    !Micro Physics
    call mpp_clock_begin(clck_mp)
    dtdt1 = dtdt
    dqdt1 = dqdt 

    if (id_inmpclw>0) used = send_data(id_inmpclw, dqdt1(:,:,:,ind_clw), Time)
    if (id_inmpq>0) used = send_data(id_inmpq, dqdt1(:,:,:,ind_q), Time)

    call micro_phys(dt_phys, plyr, p1, plyrk, delp, dqdt1(:,:,:,ind_q), &
                    dqdt1(:,:,:,ind_clw), dtdt1, rain1, snow1)

    if (id_dtmp>0) used = send_data(id_dtmp, rdt_phys*(dtdt1-dtdt), Time)
    do n = 1, ntrac
        if (id_dtrmp(n)>0) used = send_data(id_dtrmp(n), rdt_phys*(dqdt1(:,:,:,n) &
                                                         - dqdt(:,:,:,n)), Time)
    enddo
    if (id_lprmp>0) used = send_data(id_lprmp, rdt_phys*rain1*rhoh2o, Time)
    if (id_fprmp>0) used = send_data(id_fprmp, rdt_phys*snow1*rhoh2o, Time)

    dtdt = dtdt1
    dqdt = dqdt1
    rain = rain + rain1
    snow = snow + snow1
    call mpp_clock_end(clck_mp)
    if (debug) then
        call mpp_sync()
        call mpp_error(NOTE,'after micro_phys')
    endif
    !End Micro Physics

    lprcp = rain * rdt_phys
    fprcp = snow * rdt_phys

    snlhflx(:,:) = snow(:,:)*HLF*rdt_phys*rhoh2o ! Latent heat release due to snow formation (W/m2)

    totene = totene + snlhflx


    if (id_snlhflx>0) used = send_data(id_snlhflx, snlhflx, Time)
    if (id_lpr>0) used = send_data(id_lpr, rdt_phys*rain*rhoh2o, Time)
    if (id_fpr>0) used = send_data(id_fpr, rdt_phys*snow*rhoh2o, Time)
    if (id_pr>0) used = send_data(id_pr, rdt_phys*(rain+snow)*rhoh2o, Time)

    if (id_dtphy>0) used = send_data(id_dtphy, rdt_phys*(dtdt-tlyr1), Time)
    if (id_duphy>0) used = send_data(id_duphy, rdt_phys*(dudt-u1), Time)
    if (id_dvphy>0) used = send_data(id_dvphy, rdt_phys*(dvdt-v1), Time)

    do n = 1, ntrac
        if (id_dtrphy(n)>0) used = send_data(id_dtrphy(n), rdt_phys*(dqdt(:,:,:,n) &
                                                 - tr1(:,:,:,n)), Time)
        if (id_trin(n)>0) used = send_data(id_trin(n), tr1(:,:,:,n), Time)
    enddo

    if (id_totene>0) used = send_data(id_totene, totene, Time)

    if (present(enetot)) enetot=totene

    call mpp_clock_end(clck_phys)

    return
end subroutine phys


!--------------------------------------------------------------------------------   
subroutine get_phi(topo,t,q,prsi,prki,prsl,prkl,phii,phil)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(in), dimension(:,:) :: topo
    real, intent(in), dimension(:,:,:) :: t, q, prsi, prki, prsl, prkl
    real, intent(out), dimension(:,:,:) :: phii, phil

    real, dimension(size(topo,1),size(topo,2)) :: tem, dphib, dphit
    integer :: k
    real, parameter :: fvirt = RVGAS/RDGAS - 1.

    phii(1,:,:) = topo

    do k = 1, size(t,1)
        tem = CP_AIR * t(k,:,:) * (1. + fvirt * q(k,:,:)) / prkl(k,:,:)
        dphib = (prki(k,:,:) - prkl(k,:,:)) * tem 
        dphit = (prkl(k,:,:) - prki(k+1,:,:)) * tem
        phil(k,:,:) = phii(k,:,:) + dphib
        phii(k+1,:,:) = phil(k,:,:) + dphit
    enddo

    return

end subroutine get_phi


! ===================================================================== !
!  subroutine adjust_rad(tskin, t1, sfcdsw, sfcnsw, sfcdlw, swh, hlw, & !
!                dtdt, adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw)      !
!    dcyc2t3 fits radiative fluxes and heating rates from a coarse      !
!    radiation calc time interval into model's more frequent time steps.!
!    solar heating rates and fluxes are scaled by the ratio of cosine   !
!    of zenith angle at the current time to the mean value used in      !
!    radiation calc.  surface downward lw flux is scaled by the ratio   !
!    of current surface air temperature (temp**4) to the corresponding  !
!    temperature saved during lw radiation calculation. upward lw flux  !
!    at the surface is computed by current ground surface temperature.  !
!                                                                       !
!--------------------------------------------------------------------------------   
subroutine adjust_rad(coszen, tskin, t1, sfcdsw, sfcusw, sfcdlw, swh, hlw, emis, &
                      dtdt, rsds, rsus, rlds, rlus)
!--------------------------------------------------------------------------------   
    use constants_mod, only : con_sbc => STEFAN
    implicit none
    real, dimension(:,:), intent(in) :: coszen, tskin, t1, sfcdlw, sfcdsw, sfcusw
    real, dimension(:,:), intent(in) :: emis
    real, dimension(:,:,:), intent(in) :: swh, hlw
    real, dimension(:,:,:), intent(out) :: dtdt
    real, dimension(:,:), intent(out) :: rsds, rsus, rlds, rlus 

    integer :: k
!  --- ...  adjust sfc downward lw flux to account for t changes in layer 1.
!           1st line is for original version, 2nd one is for updated version
    rlds = sfcdlw * (t1**4.)

!  --- ...  compute sfc upward lw flux from current temp,
!      note: sfc emiss effect is not appied at this time
    rlus = con_sbc * (tskin**4.) * emis + (1.-emis) * rlds

!  --- ...  adjust sfc net and downward sw fluxes for zenith angle changes
    rsus = sfcusw * coszen
    rsds = sfcdsw * coszen

!  --- ...  adjust sw heating rates with zenith angle change and
!           add with lw heating to temperature tendency

    do k = 1, size(dtdt,1)
        dtdt(k,:,:) = swh(k,:,:)*coszen + hlw(k,:,:)
    enddo

    return
end subroutine adjust_rad

subroutine save_restart_phys(tstamp)
    character(len=*), optional :: tstamp

    call save_restart_sfc(tstamp)
    call save_restart_micro_phys(tstamp)
    call save_restart(rstrt,tstamp)
    
    return
end subroutine save_restart_phys


subroutine end_phys()

    call end_sfc()
    call end_micro_phys()
    call save_restart(rstrt)

    return
end subroutine end_phys

end module phys_mod

