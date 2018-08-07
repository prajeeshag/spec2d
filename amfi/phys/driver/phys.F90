module phys_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain
use mpp_domains_mod, only : mpp_get_global_domain, mpp_global_field

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use fms_io_mod, only : restart_file_type, reg_rf => register_restart_field
use fms_io_mod, only : restore_state, save_restart

use constants_mod, only : PI, CP_AIR, RDGAS, RVGAS, GRAV, RADIUS

use time_manager_mod, only : time_type, set_time, operator(==), operator(+), assignment(=), &
                             print_date

use diag_manager_mod, only : diag_axis_init, reg_df=>register_diag_field, send_data

use radiation_mod, only : init_radiation, con_solr, radiation, NF_ALBD

use sfc_mod, only : init_sfc, get_land_frac, set_surface, do_surface

use vertdiff_mod, only : do_vertical_diffusion

use astronomy_mod, only : astronomy_init, diurnal_solar

use vertical_levels_mod, only : get_pressure_at_levels

use gwdrag_mod, only : init_gwdrag, gwdrag

use gwdrag_conv_mod, only : gwdrag_conv

use cu_conv_mod, only : init_cu_conv, cu_conv

use shallow_conv_mod, only : shallow_conv

use micro_phys_mod, only : init_micro_phys, micro_phys

implicit none
private

public :: init_phys, phys

real, parameter :: rhoh2o = 1000.

type(domain2D), pointer :: domain
integer :: ntrac, nlev
integer :: is, ie, ilen, js, je, jlen
integer :: nlon, nlat
integer :: ind_q = 1, ind_clw = 2, ind_cli = 0

integer :: dt_rad
real :: dt_phys, rdt_phys, dt_atmos, rdt_atmos
type(time_type) :: time_step_rad, rad_time, time_step

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg, gdlen
real, allocatable, dimension(:,:,:) :: htsw
real, allocatable, dimension(:,:,:) :: htlw
real, allocatable, dimension(:,:) :: rldsz
real, allocatable, dimension(:,:) :: rsdsz, rsusz
real, allocatable, dimension(:,:) :: fprcp, lprcp
real, allocatable, dimension(:,:) :: slmsk

real :: pdryini = 0.

character(len=16) :: resfnm = 'phys_res'
character (len=8) :: rou='am_phys'

integer :: id_rsds, id_rsus, id_rsns, id_dtlw, id_dtrd, id_shflx, id_lhflx, id_taux, &
           id_tauy, id_dtvd, id_hpbl, id_dqvd, id_duvd, id_dvvd, id_rlds, id_rlus, id_tskin, &
           id_sfcemis
integer :: id_dugwd, id_dvgwd, id_ducgwd, id_dvcgwd
integer :: id_dtcu, id_dqcu, id_ducu, id_dvcu, id_lprcu, id_kcnv, id_fprcu
integer :: id_dtsc, id_dqsc
integer :: id_dtmp, id_dqmp, id_lprmp, id_fprmp
integer :: id_dtphy, id_dqphy, id_lpr, id_fpr, id_pr, id_duphy, id_dvphy
integer :: id_tr1

logical :: initialized=.false.

namelist/phys_nml/dt_rad

contains

!--------------------------------------------------------------------------------   
subroutine init_phys(Time, dt_phys_in, dt_atmos_in, domain_in, ntrac_in, nlev_in, &
                     lat_deg_in, lon_deg_in, rstrt)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2d), target :: domain_in
    real, intent(in) :: dt_phys_in, dt_atmos_in 
    integer, intent(in) :: ntrac_in, nlev_in
    real, intent(in) :: lat_deg_in(:), lon_deg_in(:)
    type(restart_file_type), intent(inout) :: rstrt
    
    real, allocatable :: fland(:,:), ilevs(:), ilevsp(:)
    integer :: axes(4), id_lev, unit, i, j
    integer :: jsg,jeg,isg,ieg
    integer :: indx
    real :: tem1, tem2

    unit = open_namelist_file()

    read(unit,nml=phys_nml)

    call close_file(unit)

    dt_phys = dt_phys_in
    rdt_phys = 1./dt_phys

    dt_atmos = dt_atmos_in
    rdt_atmos = 1./dt_atmos

    if (mod(real(dt_rad),dt_atmos)/=0.) &
             call mpp_error('phys_mod', 'dt for radiation ' &
                  //'should be a multiple of model timestep', FATAL)

    time_step_rad = set_time(dt_rad)
    time_step = set_time(int(dt_rad))
    rad_time = Time

    domain => domain_in

    ntrac = ntrac_in
    
    nlev = nlev_in

    call mpp_get_compute_domain(domain,js,je,is,ie)
    call mpp_get_global_domain(domain,jsg,jeg,isg,ieg)

    jlen = je - js + 1
    ilen = ie - is + 1

    allocate(ilevs(nlev))
    allocate(ilevsp(nlev+1))

    forall(i=1:nlev) ilevs(i) = i
    forall(i=1:nlev+1) ilevsp(i) = i

    axes(1) = diag_axis_init('lat',lat_deg_in(js:je),'degrees_N','Y', & 
              long_name='latitude',domain_decomp=[jsg,jeg,js,je])
    axes(2) = diag_axis_init('lon',lon_deg_in(is:ie),'degrees_E','X', &
              long_name='longitude',domain_decomp=[isg,ieg,is,ie])
    axes(3) = diag_axis_init('lev',ilevs,'','Z',long_name='')
    axes(4) = diag_axis_init('levp',ilevsp,'','Z',long_name='')

    deallocate(ilevs,ilevsp)

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))
    allocate(gdlen(js:je,is:ie))

    do i = is, ie
        lat_deg(js:je,i) = lat_deg_in(js:je)
        lat_rad(js:je,i) = lat_deg_in(js:je)*PI/180.
    enddo

    do j = js, je
        lon_deg(j,is:ie) = lon_deg_in(is:ie)
        lon_rad(j,is:ie) = lon_deg_in(is:ie)*PI/180.
    enddo

    do i = is, ie
        do j = js, je
            tem1 = 2. * PI * RADIUS * cos(lat_rad(j,i)) / (ieg-isg+1)
            tem2 = PI * RADIUS / (jeg-jsg+1)
            gdlen(j,i) = sqrt(tem1**2 + tem2**2)
        enddo
    enddo

    allocate(htsw(nlev,js:je,is:ie))
    allocate(htlw(nlev,js:je,is:ie))
    allocate(rsdsz(js:je,is:ie))
    allocate(rsusz(js:je,is:ie))
    allocate(rldsz(js:je,is:ie))
    allocate(fprcp(js:je,is:ie))
    allocate(lprcp(js:je,is:ie))

    htsw  = 0. 
    htlw  = 0.
    rsdsz = 0.
    rsusz = 0.
    rldsz = 0.

    fprcp = 0.
    indx = reg_rf(rstrt, '', 'fprcp', fprcp, domain, mandatory=.false.)

    lprcp = 0.
    indx = reg_rf(rstrt, '', 'lprcp', lprcp, domain, mandatory=.false.)

    pdryini = 0.
    indx = reg_rf(rstrt, '', 'pdryini', pdryini, mandatory=.false.)

    call init_sfc(Time,dt_atmos,domain,axes(1:2),rstrt)
    allocate(fland(js:je,is:ie))
    allocate(slmsk(js:je,is:ie))
    call get_land_frac(fland)
    slmsk = 0.
    where(fland>=0.5) slmsk = 1.

    call astronomy_init()
    call init_radiation(Time, domain, ntrac, nlev, lat_deg_in, &
                        lon_deg_in, fland, axes, ind_q_in=ind_q, ind_clw_in=ind_clw) 

    call init_gwdrag(domain,fland)

    call init_cu_conv()

    call init_micro_phys(domain, nlev, rstrt, lat_deg_in*PI/180.) 

    deallocate(fland)

    ! Diag out
    ! -------------------------------------------------------------------------------- 
    
    id_rsds = reg_df(rou, 'rsds', axes(1:2), time, 'surface downwelling shortwave',  'w m-2')

    id_rsus = reg_df(rou, 'rsus', axes(1:2), time, 'surface upwelling shortwave',  'w m-2')

    id_rsns = reg_df(rou, 'rsns', axes(1:2), time, 'surface net shortwave',  'w m-2')

    id_sfcemis = reg_df(rou, 'sfcemis', axes(1:2), time, 'surface emisivity',  'w m-2')

    id_rlds = reg_df(rou, 'rlds', axes(1:2), time, 'surface downwelling longwave',  'w m-2')
    id_rlus = reg_df(rou, 'rlus', axes(1:2), time, 'surface upwelling longwave',  'w m-2')
    id_tskin = reg_df(rou, 'tskin', axes(1:2), time, 'Surface Temperature',  'K')

    id_dtrd = reg_df(rou, 'dtrd', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (radiation)',  'k s-1')

    id_dtlw = reg_df(rou, 'dtlw', [axes(3),axes(1),axes(2)], time, &
               'temperature tendency (lw)',  'ks-1')
    id_shflx = reg_df(rou, 'shflx', [axes(1),axes(2)], time, &
               'sensible heat flux', 'wm-2')
    id_lhflx = reg_df(rou, 'lhflx', [axes(1),axes(2)], time, &
               'latent heat flux', 'wm-2')
    id_taux = reg_df(rou, 'taux', [axes(1),axes(2)], time, &
               'u-wind stress', 'pa')
    id_tauy = reg_df(rou, 'tauy', [axes(1),axes(2)], time, &
               'v-wind stress', 'pa')
    id_hpbl = reg_df(rou, 'hpbl', [axes(1),axes(2)], time, &
               'height of planetary boundary layer', 'm')
    id_dtvd = reg_df(rou, 'dtvd', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (vertical diffusion)', 'k s-1')
    id_dqvd = reg_df(rou, 'dqvd', [axes(3),axes(1),axes(2)], time,  &
               'moisture tendency (vertical diffusion)', 's-1')
    id_duvd = reg_df(rou, 'duvd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (vertical diffusion)', 'm/s-2')
    id_dvvd = reg_df(rou, 'dvvd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (vertical diffusion)', 'm/s-2')
    id_dugwd = reg_df(rou, 'dugwd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (gravity wave drag)', 'm/s-2')
    id_dvgwd = reg_df(rou, 'dvgwd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (gravity wave drag)', 'm/s-2')

    id_ducgwd = reg_df(rou, 'ducgwd', [axes(3),axes(1),axes(2)], time,  &
               'u-velocity tendency (convective gravity wave drag)', 'm/s-2')
    id_dvcgwd = reg_df(rou, 'dvcgwd', [axes(3),axes(1),axes(2)], time,  &
               'v-velocity tendency (convective gravity wave drag)', 'm/s-2')

    id_dtcu = reg_df(rou, 'dtcu', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (cumulus conv)', 'k s-1')
    id_dqcu = reg_df(rou, 'dqcu', [axes(3),axes(1),axes(2)], time,  &
               'moisture tendency (cumulus conv)', 's-1')
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

    id_dtsc = reg_df(rou, 'dtsc', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (shallow conv)', 'k s-1')
    id_dqsc = reg_df(rou, 'dqsc', [axes(3),axes(1),axes(2)], time,  &
               'moisture tendency (shallow conv)', 's-1')

    id_dtmp = reg_df(rou, 'dtmp', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (micro phys)', 'k s-1')
    id_dqmp = reg_df(rou, 'dqmp', [axes(3),axes(1),axes(2)], time,  &
               'moisture tendency (micro phys)', 's-1')
    id_lprmp = reg_df(rou, 'lprmp', [axes(1),axes(2)], time, &
                'rainfall rate (micro phys)', 'kg m-2 s-1')
    id_fprmp = reg_df(rou, 'fprmp', [axes(1),axes(2)], time, &
                'snowfall rate (micro phys)', 'kg m-2 s-1')

    id_dtphy = reg_df(rou, 'dtphy', [axes(3),axes(1),axes(2)], time,  &
               'temperature tendency (phys)', 'k s-1')
    id_dqphy = reg_df(rou, 'dqphy', [axes(3),axes(1),axes(2)], time,  &
               'moisture tendency (phys)', 's-1')
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

    id_tr1 = reg_df(rou, 'tr1', [axes(3),axes(1),axes(2)], time, &
               'temperature tendency (lw)',  'ks-1')
    !--------------------------------------------------------------------------------    

    initialized = .true.

end subroutine init_phys


!--------------------------------------------------------------------------------   
subroutine phys(Time,tlyr,tr,p,tlyr1,tr1,p1,u1,v1,vvel1,dtdt,dqdt,dudt,dvdt,topo)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(1:nlev,js:je,is:ie) :: tlyr, tlyr1, u1, v1, vvel1
    real, intent(in), dimension(1:nlev,js:je,is:ie,1:ntrac) :: tr, tr1
    real, intent(in), dimension(js:je,is:ie) :: p, p1
    real, intent(out), dimension(1:nlev,js:je,is:ie) :: dtdt, dudt, dvdt 
    real, intent(out), dimension(1:nlev,js:je,is:ie,ntrac) :: dqdt
    real, intent(in), optional, dimension(js:je,is:ie) :: topo

    real, dimension(1:nlev,js:je,is:ie) :: plyr, plyrk, prslki, delp, phil, tmp3d
    real, dimension(1:nlev+1,js:je,is:ie) :: plvl, plvlk, phii
    real, dimension(1:nlev,js:je,is:ie) :: dtdt1, dudt1, dvdt1, clw, cli
    real, dimension(1:nlev,js:je,is:ie,ntrac) :: dqdt1

    real, dimension(js:je,is:ie) :: tskin, fracday, coszen, rcoszen, rsds, rsus, rsns, &
                                    rlds, rlus, rldsg, rb, ffmm, ffhh, qss, hflx, evap, &
                                    stress, wind, dusfc1, dvsfc1, dtsfc1, dqsfc1, hpbl, &
                                    gamt, gamq, xkzm, topo1, sfcemis, cldwrk, rain1, rain, &
                                    snow, snow1

    integer, dimension(js:je,is:ie) :: kpbl, kbot, ktop, kcnv

    real, dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real :: solcon, rrsun
    integer :: k, imax
    logical :: used

    imax = size(p,1)*size(p,2)

    dtdt = 0.; dudt = 0.; dvdt = 0.; dqdt = 0.
    dtdt1 = 0.; dudt1 = 0.; dvdt1 = 0.; dqdt1 = 0.
    rain = 0.; snow = 0.; rain1 = 0.; snow1 = 0.

    topo1 = 0.
    if(present(topo)) topo1 = topo

    !Radiation
    !--------------------------------------------------------------------------------   
    if (rad_time==Time) then
        call print_date(Time,'Calling Radiation at:')
        coszen(:,:) = 0.0
        rcoszen(:,:) = 0.0
        call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step_rad)
        call set_surface(Time,tskin,coszen,sfcalb,sfcemis)
        solcon = con_solr * rrsun
        where(coszen>0.) rcoszen = 1./coszen 
        call radiation(Time, tlyr, tr, p, tskin, coszen, fracday, sfcalb, sfcemis, &
                       solcon, htsw, rsdsz, rsusz, htlw, rldsz, rlus)
        do k = 1, size(htsw,1)
            htsw(k,:,:) = htsw(k,:,:) * rcoszen
        enddo
        rsdsz = rsdsz * rcoszen
        rsusz = rsusz * rcoszen
        rldsz = rldsz / (tlyr(1,:,:)**4.)
        rad_time = rad_time + time_step_rad
    else
        call set_surface(Time,tskin,sfcemis=sfcemis)
    endif

    call diurnal_solar(lat_rad, lon_rad, Time, coszen, fracday, rrsun, time_step)
    coszen = coszen * fracday

    call adjust_rad(coszen, tskin, tlyr(1,:,:), rsdsz, rsusz, rldsz, htsw,  &
                    htlw, sfcemis, dtdt1, rsds, rsus, rlds, rlus)
    dtdt = dtdt + dtdt1

    rsns = rsds - rsus

    rldsg = rlds * sfcemis !-> rldsg for land,ocn,sea-ice

    used = send_data(id_rsds, rsds, Time) 
    used = send_data(id_rsus, rsus, Time) 
    used = send_data(id_rsns, rsns, Time) 
    used = send_data(id_rlds, rlds, Time) 
    used = send_data(id_rlus, rlus, Time) 
    used = send_data(id_tskin, tskin, Time) 
    used = send_data(id_sfcemis, sfcemis, Time) 
    used = send_data(id_dtrd, dtdt1, Time) 
    used = send_data(id_dtlw, htlw, Time) 
    !End Radiation
    !--------------------------------------------------------------------------------   

    call get_pressure_at_levels(p1,plvl,plyr,plvlk,plyrk)
    prslki = plvlk(1:nlev,:,:)/plyrk
    delp = plvl(1:nlev,:,:) - plvl(2:nlev+1,:,:)

    call get_phi(topo1, tlyr1, tr1(:,:,:,ind_q), plvl, plvlk, plyr, plyrk, phii, phil)

    !Surface Fluxes
    call do_surface(Time, p1(:,:), u1(1,:,:), v1(1,:,:), tlyr1(1,:,:), &
            tr1(1,:,:,ind_q), plyr(1,:,:), prslki(1,:,:), rsds, rsns, rldsg, &
            fprcp, lprcp, rb, ffmm, ffhh, qss, hflx, evap, stress, wind)
    !End surface Fluxes

    !Vertical Diffusion
    dtdt1 = 0.; dudt1 = 0.; dvdt1 = 0.; dqdt1 = 0.
    dusfc1 = 0.; dvsfc1 = 0.; dtsfc1 = 0.; dqsfc1 = 0.
    call do_vertical_diffusion(imax, nlev, ntrac, dvdt1, dudt1, dtdt1, dqdt1, u1, &
                    v1, tlyr1, tr1, plvlk, rb, ffmm, ffhh, qss, hflx, evap, stress, &
                    wind, kpbl, plvl, delp, plyr, plyrk, phii, phil, dt_phys, dusfc1, &
                    dvsfc1, dtsfc1, dqsfc1, hpbl, gamt, gamq, xkzm)

    used = send_data(id_shflx, dtsfc1, Time)
    used = send_data(id_lhflx, dqsfc1, Time)
    used = send_data(id_hpbl, hpbl, Time)
    used = send_data(id_taux, dusfc1, Time)
    used = send_data(id_tauy, dvsfc1, Time)
    used = send_data(id_dtvd, dtdt1, Time)
    used = send_data(id_dqvd, dqdt1(:,:,:,ind_q), Time)
    used = send_data(id_duvd, dudt1, Time)
    used = send_data(id_dvvd, dvdt1, Time)


    dtdt = dtdt + dtdt1 
    dudt = dudt + dudt1 
    dvdt = dvdt + dvdt1 
    dqdt = dqdt + dqdt1 
    !End Vertical Diffusion

    !Gravity Wave Drag
    dudt1 = 0.; dvdt1 = 0.
    call gwdrag(dvdt1, dudt1, u1, v1, tlyr1, tr1(:,:,:,ind_q), kpbl, &
                plvl, delp, plyr, plyrk, phii, phil, dt_phys, dusfc1, dvsfc1)

    used = send_data(id_dugwd, dudt1, Time)
    used = send_data(id_dvgwd, dvdt1, Time)

    dudt = dudt + dudt1 
    dvdt = dvdt + dvdt1 
    !End Gravity Wave Drag

    dtdt = tlyr1 + dtdt*dt_phys 
    dudt = u1 + dudt*dt_phys 
    dvdt = v1 + dvdt*dt_phys 
    dqdt = tr1 + dqdt*dt_phys 
    
    call get_phi(topo1, dtdt, dqdt(:,:,:,ind_q), plvl, plvlk, plyr, plyrk, phii, phil)

    !Cumulus convection
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
    if (id_dqcu>0) used = send_data(id_dqcu, rdt_phys*(dqdt1(:,:,:,ind_q) &
                                            - dqdt(:,:,:,ind_q)), Time)
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
    !End Cumulus Convection


    ! Convective Gravity Wave Drag
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
    !End Convective Gravity Wave Drag


    ! Shallow Convection
    dtdt1 = dtdt
    dqdt1 = dqdt !Initialize
    call shallow_conv(imax, nlev, dt_phys, delp, plvl, plyr, plyrk, kcnv, &
                      dqdt1(:,:,:,ind_q), dtdt1)
    if (id_dtsc>0) used = send_data(id_dtsc, rdt_phys*(dtdt1-dtdt), Time)
    if (id_dqsc>0) used = send_data(id_dqsc, rdt_phys*(dqdt1(:,:,:,ind_q) &
                                            - dqdt(:,:,:,ind_q)), Time)
    dtdt = dtdt1
    dqdt = dqdt1
    !End Shallow Convection

    !Micro Physics
    dtdt1 = dtdt
    dqdt1 = dqdt 

    call micro_phys(dt_phys, plyr, p1, plyrk, delp, dqdt1(:,:,:,ind_q), &
                    dqdt1(:,:,:,ind_clw), dtdt1, rain1, snow1)

    if (id_dtmp>0) used = send_data(id_dtmp, rdt_phys*(dtdt1-dtdt), Time)
    if (id_dqmp>0) used = send_data(id_dqmp, rdt_phys*(dqdt1(:,:,:,ind_q) &
                                            - dqdt(:,:,:,ind_q)), Time)
    if (id_lprmp>0) used = send_data(id_lprmp, rdt_phys*rain1*rhoh2o, Time)
    if (id_fprmp>0) used = send_data(id_fprmp, rdt_phys*snow1*rhoh2o, Time)

    dtdt = dtdt1
    dqdt = dqdt1
    rain = rain + rain1
    snow = snow + snow1
    !End Micro Physics

    lprcp = rain * rdt_phys
    fprcp = snow * rdt_phys

    if (id_lpr>0) used = send_data(id_lpr, rdt_phys*rain*rhoh2o, Time)
    if (id_fpr>0) used = send_data(id_fpr, rdt_phys*snow*rhoh2o, Time)
    if (id_pr>0) used = send_data(id_pr, rdt_phys*(rain+snow)*rhoh2o, Time)

    if (id_dtphy>0) used = send_data(id_dtphy, rdt_phys*(dtdt-tlyr1), Time)
    if (id_dqphy>0) used = send_data(id_dqphy, rdt_phys*(dqdt(:,:,:,ind_q) &
                                            - tr1(:,:,:,ind_q)), Time)
    if (id_duphy>0) used = send_data(id_duphy, rdt_phys*(dudt-u1), Time)
    if (id_dvphy>0) used = send_data(id_dvphy, rdt_phys*(dvdt-v1), Time)

    if (id_tr1>0) used = send_data(id_tr1, tr1(:,:,:,1), Time)

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

end module phys_mod

