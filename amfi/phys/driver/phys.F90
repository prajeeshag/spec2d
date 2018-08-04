module phys_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain
use mpp_domains_mod, only : mpp_get_global_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use fms_io_mod, only : restart_file_type, reg_rf => register_restart_field
use fms_io_mod, only : restore_state, save_restart

use constants_mod, only : PI, CP_AIR, RDGAS, RVGAS

use time_manager_mod, only : time_type, set_time, operator(==), operator(+), assignment(=)

use diag_manager_mod, only : diag_axis_init, reg_df=>register_diag_field, send_data

use radiation_mod, only : init_radiation, con_solr, radiation, NF_ALBD

use sfc_mod, only : init_sfc, get_land_frac, set_surface, do_surface

use vertdiff_mod, only : do_vertical_diffusion

use astronomy_mod, only : astronomy_init, diurnal_solar

use vertical_levels_mod, only : get_pressure_at_levels

implicit none
private

public :: init_phys, phys

type(domain2D), pointer :: domain
integer :: ntrac, nlev
integer :: is, ie, ilen, js, je, jlen
integer :: nlon, nlat
integer :: ind_q = 1, ind_clw = 2

integer :: dt_rad
real :: deltim, deltimr
type(time_type) :: time_step_rad, rad_time, time_step

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg
real, allocatable, dimension(:,:,:) :: htsw
real, allocatable, dimension(:,:,:) :: htlw
real, allocatable, dimension(:,:) :: rldsz
real, allocatable, dimension(:,:) :: rsdsz, rsusz
real, allocatable, dimension(:,:) :: fprcp, lprcp

character(len=16) :: resfnm = 'phys_res'
character (len=8) :: rou='am_phys'

integer :: id_rsds, id_rsus, id_rsns, id_htlw, id_htrd, id_shflx, id_lhflx, id_taux, &
           id_tauy, id_htvd, id_hpbl

logical :: initialized=.false.

namelist/phys_nml/dt_rad

contains

!--------------------------------------------------------------------------------   
subroutine init_phys(Time, deltim_in, domain_in, ntrac_in, nlev_in, &
                     lat_deg_in, lon_deg_in, rstrt)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2d), target :: domain_in
    real, intent(in) :: deltim_in 
    integer, intent(in) :: ntrac_in, nlev_in
    real, intent(in) :: lat_deg_in(:), lon_deg_in(:)
    type(restart_file_type), intent(inout) :: rstrt
    
    real, allocatable :: fland(:,:), ilevs(:), ilevsp(:)
    integer :: axes(4), id_lev, unit, i, j
    integer :: jsg,jeg,isg,ieg
    integer :: indx

    unit = open_namelist_file()

    read(unit,nml=phys_nml)

    call close_file(unit)

    deltim = deltim_in

    deltimr = real(dt_rad)

    if (mod(deltimr,deltim)/=0.) call mpp_error('phys_mod', 'deltim for radiation ' &
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

    call init_sfc(Time,deltim,domain,axes(1:2),rstrt)

    allocate(fland(js:je,is:ie))

    call get_land_frac(fland)

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))

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

    call restore_state(rstrt)

    do i = is, ie
        lat_deg(js:je,i) = lat_deg_in(js:je)
        lat_rad(js:je,i) = lat_deg_in(js:je)*PI/180.
    enddo

    do j = js, je
        lon_deg(j,is:ie) = lon_deg_in(is:ie)
        lon_rad(j,is:ie) = lon_deg_in(is:ie)*PI/180.
    enddo

    call astronomy_init()
    call init_radiation(Time, deltim, domain, ntrac, nlev, lat_deg_in, &
                        lon_deg_in, fland, axes, ind_q_in=ind_q, ind_clw_in=ind_clw) 


    deallocate(fland)

    ! Diag out
    ! -------------------------------------------------------------------------------- 
    
    id_rsds = reg_df(rou, 'rsds', axes(1:2), Time, 'Surface Downwelling Shortwave',  'W m-2')

    id_rsus = reg_df(rou, 'rsus', axes(1:2), Time, 'Surface Upwelling Shortwave',  'W m-2')

    id_rsns = reg_df(rou, 'rsns', axes(1:2), Time, 'Surface Net Shortwave',  'W m-2')

    id_htrd = reg_df(rou, 'htrd', [axes(3),axes(1),axes(2)], Time,  &
               'Heating Rate Due to Radiation',  'K s-1')

    id_htlw = reg_df(rou, 'htlw', [axes(3),axes(1),axes(2)], Time, &
               'Heating Rate Due to LW',  'Ks-1')
    id_shflx = reg_df(rou, 'shflx', [axes(1),axes(2)], Time, &
               'Sensible Heat Flux', 'Wm-2')
    id_lhflx = reg_df(rou, 'lhflx', [axes(1),axes(2)], Time, &
               'Latent Heat Flux', 'Wm-2')
    id_taux = reg_df(rou, 'taux', [axes(1),axes(2)], Time, &
               'U-Wind Stress', 'Pa')
    id_tauy = reg_df(rou, 'tauy', [axes(1),axes(2)], Time, &
               'V-Wind Stress', 'Pa')
    id_hpbl = reg_df(rou, 'hpbl', [axes(1),axes(2)], Time, &
               'Height of Planetary Boundary Layer', 'm')
    id_htvd = reg_df(rou, 'htvd', [axes(3),axes(1),axes(2)], Time,  &
               'Heating Rate Due to Vertical Diffusion', 'K s-1')
    !--------------------------------------------------------------------------------    

    initialized = .true.

end subroutine init_phys


!--------------------------------------------------------------------------------   
subroutine phys(Time,tlyr,tr,p,tlyr1,tr1,p1,u1,v1,topo)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(1:nlev,js:je,is:ie) :: tlyr, tlyr1, u1, v1
    real, intent(in), dimension(1:nlev,js:je,is:ie,1:ntrac) :: tr, tr1
    real, intent(in), dimension(js:je,is:ie) :: p, p1
    real, intent(in), optional, dimension(js:je,is:ie) :: topo

    real, dimension(1:nlev,js:je,is:ie) :: plyr, plyrk, prslki, delp, phil
    real, dimension(1:nlev+1,js:je,is:ie) :: plvl, plvlk, phii
    real, dimension(1:nlev,js:je,is:ie) :: dtdt, dudt, dvdt
    real, dimension(1:nlev,js:je,is:ie) :: dtdt1, dudt1, dvdt1
    real, dimension(1:nlev,js:je,is:ie,ntrac) :: dqdt, dqdt1

    real, dimension(js:je,is:ie) :: tskin, fracday, coszen, rcoszen, rsds, rsus, rsns, &
                                    rlds, rlus, rldsg, rb, ffmm, ffhh, qss, hflx, evap, &
                                    stress, wind, dusfc1, dvsfc1, dtsfc1, dqsfc1, hpbl, &
                                    gamt, gamq, xkzm, topo1, sfcemis

    integer, dimension(js:je,is:ie) :: kpbl

    real, dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real :: solcon, rrsun
    integer :: k, imax
    logical :: used

    imax = size(p,1)*size(p,2)

    dtdt = 0.; dudt = 0.; dvdt = 0.; dqdt = 0.
    dtdt1 = 0.; dudt1 = 0.; dvdt1 = 0.; dqdt1 = 0.

    topo1 = 0.
    if(present(topo)) topo1 = topo

    if (rad_time==Time) then
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
    used = send_data(id_htrd, dtdt1, Time) 
    used = send_data(id_htlw, htlw, Time) 

    call get_pressure_at_levels(p1,plvl,plyr,plvlk,plyrk)
    prslki = plvlk(1:nlev,:,:)/plyrk
    delp = plvl(1:nlev-1,:,:) - plvl(2:nlev,:,:)

    call get_phi(topo1, tlyr1, tr1(:,:,:,ind_q), plvl, plvlk, plyr, plyrk, phii, phil)
    
    call do_surface(Time, p1(:,:), u1(1,:,:), v1(1,:,:), tlyr1(1,:,:), &
            tr1(1,:,:,ind_q), plyr(1,:,:), prslki(1,:,:), rsds, rsns, rldsg, &
            fprcp, lprcp, rb, ffmm, ffhh, qss, hflx, evap, stress, wind)

    call do_vertical_diffusion(imax, nlev, ntrac, dvdt1, dudt1, dtdt1, dqdt1, u1, &
                    v1, tlyr1, tr1, plvlk, rb, ffmm, ffhh, qss, hflx, evap, stress, &
                    wind, kpbl, plvl, delp, plyr, plyrk, phii, phil, deltim, dusfc1, &
                    dvsfc1, dtsfc1, dqsfc1, hpbl, gamt, gamq, xkzm)

    dtdt = dtdt + dtdt1 
    dudt = dudt + dudt1 
    dvdt = dvdt + dvdt1 
    dqdt = dqdt + dqdt1 

    used = send_data(id_shflx, dtsfc1, Time)
    used = send_data(id_lhflx, dqsfc1, Time)
    used = send_data(id_hpbl, hpbl, Time)
    used = send_data(id_taux, dusfc1, Time)
    used = send_data(id_tauy, dvsfc1, Time)
    used = send_data(id_htvd, dtdt1, Time)

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

