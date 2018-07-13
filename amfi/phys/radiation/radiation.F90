module radiation_mod

use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain

use fms_mod, only : open_namelist_file, close_file, read_data, write_data

use time_manager_mod, only : time_type, set_time

use data_override_mod, only : data_override

use swrad_mod, only : init_swrad, swrad, NBDSW

use vertical_levels_mod, only : get_pressure_at_levels

use astronomy_mod, only : astronomy_init, diurnal_solar

use sfc_mod, only : get_albedo

implicit none
private

public :: init_radiation

real :: deltim=600.
real :: deltimr=3600.

type(time_type) :: dt, rad_dt

integer :: icwp=1, iovr=1, isubc=2

type(domain2D), pointer :: domain

integer :: is, ie, js, je, nlev, ntrac, nlevp1
integer :: ilen, jlen, imax

integer :: ngases=9, ncldparms=9, nsfcaldparms=4

real, parameter :: con_solr = 1.3610e+3      ! solar constant (W/m2)-liu(2002)

integer, parameter :: NF_AESW=3, NF_AELW=3
integer, parameter :: NF_CLDS=9, NF_VGAS=9
integer, parameter :: NF_ALBD=4

logical :: initialized=.true.

contains

!--------------------------------------------------------------------------------   
subroutine init_radiation(Time,deltim_in,domain_in,ntrac_in,nlev_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in) :: deltim_in
    type(domain2D), target :: domain_in
    integer, intent(in) :: nlev_in, ntrac_in
    integer :: unit
  
    namelist/radiation_nml/icwp, iovr, isubc, deltimr

    unit = open_namelist_file()

    read(unit,nml=radiation_nml)

    call close_file(unit)

    deltim = deltim_in

    dt = set_time(seconds=int(deltim))
    rad_dt = set_time(seconds=int(deltimr))

    domain => domain_in

    nlev = nlev_in
    nlevp1 = nlev_in + 1
    
    ntrac = ntrac_in

    if (mod(deltimr,deltim)/=0.) call mpp_error('radiation_mod', 'deltim for radiation ' &
                                        //'should be a multiple of model timestep', FATAL)

    call mpp_get_compute_domain(domain,js,je,is,ie)

    jlen = je - js + 1
    ilen = ie - is + 1
    imax = jlen*ilen

    call astronomy_init()

    call init_swrad(icwp,iovr,isubc)

    initialized = .true.

end subroutine init_radiation

!--------------------------------------------------------------------------------   
subroutine radiation(Time, tlyr, tr, p, tslnd, tsocn, tsice, lfrac, &
                     ifrac, ofrac, snowf, zorl, hice, hsnow, xlat, xlon)
!--------------------------------------------------------------------------------   

    type(time_type), intent(in) :: Time
    real, intent(in) :: tlyr(1:nlev,js:je,is:ie)
    real, intent(in) :: tr(1:nlev,js:je,is:ie,1:ntrac)
    real, intent(in), dimension(js:je,is:ie) :: p   !-> is in cb, should be converted to mb
    real, intent(in), dimension(js:je,is:ie) :: tslnd, tsocn, tsice
    real, intent(in), dimension(js:je,is:ie) :: lfrac, ifrac, ofrac
    real, intent(in), dimension(js:je,is:ie) :: snowf, zorl, hice, hsnow
    real, intent(in), dimension(js:je,is:ie) :: xlat, xlon

    real, dimension(1:nlev,js:je,is:ie) :: plyr, qlyr, olyr
    real, dimension(1:nlevp1,js:je,is:ie) :: plvl, tlvl

    real, dimension(1:nlev,js:je,is:ie,NF_CLDS) :: clouds
    real, dimension(1:nlev,js:je,is:ie,NF_VGAS) :: gasvmr
    real, dimension(js:je,is:ie,NF_ALBD)        :: sfcalb
    real, dimension(js:je,is:ie,NF_ALBD)        :: sfcalblnd
    real, dimension(js:je,is:ie,NF_ALBD)        :: sfcalbocn
    real, dimension(js:je,is:ie,NF_ALBD)        :: sfcalbice
    real, dimension(1:nlev,js:je,is:ie,NBDSW,NF_AESW) :: faersw

    real, dimension(1:nlev,js:je,is:ie) :: hswc
    real, dimension(1:nlev,js:je,is:ie) :: hlwc

    real, dimension(js:je,is:ie) :: rsdt, rsut, rsutc
    real, dimension(js:je,is:ie) :: rsds, rsdsc, rsus, rsusc
    real, dimension(js:je,is:ie) :: ruvds, ruvdsc
    real, dimension(js:je,is:ie) :: rnirdsbm, rnirdsdf
    real, dimension(js:je,is:ie) :: rvisdsbm, rvisdsdf

    real, dimension(js:je,is:ie) :: coszen, coszdg
    integer, dimension(js:je,is:ie) :: icseed

    real, dimension(js:je,is:ie) :: tskin, alvsf, alnsf, alvwf, alnwf
    real, dimension(js:je,is:ie) :: facsf, facwf

    real :: tem2db(nlevp1,js:je,is:ie), tem2da(nlev,js:je,is:ie)

    logical, dimension(js:je,is:ie) :: lmask, imask, omask

    real :: solcon, rrsun

    integer :: k

    if(.not.initialized) call mpp_error(FATAL,'radiation_mod: module not initialized')
    
    call get_pressure_at_levels(p,plvl,plyr)

    lmask = .false.; imask = .false.; omask = .false.

    where(lfrac>0.) lmask = .true.
    where(ofrac>0.) omask = .true.
    where(ifrac>0.) imask = .true.

    plvl = plvl * 10.0 !cb to mb
    plyr = plyr * 10.0 !cb to mb

    tlvl(nlevp1,:,:) = tlyr(nlev,:,:)

    tskin = tslnd*lfrac + ofrac*tsocn + ifrac*tsice

    tlvl(1,:,:) = tskin(:,:)

    tem2db = log(plvl)
    tem2da = log(plyr)

    do k = 1, nlev-1
        tlvl(k+1,:,:) = tlyr(k,:,:) + (tlyr(k+1,:,:) - tlyr(k,:,:)) &
                                    * (tem2db(k+1,:,:) - tem2da(k,:,:)) &
                                    / (tem2da(k+1,:,:) - tem2da(k,:,:))
    enddo

    qlyr = tr(:,:,:,1)

    call data_override('ATM','alvsf',alvsf,Time) 
    call data_override('ATM','alnsf',alnsf,Time) 
    call data_override('ATM','alvwf',alvwf,Time) 
    call data_override('ATM','alnwf',alnwf,Time) 
    call data_override('ATM','facsf',facsf,Time) 
    call data_override('ATM','facwf',facwf,Time) 

    coszen(:,:) = 0.0
    coszdg(:,:) = 0.0

    call diurnal_solar(xlat, xlon, Time, coszen, coszdg, rrsun, rad_dt)
    solcon = con_solr * rrsun
    coszdg = coszen * coszdg

    call get_albedo(Time,sfcalb)

    faersw = 0.

    call swrad_drv(plyr,plvl,tlyr,tlvl,qlyr, &
                   olyr,gasvmr,clouds,icseed,faersw,sfcalb,coszen,solcon, &
                   hswc,rsdt,rsut,rsutc,rsds,rsdsc,rsus,rsusc,ruvds, &
                   ruvdsc,rnirdsbm,rnirdsdf,rvisdsbm,rvisdsdf)

end subroutine radiation



!--------------------------------------------------------------------------------   
subroutine swrad_drv(plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr, &
                     clouds,icseed,aerosols,sfcalb,cosz,solcon, &
                     hswc,topdnfxc,topupfxc,topupfx0, &
                     sfcdnfxc,sfcdnfx0,sfcupfxc,sfcupfx0, &
                     sfcuvbfc,sfcuvbf0,sfcnirbm,sfcnirdf, &
                     sfcvisbm,sfcvisdf)
!--------------------------------------------------------------------------------   

    real, dimension(nlev,imax), intent(in) :: plyr, tlyr, qlyr, olyr 
    real, dimension(nlevp1,imax), intent(in) :: plvl, tlvl
    real, intent(in) :: gasvmr(nlev,imax,NF_VGAS), clouds(nlev,imax,NF_CLDS)
    integer, intent(in) :: icseed(imax)
    real, intent(in) :: aerosols(nlev,imax,NBDSW,NF_AESW)
    real, intent(in) :: sfcalb(imax,NF_ALBD), cosz(imax), solcon
    real, intent(out) :: hswc(nlev,imax)
    real, dimension(imax), intent(out) :: topdnfxc,topupfxc,topupfx0
    real, dimension(imax), intent(out) :: sfcdnfxc,sfcdnfx0,sfcupfxc,sfcupfx0
    real, dimension(imax), intent(out) :: sfcuvbfc,sfcuvbf0,sfcnirbm,sfcnirdf
    real, dimension(imax), intent(out) :: sfcvisbm,sfcvisdf

    integer :: i


    do i = 1, imax
        call swrad(nlev, nlevp1, plyr(:,i), plvl(:,i), tlyr(:,i), tlvl(:,i), &
                qlyr(:,i), olyr(:,i), gasvmr(:,i,:), clouds(:,i,:), icseed(i), &
                aerosols(:,i,:,:), sfcalb(i,:), cosz(i), solcon, &
                hswc(:,i), topdnfxc(i), topupfxc(i), topupfx0(i), &
                sfcdnfxc(i), sfcdnfx0(i), sfcupfxc(i), sfcupfx0(i), &
                sfcuvbfc(i), sfcuvbf0(i), sfcnirbm(i), sfcnirdf(i), &
                sfcvisbm(i), sfcvisdf(i))

    enddo


end subroutine swrad_drv


end module radiation_mod
