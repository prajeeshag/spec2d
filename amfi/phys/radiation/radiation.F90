module radiation_mod

use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain

use fms_mod, only : open_namelist_file, close_file, read_data, write_data

use time_manager_mod, only : time_type, set_time, get_date

use data_override_mod, only : data_override

use diag_manager_mod, only : diag_axis_init, register_diag_field, send_data

use swrad_mod, only : init_swrad, swrad, NBDSW

use vertical_levels_mod, only : get_pressure_at_levels

use sat_vapor_pres_mod, only : compute_qs

use constants_mod, only : PI

use mersenne_twister, only : random_setseed, random_index, random_stat

use module_radiation_clouds, only : progcld1

implicit none
private

public :: init_radiation, radiation

real :: deltim=0.

type(time_type) :: dt, rad_dt

type(domain2D), pointer :: domain

integer :: is, ie, js, je, nlev, ntrac, nlevp1
integer :: ilen, jlen, imax

real, public, parameter :: con_solr = 1.3610e+3      ! solar constant (W/m2)-liu(2002)

integer, public, parameter :: NF_AESW=3, NF_AELW=3
integer, public, parameter :: NF_CLDS=9, NF_VGAS=10
integer, public, parameter :: NF_ALBD=4
integer :: icwp=1, iovr=1, isubc=2

integer :: ind_q=-1, ind_clw=-1, ind_oz=-1

real :: co2vmr = 284.7E-6
real :: n2ovmr =  0.275E-6
real :: ch4vmr =  0.790E-6
real :: o2vmr  =  0.209
real :: covmr  = 1.50E-8
real :: f11vmr = 0.0
real :: f12vmr = 0.0
real :: f22vmr = 0.0
real :: cl4vmr = 0.0
real :: f113vmr = 0.0

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg
real, allocatable, dimension(:,:) :: fland

integer :: id_plvl, id_plyr, id_tskin, id_albedo, id_tlyr, id_tlvl, id_coszen, &
           id_coszdg, id_ozone

logical :: initialized=.true.

namelist/radiation_nml/icwp, iovr, isubc, co2vmr, n2ovmr, ch4vmr, &
                       o2vmr, covmr, f11vmr, f12vmr, f22vmr, cl4vmr, f113vmr

contains

!--------------------------------------------------------------------------------   
subroutine init_radiation(Time,deltim_in,domain_in,ntrac_in,nlev_in, &
                          lat_deg_in,lon_deg_in,fland_in, axes_in, &
                          ind_q_in,ind_clw_in,ind_oz_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in) :: deltim_in
    type(domain2D), target :: domain_in
    integer, intent(in) :: nlev_in, ntrac_in
    real, intent(in) :: lat_deg_in(:), lon_deg_in(:)
    real, intent(in) :: fland_in(:,:)
    integer, intent(in) :: axes_in(:)
    integer, intent(in), optional :: ind_q_in, ind_oz_in, ind_clw_in
    integer :: unit, i, j
  
    unit = open_namelist_file()

    read(unit,nml=radiation_nml)

    call close_file(unit)

    deltim = deltim_in

    domain => domain_in

    nlev = nlev_in
    nlevp1 = nlev_in + 1
    
    ntrac = ntrac_in

    if(present(ind_q_in)) ind_q = ind_q_in
    if(present(ind_oz_in)) ind_oz = ind_oz_in
    if(present(ind_clw_in)) ind_clw = ind_clw_in

    call mpp_get_compute_domain(domain,js,je,is,ie)

    jlen = je - js + 1
    ilen = ie - is + 1
    imax = jlen*ilen

    allocate(lat_rad(js:je,is:ie))
    allocate(lat_deg(js:je,is:ie))
    allocate(lon_rad(js:je,is:ie))
    allocate(lon_deg(js:je,is:ie))
    allocate(fland(js:je,is:ie))

    fland(js:je,is:ie) = fland_in

    do i = is, ie
        lat_deg(js:je,i) = lat_deg_in(js:je)
        lat_rad(js:je,i) = lat_deg_in(js:je)*PI/180.
    enddo

    do j = js, je
        lon_deg(j,is:ie) = lon_deg_in(is:ie)
        lon_rad(j,is:ie) = lon_deg_in(is:ie)*PI/180.
    enddo

    call init_swrad(icwp,iovr,isubc)

    call rad_diag_init(Time,axes_in)

    initialized = .true.

end subroutine init_radiation


!--------------------------------------------------------------------------------   
subroutine rad_diag_init(Time,axes_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer, intent(in) :: axes_in(:)

    integer :: id_lev, id_levp, id_lat, id_lon
    integer :: axes(3), axesp(3), i
    real, allocatable :: ilevs(:), ilevsp(:)

    allocate(ilevs(nlev)) 
    allocate(ilevsp(nlev+1)) 

    forall(i=1:nlev) ilevs(i) = i
    forall(i=1:nlev+1) ilevsp(i) = i

    id_lat = axes_in(1) 
    id_lon = axes_in(2) 
    id_lev = diag_axis_init('lev',ilevs,'','Z',long_name='')
    id_levp = diag_axis_init('levp',ilevsp,'','Z',long_name='')

    axes = [id_lev,id_lat,id_lon] 
    axesp = [id_levp,id_lat,id_lon] 

    id_plvl = register_diag_field('amfi_rad','plvl',axesp,Time,'Pressure','mb')
    id_plyr = register_diag_field('amfi_rad','plyr',axes,Time,'Pressure','mb')
    id_tskin = register_diag_field('amfi_rad','tskin',axes(2:3),Time,'Surface Temperature','K')
    id_albedo = register_diag_field('amfi_rad','albedo',axes(2:3),Time,'Surface Albedo 1','1')
    id_tlyr = register_diag_field('amfi_rad','tlyr',axes,Time,'Temperature','K')
    id_tlvl = register_diag_field('amfi_rad','tlvl',axesp,Time,'Temperature','K')
    id_coszen = register_diag_field('amfi_rad','coszen',axes(2:3),Time,'coszen','1')
    id_coszdg = register_diag_field('amfi_rad','coszdg',axes(2:3),Time,'coszdg','1')
    id_ozone = register_diag_field('amfi_rad','ozone',axes,Time,'Ozone','1')

end subroutine rad_diag_init

!--------------------------------------------------------------------------------   
subroutine radiation(Time, tlyr, tr, p, tskin, coszen, coszdg, sfcalb, solcon)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in) :: tlyr(1:nlev,js:je,is:ie)
    real, intent(in) :: tr(1:nlev,js:je,is:ie,1:ntrac)
    real, intent(in), dimension(js:je,is:ie) :: p   !-> is in cb, should be converted to mb
    real, intent(in), dimension(js:je,is:ie) :: tskin, coszen, coszdg
    real, intent(in), dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real, intent(in) :: solcon

    real, dimension(1:nlev,js:je,is:ie) :: plyr, qlyr, olyr, clw
    real, dimension(1:nlevp1,js:je,is:ie) :: plvl, tlvl

    real, dimension(1:nlev,NF_CLDS,js:je,is:ie) :: clouds
    real, dimension(1:nlev,js:je,is:ie,NF_VGAS) :: gasvmr
    real, dimension(1:nlev,NBDSW,NF_AESW,js:je,is:ie) :: faersw

    real, dimension(1:nlev,js:je,is:ie) :: hswc
    real, dimension(1:nlev,js:je,is:ie) :: hlwc

    real, dimension(js:je,is:ie) :: rsdt, rsut, rsutc
    real, dimension(js:je,is:ie) :: rsds, rsdsc, rsus, rsusc
    real, dimension(js:je,is:ie) :: ruvds, ruvdsc
    real, dimension(js:je,is:ie) :: rnirdsbm, rnirdsdf
    real, dimension(js:je,is:ie) :: rvisdsbm, rvisdsdf

    integer, dimension(js:je,is:ie,2) :: icseed

    real :: tem2db(nlevp1,js:je,is:ie), tem2da(nlev,js:je,is:ie)

    logical :: used
    integer :: k

    if(.not.initialized) call mpp_error(FATAL,'radiation_mod: module not initialized')
    
    call get_pressure_at_levels(p,plvl,plyr)

    plvl = plvl * 10.0 !cb to mb
    plyr = plyr * 10.0 !cb to mb

    used = send_data(id_plvl,plvl,Time)
    used = send_data(id_plyr,plyr,Time)
    used = send_data(id_tskin,tskin,Time)
    used = send_data(id_albedo,sum(sfcalb,1)/NF_ALBD,Time)
    used = send_data(id_tlyr,tlyr,Time)
    used = send_data(id_coszen,coszen,Time)
    used = send_data(id_coszdg,coszdg,Time)

    tlvl(nlevp1,:,:) = tlyr(nlev,:,:)

    tlvl(1,:,:) = tskin(:,:)
    
    tem2db(nlevp1,:,:) = 0.
    tem2db(1:nlev,:,:) = log(plvl(1:nlev,:,:))
    tem2da = log(plyr)

    do k = 1, nlev-1
        tlvl(k+1,:,:) = tlyr(k,:,:) + (tlyr(k+1,:,:) - tlyr(k,:,:)) &
                                    * (tem2db(k+1,:,:) - tem2da(k,:,:)) &
                                    / (tem2da(k+1,:,:) - tem2da(k,:,:))
    enddo

    used = send_data(id_tlvl,tlvl,Time)

    qlyr = 0.
    olyr = 0.
    clw = 0.
    if (ind_q>0) qlyr = tr(:,:,:,ind_q)
    if (ind_oz>0) olyr = tr(:,:,:,ind_oz)
    if (ind_clw>0) clw = tr(:,:,:,ind_clw)

    faersw = 0.

    call data_override('ATM','ozone',olyr,Time,kxy=1)

    used = send_data(id_ozone,olyr,Time)

    !call get_gases(Time,gasvmr)

    !call get_clouds(plyr,plvl,tlyr,qlyr,clw,clouds,Time)

    !if (isubc==2) call get_icseed(Time,icseed) 

    !call swrad_drv(plyr,plvl,tlyr,tlvl,qlyr, &
    !               olyr,gasvmr,clouds,icseed(:,:,1),faersw,sfcalb,coszen,solcon, &
    !               hswc,rsdt,rsut,rsutc,rsds,rsdsc,rsus,rsusc,ruvds, &
    !               ruvdsc,rnirdsbm,rnirdsdf,rvisdsbm,rvisdsdf)

end subroutine radiation


!--------------------------------------------------------------------------------   
subroutine get_icseed(Time,icseed)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    integer, intent(out) :: icseed(:,:,:)

    integer :: numrdm(size(icseed,1)*size(icseed,2)*size(icseed,3))
    integer :: ipseed, jdat(6)
    integer, parameter :: ipsdlim = 1.0E8
    integer, parameter :: ipsd0 = 17*2009 + 43*12 + 37*12 + 23*0
    type(random_stat) :: stat

    call get_date(Time,jdat(1),jdat(2),jdat(3),jdat(4),jdat(5),jdat(6))

    ipseed = mod(nint(100.0*sqrt(real(sum(jdat))*3600)), ipsdlim) + 1 + ipsd0

    call random_setseed(ipseed,stat)
    call random_index(ipsdlim,numrdm,stat)

    icseed = reshape(numrdm,[size(icseed,1),size(icseed,2),size(icseed,3)])

end subroutine get_icseed


!--------------------------------------------------------------------------------   
subroutine get_clouds(plyr,plvl,tlyr,qlyr,clw,clouds,Time)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, dimension(:,js:,is:), intent(in) :: plyr, plvl, tlyr, qlyr, clw
    real, dimension(:,:,js:,is:), intent(out) :: clouds

    real, dimension(5,js:je,is:ie) :: clds
    integer, dimension(3,js:je,is:ie) :: mtop, mbot
    real, dimension(size(plyr,1))  :: qstl, rhly
    integer                        :: i, j, nlay


    nlay = size(plyr,1)
 
    do i = is, ie
        do j = js, je
            call compute_qs(tlyr(:,j,i),100*plyr(:,j,i),qstl(:))
            rhly(:) = qlyr(:,j,i)/qstl(:)

            call progcld1(plyr(:,j,i),plvl(:,j,i),tlyr(:,j,i),qlyr(:,j,i), &
                          qstl(:),rhly(:),clw(:,j,i),lat_rad(j,i),fland(j,i),1,nlay, &
                          nlay+1, iovr, .false., .false., .true., clouds(:,:,j,i), &
                          clds(:,j,i), mtop(:,j,i), mbot(:,j,i)) 
        enddo
    enddo
      

end subroutine get_clouds


!--------------------------------------------------------------------------------   
subroutine get_gases(Time,gasvmr)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(out) :: gasvmr(:,:,:,:)

    call data_override('ATM','co2vmr',co2vmr,Time)
    call data_override('ATM','n2ovmr',n2ovmr,Time)
    call data_override('ATM','ch4vmr',ch4vmr,Time)
    call data_override('ATM','o2vmr',o2vmr,Time)
    call data_override('ATM','covmr',covmr,Time)
    call data_override('ATM','f11vmr',f11vmr,Time)
    call data_override('ATM','f12vmr',f12vmr,Time)
    call data_override('ATM','f22vmr',f22vmr,Time)
    call data_override('ATM','cl4vmr',cl4vmr,Time)
    call data_override('ATM','f113vmr',f113vmr,Time)

    gasvmr(:,:,:,1)  = co2vmr
    gasvmr(:,:,:,2)  = n2ovmr
    gasvmr(:,:,:,3)  = ch4vmr
    gasvmr(:,:,:,4)  = o2vmr
    gasvmr(:,:,:,5)  = covmr
    gasvmr(:,:,:,6)  = f11vmr
    gasvmr(:,:,:,7)  = f12vmr
    gasvmr(:,:,:,8)  = f22vmr
    gasvmr(:,:,:,9)  = cl4vmr
    gasvmr(:,:,:,10) = f113vmr

end subroutine get_gases

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
    real, intent(in) :: gasvmr(nlev,NF_VGAS,imax), clouds(nlev,NF_CLDS,imax)
    integer, intent(in) :: icseed(imax)
    real, intent(in) :: aerosols(nlev,NBDSW,NF_AESW,imax)
    real, intent(in) :: sfcalb(NF_ALBD,imax), cosz(imax), solcon
    real, intent(out) :: hswc(nlev,imax)
    real, dimension(imax), intent(out) :: topdnfxc,topupfxc,topupfx0
    real, dimension(imax), intent(out) :: sfcdnfxc,sfcdnfx0,sfcupfxc,sfcupfx0
    real, dimension(imax), intent(out) :: sfcuvbfc,sfcuvbf0,sfcnirbm,sfcnirdf
    real, dimension(imax), intent(out) :: sfcvisbm,sfcvisdf

    integer :: i


    do i = 1, imax
        call swrad(nlev, nlevp1, plyr(:,i), plvl(:,i), tlyr(:,i), tlvl(:,i), &
                qlyr(:,i), olyr(:,i), gasvmr(:,:,i), clouds(:,:,i), icseed(i), &
                aerosols(:,:,:,i), sfcalb(:,i), cosz(i), solcon, &
                hswc(:,i), topdnfxc(i), topupfxc(i), topupfx0(i), &
                sfcdnfxc(i), sfcdnfx0(i), sfcupfxc(i), sfcupfx0(i), &
                sfcuvbfc(i), sfcuvbf0(i), sfcnirbm(i), sfcnirdf(i), &
                sfcvisbm(i), sfcvisdf(i))

    enddo


end subroutine swrad_drv


end module radiation_mod
