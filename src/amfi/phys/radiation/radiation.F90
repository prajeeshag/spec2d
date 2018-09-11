module radiation_mod

use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, &
                    mpp_clock_id, mpp_clock_begin, mpp_clock_end

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain

use fms_mod, only : open_namelist_file, close_file, read_data, write_data, field_size

use time_manager_mod, only : time_type, set_time, get_date

use data_override_mod, only : data_override

use diag_manager_mod, only : diag_axis_init, reg_df=>register_diag_field, send_data

use sat_vapor_pres_mod, only : compute_qs

use constants_mod, only : PI

use mersenne_twister, only : random_setseed, random_index, random_stat

use module_radiation_clouds, only : progcld1

use module_radsw_main, only : init_swrad, swrad, NBDSW

use module_radlw_main, only : init_lwrad, lwrad, NBDLW=>nbands

use set_aerosols_mod, only : init_set_aerosols, set_aerosols

#ifdef AQUAPLANET
use aqua_planet_mod, only : aquape_o3
#endif

implicit none
private

public :: init_radiation, radiation

type(time_type) :: dt, rad_dt

type(domain2D), pointer :: domain

integer :: is, ie, js, je, nlev, ntrac, nlevp1
integer :: ilen, jlen, imax

real, public, parameter :: con_solr = 1.3610e+3      ! solar constant (W/m2)-liu(2002)

integer, public, parameter :: NF_AESW=3, NF_AELW=3
integer, public, parameter :: NF_CLDS=9, NF_VGAS=10
integer, public, parameter :: NF_ALBD=4
integer :: icwp=1, iovr=1, isubc=2
character(len=512) :: ozone_fnm="NOOZON" ! ozone file name, if no ozone file 
                                   ! name given, assumes no ozone

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

real, allocatable, dimension(:,:) :: lat_rad, lon_rad, lat_deg, lon_deg, fland
real, allocatable, dimension(:,:,:) :: o3tmp
real, allocatable, dimension(:) :: o3lev

integer :: id_plvl, id_plyr, id_tskin, id_albedo, id_tlyr, id_tlvl, id_coszen, &
           id_fracday, id_ozone, id_qlyr, id_clw, id_clouds(NF_CLDS), id_gasvmr(NF_VGAS), &
           id_icseed(2), id_sfcemis
integer :: id_htsw, id_rsdt, id_rsut, id_rsutc, id_rsds, id_rsdsc, id_rsus, id_rsusc, &
           id_ruvds, id_ruvdsc, id_rnirdsbm, id_rnirdsbf, id_rvisdsbm, id_rvisdsdf 
integer :: id_htlw, id_rlut, id_rlutc, id_rlds, id_rldsc, id_rlus, id_rlusc, &
           id_aod, id_asy, id_ssa

integer :: clck_swrad, clck_lwrad

logical :: initialized=.false., debug=.false.

namelist/radiation_nml/icwp, iovr, isubc, ozone_fnm, debug

contains

!--------------------------------------------------------------------------------   
subroutine init_radiation(Time,domain_in,ntrac_in,nlev_in, &
                          lat_deg_in,lon_deg_in,fland_in, axes_in, &
                          ind_q_in,ind_clw_in,ind_oz_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2D), target :: domain_in
    integer, intent(in) :: nlev_in, ntrac_in
    real, intent(in) :: lat_deg_in(:,:), lon_deg_in(:,:)
    real, intent(in) :: fland_in(:,:)
    integer, intent(in) :: axes_in(:)
    integer, intent(in), optional :: ind_q_in, ind_oz_in, ind_clw_in
    integer :: unit, i, j, me
  
    unit = open_namelist_file()

    read(unit,nml=radiation_nml)

    call close_file(unit)

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

    lat_deg(js:je,is:ie) = lat_deg_in(js:je,is:ie)
    lat_rad(js:je,is:ie) = lat_deg_in(js:je,is:ie)*PI/180.

    lon_deg(js:je,is:ie) = lon_deg_in(js:je,is:ie)
    lon_rad(js:je,is:ie) = lon_deg_in(js:je,is:ie)*PI/180.

    me = 1
    if (mpp_pe()==mpp_root_pe()) me = 0

    call init_swrad(icwp,me,iovr,isubc)
    call init_lwrad(icwp,me,iovr,isubc)

    call rad_diag_init(Time,axes_in)

    call init_set_aerosols(is,ie,js,je,nlev)

    clck_swrad = mpp_clock_id('SW radiation')
    clck_lwrad = mpp_clock_id('LW radiation')

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
    character(len=8) :: fldnm, rou='am_rad'

    allocate(ilevs(nlev)) 
    allocate(ilevsp(nlev+1)) 

    forall(i=1:nlev) ilevs(i) = i
    forall(i=1:nlev+1) ilevsp(i) = i

    id_lat = axes_in(1) 
    id_lon = axes_in(2) 
    id_lev = axes_in(3)
    id_levp = axes_in(4) 

    axes = [id_lev,id_lat,id_lon] 
    axesp = [id_levp,id_lat,id_lon] 
    
    id_aod = reg_df(rou, 'aod', axes, Time, 'Aerosol Optical Depth at 550 nm', ' ')
    id_ssa = reg_df(rou, 'ssa', axes, Time, 'Single scattering Albedo', ' ')
    id_asy = reg_df(rou, 'asy', axes, Time, 'Asymmetry Parameter', ' ')

    id_plvl = reg_df(rou,'plvl',axesp,Time,'Pressure','mb')
    id_plyr = reg_df(rou,'plyr',axes,Time,'Pressure','mb')
    id_tskin = reg_df(rou,'tskin',axes(2:3),Time,'Surface Temperature','K')
    id_albedo = reg_df(rou,'albedo',axes(2:3),Time,'Surface Albedo 1','1')
    id_tlyr = reg_df(rou,'tlyr',axes,Time,'Temperature','K')
    id_tlvl = reg_df(rou,'tlvl',axesp,Time,'Temperature','K')
    id_coszen = reg_df(rou,'coszen',axes(2:3),Time,'coszen','1')
    id_fracday = reg_df(rou,'fracday',axes(2:3),Time,'fracday','1')
    id_ozone = reg_df(rou,'ozone',axes,Time,'Ozone','?')
    id_qlyr = reg_df(rou,'qlyr',axes,Time,'Sp Humidity','?')
    id_clw = reg_df(rou,'clw',axes,Time,'Cloud Condensate','?')

    id_sfcemis = reg_df(rou,'sfcemis',axes(2:3),Time,'Surface Emissivity','?')

    do i = 1, NF_CLDS
        write(fldnm,'(A,I2.2)') 'clouds',i
        id_clouds(i) = reg_df(rou,fldnm,axes,Time,fldnm,'?')
    enddo

    do i = 1, NF_VGAS
        write(fldnm,'(A,I2.2)') 'gas',i
        id_gasvmr(i) = reg_df(rou,fldnm,axes,Time,fldnm,'?')
    enddo

    do i = 1, 2
        write(fldnm,'(A,I2.2)') 'icsd',i
        id_icseed(i) = reg_df(rou,fldnm,axes,Time,fldnm,'?')
    enddo

    !Shortwave ----------------------------------------------    
    id_htsw = reg_df(rou, 'htsw', axes, Time,  &
              'Shortwave Heating Rate', '?')
    id_rsdt = reg_df(rou, 'rsdt', axes(2:3), Time,  &
              'TOA Incident Shortwave Radiation',  'W m-2')
    id_rsut = reg_df(rou, 'rsut', axes(2:3), Time,  &
              'TOA Upward Shortwave',  'W m-2')
    id_rsutc = reg_df(rou, 'rsutc', axes(2:3), Time, &
              'TOA Upward Shortwave (Clear)',  'W m-2')
    id_rsds = reg_df(rou, 'rsds', axes(2:3), Time, &
              'Surface Downwelling Shortwave',  'W m-2')
    id_rsdsc = reg_df(rou, 'rsdsc', axes(2:3), Time, &
              'Surface Downwelling Shortwave (Clear)',  'W m-2')
    id_rsus = reg_df(rou, 'rsus', axes(2:3), Time, &
              'Surface Upwelling Shortwave',  'W m-2')
    id_rsusc = reg_df(rou, 'rsusc', axes(2:3), Time, &
              'Surface Upwelling Shortwave (Clear)',  'W m-2')
    id_ruvds = reg_df(rou, 'ruvds', axes(2:3), Time, &
              'Surface Downwelling UV',  'W m-2')
    id_ruvdsc = reg_df(rou, 'ruvdsc', axes(2:3), Time, &
              'Surface Downwelling UV (Clear)',  'W m-2')
    id_rnirdsbm = reg_df(rou, 'rnirdsbm', axes(2:3), Time, &
              'Surface Downwelling Near-IR (dir)',  'W m-2')
    id_rnirdsbf = reg_df(rou, 'rnirdsdf', axes(2:3), Time, &
              'Surface Downwelling Near-IR (dif)',  'W m-2')
    id_rvisdsbm = reg_df(rou, 'rvisdsbm', axes(2:3), Time, &
              'Surface Downwelling visible (dir)',  'W m-2')
    id_rvisdsdf = reg_df(rou, 'rvisdsdf ', axes(2:3), Time, &
              'Surface Downwelling visible (dif)',  'W m-2')
    !-----------------------------------------------------------
    
    !Longwave ----------------------------------------------    
    id_htlw = reg_df(rou, 'htlw', axes, Time,  &
              'Longwave Heating Rate', 'K sec-1')
    id_rlut = reg_df(rou, 'rlut', axes(2:3), Time,  &
              'TOA Upward Longwave',  'W m-2')
    id_rlutc = reg_df(rou, 'rlutc', axes(2:3), Time, &
              'TOA Upward Longwave (Clear)',  'W m-2')
    id_rlds = reg_df(rou, 'rlds', axes(2:3), Time, &
              'Surface Downwelling Longwave',  'W m-2')
    id_rldsc = reg_df(rou, 'rldsc', axes(2:3), Time, &
              'Surface Downwelling Longwave (Clear)',  'W m-2')
    id_rlus = reg_df(rou, 'rlus', axes(2:3), Time, &
              'Surface Upwelling Longwave',  'W m-2')
    id_rlusc = reg_df(rou, 'rlusc', axes(2:3), Time, &
              'Surface Upwelling Longwave (Clear)',  'W m-2')
    !-----------------------------------------------------------
end subroutine rad_diag_init


!--------------------------------------------------------------------------------   
subroutine radiation(Time, tlyr, tr, prsl, prsi, phil, oro, tskin, coszen, fracday, &
                     sfcalb, sfcemis, solcon, htsw, rsds, rsus, htlw, rlds, rlus)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(1:nlev,js:je,is:ie) :: prsl, tlyr, phil
    real, intent(in), dimension(1:nlev+1,js:je,is:ie) :: prsi
    real, intent(in)  :: tr(1:nlev,js:je,is:ie,1:ntrac)
    real, intent(in), dimension(js:je,is:ie) :: tskin, coszen, fracday, sfcemis, oro
    real, intent(in), dimension(NF_ALBD,js:je,is:ie) :: sfcalb
    real, intent(in) :: solcon

    real, intent(out), dimension(1:nlev,js:je,is:ie) :: htsw
    real, intent(out), dimension(1:nlev,js:je,is:ie) :: htlw
    real, intent(out), dimension(js:je,is:ie) :: rsds, rsus
    real, intent(out), dimension(js:je,is:ie) :: rlds, rlus

    real, dimension(1:nlev,js:je,is:ie) :: plyr, qlyr, olyr, clw
    real, dimension(1:nlevp1,js:je,is:ie) :: plvl, tlvl

    real, dimension(1:nlev,NF_CLDS,js:je,is:ie) :: clouds
    real, dimension(1:nlev,NF_VGAS,js:je,is:ie) :: gasvmr
    real, dimension(1:nlev,NBDSW,NF_AESW,js:je,is:ie) :: faersw
    real, dimension(1:nlev,NBDLW,NF_AELW,js:je,is:ie) :: faerlw

    real, dimension(js:je,is:ie) :: rsdt, rsut, rsutc
    real, dimension(js:je,is:ie) :: rlut, rlutc
    real, dimension(js:je,is:ie) :: rsdsc, rsusc
    real, dimension(js:je,is:ie) :: rldsc, rlusc
    real, dimension(js:je,is:ie) :: ruvds, ruvdsc
    real, dimension(js:je,is:ie) :: rnirdsbm, rnirdsdf
    real, dimension(js:je,is:ie) :: rvisdsbm, rvisdsdf

    integer, dimension(js:je,is:ie,2) :: icseed

    real :: tem2db(nlevp1,js:je,is:ie), tem2da(nlev,js:je,is:ie)

    logical :: used
    integer :: k, i, j

    if(.not.initialized) call mpp_error(FATAL,'radiation_mod: module not initialized')
    
    plvl = prsi * 10.0 !cb to mb
    plyr = prsl * 10.0 !cb to mb

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

    qlyr = 0.
    olyr = 0.
    clw = 0.
   

    if (ind_q>0) qlyr = tr(:,:,:,ind_q)

    !olyr in mass mixing ratio
    if (ind_oz>0) then
        olyr = tr(:,:,:,ind_oz)
    else
#ifndef AQUAPLANET
        call get_ozone(time,plyr,olyr) 
#else
        do i = is, ie
            do j = js, je
                call aquape_o3(lat_deg(j,i),plyr(:,j,i),olyr(:,j,i))
            end do
        end do
#endif
    endif
    if (debug) call mpp_error(NOTE,"radiation: after get_ozone")
        
    if (ind_clw>0) clw = tr(:,:,:,ind_clw)

    faersw = 0.
    faerlw = 0.

    call set_aerosols(Time, prsi, prsl, tlyr, phil, oro, faersw, faerlw)
    if (debug) call mpp_error(NOTE,"radiation: after set_aerosols")

    if (id_aod>0) used = send_data(id_aod, faersw(:,10,1,:,:), Time)
    if (id_ssa>0) used = send_data(id_ssa, faersw(:,10,2,:,:), Time)
    if (id_asy>0) used = send_data(id_asy, faersw(:,10,3,:,:), Time)

    call get_gases(Time,gasvmr)

    call get_clouds(plyr,plvl,tlyr,qlyr,clw,clouds,Time)
    if (debug) call mpp_error(NOTE,"radiation: after get_clouds")

    if (isubc==2) call get_icseed(Time,icseed) 

    call swrad_drv(plyr, plvl, tlyr, tlvl, qlyr, olyr, gasvmr, clouds, icseed(:,:,1), faersw, &
                   sfcalb, coszen, solcon, htsw, rsdt, rsut, rsutc, rsds, rsdsc, rsus, rsusc, &
                   ruvds, ruvdsc, rnirdsbm, rnirdsdf, rvisdsbm, rvisdsdf)
    if (debug) call mpp_error(NOTE,"radiation: after swrad_drv")

    call lwrad_drv(plyr, plvl, tlyr, tlvl, qlyr, olyr, gasvmr, clouds, icseed(:,:,2), faerlw, &
                   sfcemis, tskin, htlw, rlut, rlutc, rlds, rldsc, rlus, rlusc)
    if (debug) call mpp_error(NOTE,"radiation: after lwrad_drv")

    used = send_data(id_ozone,   olyr, Time)
    used = send_data(id_tlvl,    tlvl, Time)
    used = send_data(id_plvl,    plvl, Time)
    used = send_data(id_plyr,    plyr, Time)
    used = send_data(id_tskin,   tskin, Time)
    if (id_albedo>0) used = send_data(id_albedo,  sum(sfcalb,1)/NF_ALBD, Time)
    used = send_data(id_tlyr,    tlyr, Time)
    used = send_data(id_coszen,  coszen, Time)
    used = send_data(id_fracday, fracday, Time)
    used = send_data(id_qlyr,    qlyr, Time)
    used = send_data(id_clw,     clw, Time)

    ! Shortwave ------------------------------
    if (id_htsw > 0)     used = send_data(id_htsw,     htsw,             Time)
    if (id_rsdt > 0)     used = send_data(id_rsdt,     fracday*rsdt,     Time)
    if (id_rsut > 0)     used = send_data(id_rsut,     fracday*rsut,     Time)
    if (id_rsutc > 0)    used = send_data(id_rsutc,    fracday*rsutc,    Time)
    if (id_rsds > 0)     used = send_data(id_rsds,     fracday*rsds,     Time)
    if (id_rsdsc > 0)    used = send_data(id_rsdsc,    fracday*rsdsc,    Time)
    if (id_rsus > 0)     used = send_data(id_rsus,     fracday*rsus,     Time)
    if (id_rsusc > 0)    used = send_data(id_rsusc,    fracday*rsusc,    Time)
    if (id_ruvds > 0)    used = send_data(id_ruvds,    fracday*ruvds,    Time)
    if (id_ruvdsc > 0)   used = send_data(id_ruvdsc,   fracday*ruvdsc,   Time)
    if (id_rnirdsbm > 0) used = send_data(id_rnirdsbm, fracday*rnirdsbm, Time)
    if (id_rnirdsbf > 0) used = send_data(id_rnirdsbf, fracday*rnirdsdf, Time)
    if (id_rvisdsbm > 0) used = send_data(id_rvisdsbm, fracday*rvisdsbm, Time)
    if (id_rvisdsdf > 0) used = send_data(id_rvisdsdf, fracday*rvisdsdf, Time) 
    !------------------------------------------

    ! Longwave ------------------------------
    if (id_htlw > 0)     used = send_data(id_htlw,     htlw,     Time)
    if (id_rlut > 0)     used = send_data(id_rlut,     rlut,     Time)
    if (id_rlutc > 0)    used = send_data(id_rlutc,    rlutc,    Time)
    if (id_rlds > 0)     used = send_data(id_rlds,     rlds,     Time)
    if (id_rldsc > 0)    used = send_data(id_rldsc,    rldsc,    Time)
    if (id_rlus > 0)     used = send_data(id_rlus,     rlus,     Time)
    if (id_rlusc > 0)    used = send_data(id_rlusc,    rlusc,    Time)
    !------------------------------------------

    if (id_sfcemis > 0)    used = send_data(id_sfcemis,    sfcemis,    Time)

    do i = 1, NF_CLDS
        if(id_clouds(i) > 0) used = send_data(id_clouds(i), clouds(:,i,:,:), Time)
    enddo

    do i = 1, NF_VGAS
        if(id_gasvmr(i) > 0) used = send_data(id_gasvmr(i), gasvmr(:,i,:,:), Time)
    enddo

    do i = 1, 2
        if(id_icseed(i) > 0) used = send_data(id_icseed(i), icseed(:,:,i)*1., Time)
    enddo

end subroutine radiation

!--------------------------------------------------------------------------------   
subroutine get_ozone(time,prsl,o3)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: time
    real, intent(in) :: prsl(:,js:,is:)
    real, intent(out) :: o3(:,js:,is:)
    integer :: siz(4), i, j

    o3 = 0.
    if (trim(ozone_fnm)=='NOOZON') return

    if (.not.allocated(o3tmp)) then
        call field_size(trim(ozone_fnm),'level',siz)
        allocate(o3lev(siz(1)))
        allocate(o3tmp(siz(1),js:je,is:ie))
        call read_data(trim(ozone_fnm),'level',o3lev)
    end if
   
    call data_override('ATM','ozone',o3tmp,Time,kxy=1)

    do i = is, ie
        do j = js, je
            call interp_vert(o3tmp(:,j,i),o3(:,j,i),o3lev,prsl(:,j,i),.false.)
        end do
    end do    
   
    return 
end subroutine get_ozone


!--------------------------------------------------------------------------------   
subroutine interp_vert (fldin,fldout,axin,axout,extrap)
!--------------------------------------------------------------------------------   
    real, dimension(:), intent(in) :: fldin, axin, axout
    real, dimension(:), intent(out) :: fldout
    logical, intent(in) :: extrap

    integer :: i1, i2, ni, j
    real :: w1, w2, tmp(size(axin)), w12
   
    ni = size(axin) 

    do j = 1, size(axout)
        tmp = axin-axout(j)
        i1 = 0; i2 = 0
        if (any(tmp>=0)) i1 = minloc(tmp, 1, mask=tmp>=0.)
        if (any(tmp<=0)) i2 = maxloc(tmp, 1, mask=tmp<=0.)

        if(i1<=0.and.i2<=0) call mpp_error(FATAL,'radiation_mod:interp_vert: both i1 and i2 <= 0') 
        if (i1<=0) i1=i2
        if (i2<=0) i2=i1

        w1 = abs(axout(j)-axin(i2))
        w2 = abs(axout(j)-axin(i1))

        if (w1==0.or.w2==0.) then
            if (extrap) then
                w1 = 1.; w2 = 1.
            else
                fldout(j) = 0.
                cycle 
            endif
        endif

        w12 = w1 + w2

        w1 = w1/w12; w2=w2/w12

        fldout(j) = fldin(i1)*w1+fldin(i2)*w2
    end do

    return

end subroutine interp_vert


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

    gasvmr(:,1,:,:)  = co2vmr
    gasvmr(:,2,:,:)  = n2ovmr
    gasvmr(:,3,:,:)  = ch4vmr
    gasvmr(:,4,:,:)  = o2vmr
    gasvmr(:,5,:,:)  = covmr
    gasvmr(:,6,:,:)  = f11vmr
    gasvmr(:,7,:,:)  = f12vmr
    gasvmr(:,8,:,:)  = f22vmr
    gasvmr(:,9,:,:)  = cl4vmr
    gasvmr(:,10,:,:) = f113vmr

end subroutine get_gases

!--------------------------------------------------------------------------------   
subroutine swrad_drv(plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr, &
                     clouds,icseed,aerosols,sfcalb,cosz,solcon, &
                     htsw,rsdt,rsut,rsutc,rsds,rsdsc,rsus,rsusc, &
                     ruvds,ruvdsc,rnirbm,rnirdf,rvisbm,rvisdf)
!--------------------------------------------------------------------------------   

    real, dimension(nlev,imax), intent(in) :: plyr, tlyr, qlyr, olyr 
    real, dimension(nlevp1,imax), intent(in) :: plvl, tlvl
    real, intent(in) :: gasvmr(nlev,NF_VGAS,imax), clouds(nlev,NF_CLDS,imax)
    integer, intent(in) :: icseed(imax)
    real, intent(in) :: aerosols(nlev,NBDSW,NF_AESW,imax)
    real, intent(in) :: sfcalb(NF_ALBD,imax), cosz(imax), solcon
    real, intent(out) :: htsw(nlev,imax)
    real, dimension(imax), intent(out) :: rsdt, rsut, rsutc, rsds, rsdsc, rsus, rsusc
    real, dimension(imax), intent(out) :: ruvds, ruvdsc, rnirbm, rnirdf, rvisbm, rvisdf

    integer :: i

    call mpp_clock_begin(clck_swrad)

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(i)
    do i = 1, imax
        call swrad(plyr(:,i), plvl(:,i), tlyr(:,i), tlvl(:,i), &
                qlyr(:,i), olyr(:,i), gasvmr(:,:,i), clouds(:,:,i), icseed(i), &
                aerosols(:,:,:,i), sfcalb(:,i), cosz(i), solcon, &
                nlev, nlevp1, NF_VGAS, NF_CLDS, NF_AESW, NF_ALBD, &
                htsw(:,i),rsdt(i),rsut(i),rsutc(i),rsds(i),rsdsc(i),rsus(i),rsusc(i), &
                ruvds(i),ruvdsc(i),rnirbm(i),rnirdf(i),rvisbm(i),rvisdf(i))
    enddo
    !$OMP END PARALLEL DO

    call mpp_clock_end(clck_swrad)

    return

end subroutine swrad_drv

!--------------------------------------------------------------------------------   
subroutine lwrad_drv(plyr, plvl, tlyr, tlvl, qlyr, olyr, gasvmr, clouds, icseed, aerosols, &
                     sfcemis, tskin, htlw, rlut, rlutc, rlds, rldsc, rlus, rlusc)
!--------------------------------------------------------------------------------  
    implicit none
     
    real, dimension(nlev,imax), intent(in) :: plyr, tlyr, qlyr, olyr 
    real, dimension(nlevp1,imax), intent(in) :: plvl, tlvl
    real, intent(in) :: gasvmr(nlev,NF_VGAS,imax), clouds(nlev,NF_CLDS,imax)
    integer, intent(in) :: icseed(imax)
    real, intent(in) :: aerosols(nlev,NBDSW,NF_AELW,imax)
    real, intent(in) :: sfcemis(imax), tskin(imax)
    real, intent(out) :: htlw(nlev,imax)
    real, dimension(imax), intent(out) :: rlut, rlutc, rlds, rldsc, rlus, rlusc

    integer :: i

    call mpp_clock_begin(clck_lwrad)

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(i)
    do i = 1, imax
        call lwrad(plyr(:,i), plvl(:,i), tlyr(:,i), tlvl(:,i), qlyr(:,i), olyr(:,i), gasvmr(:,:,i), &
                   clouds(:,:,i), icseed(i), aerosols(:,:,:,i), sfcemis(i), tskin(i), &
                   nlev, nlevp1, NF_VGAS, NF_CLDS, NF_AELW, &
                   htlw(:,i), rlut(i), rlutc(i), rlds(i), rldsc(i), rlus(i), rlusc(i))
    enddo
    !$OMP END PARALLEL DO

    call mpp_clock_end(clck_lwrad)
    return
end subroutine lwrad_drv

end module radiation_mod
