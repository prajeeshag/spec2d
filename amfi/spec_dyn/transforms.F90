module transforms_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_sum
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_global_domain

use fms_mod, only : open_namelist_file, close_file, fms_init
use fms_mod, only : read_data, write_data

use fms_io_mod, only : read_data, restart_file_type, register_restart_field
use fms_io_mod, only : restore_state, save_restart
use fms_io_mod, only : fms_io_exit 

use constants_mod, only : RADIUS

use grid_fourier_mod, only : init_grid_fourier, fft_1dr2c_serial, fft_1dc2c_serial
use grid_fourier_mod, only : end_grid_fourier, grid_to_fourier, fourier_to_grid

use spherical_mod, only : init_spherical, fourier_to_spherical
use spherical_mod, only : spherical_to_fourier

use spherical_mod, only : compute_lon_deriv_cos, compute_lat_deriv_cos
use spherical_mod, only : compute_vor_div, compute_ucos_vcos
use spherical_mod, only : ev, od, do_truncation, get_wdecomp, get_spherical_wave
use spherical_mod, only : get_lats


implicit none
private

public :: compute_ucos_vcos, compute_vor_div, get_lats
public :: get_wdecomp, get_spherical_wave, get_lons
public :: spherical_to_grid, grid_to_spherical, init_transforms
public :: register_spec_restart, restore_spec_state, save_spec_restart

type(domain2d) :: domainf

integer :: ilen, jlen
integer :: nwaves_oe=0
integer :: nwaves_oe_total=0

logical :: fpe=.false. 

integer, allocatable :: pelist(:), extent(:)
integer, allocatable :: fpelist(:), fextent(:)
integer, allocatable :: Tshuff(:), spextent(:)

real, allocatable :: cosm_lat(:), cosm2_lat(:)
real, allocatable :: deg_lat(:), deg_lon(:)

real :: lon_start = 0.
real :: dlon = 0.

integer :: comm, n
real :: ref_temp
integer :: isc, iec, m, l, t, i, k
integer :: isf, ief, flen, jsc, jec, j, jsf, jef
integer :: trunc, nlon, nlat
integer :: pe
integer :: ishuff=0

integer, parameter :: max_num_dom=10
integer :: num_dom = 0
type(domain2d) :: spresdom(max_num_dom)
type(restart_file_type) :: specres

logical :: initialized=.false.

interface vor_div_to_uv_grid
    module procedure vor_div_to_uv_grid3d
    module procedure vor_div_to_uv_grid2d
end interface vor_div_to_uv_grid

interface vor_div_from_uv_grid
    module procedure vor_div_from_uv_grid3D
    module procedure vor_div_from_uv_grid2D
end interface vor_div_from_uv_grid

interface uv_grid_to_vor_div
    module procedure uv_grid_to_vor_div3D
    module procedure uv_grid_to_vor_div2D
end interface uv_grid_to_vor_div

interface grid_to_spherical
    module procedure grid_to_spherical3D
    module procedure grid_to_spherical2D
end interface

interface spherical_to_grid
    module procedure spherical_to_grid3D
    module procedure spherical_to_grid2D
end interface

contains

!--------------------------------------------------------------------------------   
subroutine init_transforms(domainl,trunc_in,nwaves)
!--------------------------------------------------------------------------------

    type(domain2d) :: domainl
    integer, intent(in) :: trunc_in
    integer, intent(out) :: nwaves
    integer :: isg, ieg, jsg, jeg

    call mpp_init() 
    call fms_init()

    trunc = trunc_in

    pe = mpp_pe()

    allocate(pelist(mpp_npes()))
    allocate(extent(mpp_npes()))

    call mpp_get_compute_domain(domainl, jsc, jec, isc, iec,ishuff=ishuff)
    call mpp_get_global_domain(domainl, jsg, jeg, isg, ieg)

    ilen = iec-isc+1
    jlen = jec-jsc+1

    nlon = ieg-isg+1
    nlat = jeg-jsg+1

    call mpp_get_current_pelist(pelist,commid=comm)

    allocate(Tshuff(0:trunc))
    call init_grid_fourier (nlon, ilen, trunc, isf, flen, comm, Tshuff)

    call mpp_gather([flen], extent)
    call mpp_broadcast(extent,size(extent), mpp_root_pe())

    allocate(fextent(count(extent>0)))
    allocate(spextent(count(extent>0)))
    allocate(fpelist(count(extent>0)))
   
    k = 0 
    do i = 1, size(extent)
        if (extent(i)>0) then
            k = k + 1
            fextent(k) = extent(i)
            fpelist(k) = pelist(i)
        endif
    enddo
    
    fpe = any(fpelist==mpp_pe())

    call mpp_declare_pelist(fpelist,'fourier_pes')
   
    isf = 0; ief = -1 
    jsf = 0; jef = -1 

    if (fpe) then

        call mpp_set_current_pelist(fpelist)

        call mpp_define_domains( [1,nlat, 0,trunc], [1,mpp_npes()], domainf, yextent=fextent)

        call mpp_get_compute_domain(domainf, jsf, jef, isf, ief)

        if(flen /= ief-isf+1) call mpp_error('gloopa', 'flen /= ief-isf+1', FATAL)

        call init_spherical(trunc, nlat, ishuff, nwaves_oe, domainf, Tshuff) 

        nwaves_oe_total=nwaves_oe

        call mpp_sum(nwaves_oe_total)

        call mpp_gather([nwaves_oe], spextent)

        call mpp_broadcast(spextent,size(spextent),mpp_root_pe())

    endif

    call mpp_set_current_pelist()

    nwaves = nwaves_oe

    allocate(cosm_lat(jsc:jec))
    allocate(cosm2_lat(jsc:jec))
    allocate(deg_lat(nlat))
    allocate(deg_lon(nlon))

    call get_lats(cosmlat=cosm_lat,cosm2lat=cosm2_lat,deglat=deg_lat)

    dlon = 360./nlon
    deg_lon(1) = lon_start
    do i = 2, nlon
       deg_lon(i) = deg_lon(i-1) + dlon
    enddo 

    initialized = .true.

    return
end subroutine init_transforms


!--------------------------------------------------------------------------------   
subroutine get_lons(deglon)
!--------------------------------------------------------------------------------   
    real, intent(out) :: deglon(:)

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')

    if (size(deglon,1)==nlon) then
        deglon = deg_lon
    elseif (size(deglon,1)==ilen) then
        deglon = deg_lon(isc:iec)
    else
        call mpp_error(FATAL,'get_lons: Wrong size, deglon should either be on' &
                             //' compute domain or global domain.')

    endif

end subroutine get_lons
    

!--------------------------------------------------------------------------------   
function register_spec_restart(filename,fieldname,data,mandatory,data_default)
!--------------------------------------------------------------------------------
    character (len=*), intent(in) :: filename, fieldname
    complex, intent(in) :: data(:,:,:)
    real, optional, intent(in) :: data_default
    logical, optional, intent(in) :: mandatory
    integer :: register_spec_restart 

    integer :: isize, jsize, idx_dom, gdom(4), layout(2)
    type(C_PTR) :: cptr
    real, pointer :: rdata(:,:,:) => NULL()
   
    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')

    register_spec_restart = 0

    if (.not.fpe) return

    idx_dom = 0

    do i = 1, num_dom
        call mpp_get_compute_domain(spresdom(i),xsize=isize)
        if (isize==size(data,1)*2) then
            idx_dom = i
            exit
        endif
    enddo

    if (idx_dom<1) then
        num_dom = num_dom + 1
        if (num_dom>max_num_dom) call mpp_error(FATAL,'register_spec_restart: num_dom>max_num_dom') 
        idx_dom = num_dom
        gdom = [1,size(data,1)*2,1,nwaves_oe_total]
        layout = [1,mpp_npes()]
        call mpp_define_domains(gdom, layout, spresdom(idx_dom),yextent=spextent)
    endif

        cptr = C_LOC(data)
        call c_f_pointer(cptr,rdata,[size(data,1)*2,size(data,2),size(data,3)])
        register_spec_restart = &
            register_restart_field(specres, filename, fieldname, rdata, &
                     spresdom(idx_dom), mandatory, data_default=data_default)

end function register_spec_restart

!--------------------------------------------------------------------------------   
subroutine restore_spec_state()
!--------------------------------------------------------------------------------   
    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    call restore_state(specres)
end subroutine restore_spec_state

!--------------------------------------------------------------------------------   
subroutine save_spec_restart(prefix)
!--------------------------------------------------------------------------------   
    character(len=*), optional :: prefix

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    call save_restart(specres,prefix)
end subroutine save_spec_restart 

!--------------------------------------------------------------------------------   
subroutine vor_div_to_uv_grid3d(vor,div,u,v,getcosuv)
!--------------------------------------------------------------------------------   
    complex,dimension(:,:,:), intent(in) :: vor, div
    real, intent(out) :: u(:,:,:), v(:,:,:)
    logical, intent(in), optional :: getcosuv
    complex,dimension(size(vor,1),size(vor,2),2) :: sucos, svcos
    logical :: getcosuv1

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    getcosuv1 = .false.

    if(present(getcosuv)) getcosuv1=getcosuv

    call compute_ucos_vcos(vor,div,sucos,svcos)

    call spherical_to_grid3D(sucos,grid=u)

    call spherical_to_grid3D(svcos,grid=v)

    if (getcosuv1) return

    do j = jsc, jec
        u(:,j,:) = u(:,j,:) * cosm_lat(j)
        v(:,j,:) = v(:,j,:) * cosm_lat(j)
    end do

    return

end subroutine vor_div_to_uv_grid3d

     
subroutine vor_div_to_uv_grid2d(vor,div,u,v,getcosuv)

    complex,dimension(:,:,:), intent(in) :: vor, div
    real, intent(out) :: u(:,:), v(:,:)
    logical, intent(in), optional :: getcosuv

    real :: u3d(1,size(u,1),size(u,2))
    real :: v3d(1,size(u,1),size(u,2))

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    call vor_div_to_uv_grid3d(vor, div, u3d, v3d,getcosuv)

    u(:,:) = u3d(1,:,:)
    v(:,:) = v3d(1,:,:)

    return

end subroutine vor_div_to_uv_grid2d

!--------------------------------------------------------------------------------   
subroutine vor_div_from_uv_grid3D(u,v,vor,div,uvcos_in)
!--------------------------------------------------------------------------------   
    real, intent(in) :: u(:,:,:), v(:,:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div
    logical, optional :: uvcos_in

    complex,dimension(size(vor,1),size(vor,2),2) :: usm, vsm, usp, vsp

    complex :: four(size(u,1)*size(u,2), flen)
    complex, pointer :: four3(:,:,:)

    real, pointer :: grd(:,:)
    type(C_PTR) :: pgrd, pfour

    real :: rradius=1./RADIUS
    integer :: nk, nj, ni, howmany
    integer :: i, j, k
    logical :: uvcos_in1

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    uvcos_in1=.false.

    if(present(uvcos_in)) uvcos_in1=uvcos_in

    nk = size(u,1); nj = size(u,2); ni = size(u,3)

    howmany = nk*nj

    pgrd = C_LOC(u)
    call c_f_pointer(pgrd, grd, [howmany, ni])

    pfour = C_LOC(four)
    call c_f_pointer(pfour, four3, [nk, nj, flen])

    call grid_to_fourier(grd,four)

    if (uvcos_in1) then
        do j = jsc, jec
            four3(:,j,:) = four3(:,j,:) * cosm2_lat(j)
        enddo
    else
        do j = jsc, jec
            four3(:,j,:) = four3(:,j,:) * cosm_lat(j)
        enddo
    endif
    
    call fourier_to_spherical(four3,usp,do_trunc=.true.)
    call compute_lon_deriv_cos(usp,usm)
    call fourier_to_spherical(four3,usp,useHnm=.true.,do_trunc=.true.)

    pgrd = C_LOC(v)
    call c_f_pointer(pgrd, grd, [howmany, ni])

    call grid_to_fourier(grd,four)

    if (uvcos_in1) then
        do j = jsc, jec
            four3(:,j,:) = four3(:,j,:) * cosm2_lat(j)
        enddo
    else
        do j = jsc, jec
            four3(:,j,:) = four3(:,j,:) * cosm_lat(j)
        enddo
    endif
    
    call fourier_to_spherical(four3,vsp,do_trunc=.true.)
    call compute_lon_deriv_cos(vsp,vsm)
    call fourier_to_spherical(four3,vsp,useHnm=.true.,do_trunc=.true.)

    rradius = 1./RADIUS
    vor = vsm + (usp*rradius)
    div = usm - (vsp*rradius)

end subroutine vor_div_from_uv_grid3D 

!--------------------------------------------------------------------------------   
subroutine uv_grid_to_vor_div3D(u,v,vor,div)
!--------------------------------------------------------------------------------   
    real, intent(in) :: u(:,:,:), v(:,:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div

    complex,dimension(size(vor,1),size(vor,2),2) :: us, vs

    complex :: four(size(u,1)*size(u,2), flen)
    complex, pointer :: four3(:,:,:)

    real, pointer :: grd(:,:)
    type(C_PTR) :: pgrd, pfour

    integer :: nk, nj, ni, howmany
    integer :: i, j, k

    if(.not.initialized) call mpp_error(FATAL,'transforms_mod: module not initialized!')
    nk = size(u,1); nj = size(u,2); ni = size(u,3)

    howmany = nk*nj

    pgrd = C_LOC(u)
    call c_f_pointer(pgrd, grd, [howmany, ni])

    pfour = C_LOC(four)
    call c_f_pointer(pfour, four3, [nk, nj, flen])

    call grid_to_fourier(grd,four)

    do j = jsc, jec
        four3(:,j,:) = four3(:,j,:) * cosm_lat(j)
    enddo
    
    call fourier_to_spherical(four3,us,do_trunc=.false.)

    pgrd = C_LOC(v)
    call c_f_pointer(pgrd, grd, [howmany, ni])

    call grid_to_fourier(grd,four)

    do j = jsc, jec
        four3(:,j,:) = four3(:,j,:) * cosm_lat(j)
    enddo
    
    call fourier_to_spherical(four3,vs,do_trunc=.false.)

    call compute_vor_div(us,vs,vor,div)

end subroutine uv_grid_to_vor_div3D 

!--------------------------------------------------------------------------------   
subroutine vor_div_from_uv_grid2D(u,v,vor,div)
!--------------------------------------------------------------------------------   
    real, intent(in) :: u(:,:), v(:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div
    real :: u3d(1,size(u,1),size(u,2))
    real :: v3d(1,size(u,1),size(u,2))

    u3d(1,:,:) = u(:,:)
    v3d(1,:,:) = v(:,:)
    call vor_div_from_uv_grid3D(u3d,v3d,vor,div)

end subroutine vor_div_from_uv_grid2D

!--------------------------------------------------------------------------------   
subroutine uv_grid_to_vor_div2D(u,v,vor,div)
!--------------------------------------------------------------------------------   
    real, intent(in) :: u(:,:), v(:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div
    real :: u3d(1,size(u,1),size(u,2))
    real :: v3d(1,size(u,1),size(u,2))

    u3d(1,:,:) = u(:,:)
    v3d(1,:,:) = v(:,:)
    call uv_grid_to_vor_div3D(u3d,v3d,vor,div)

end subroutine uv_grid_to_vor_div2D

!--------------------------------------------------------------------------------   
subroutine grid_to_spherical3D(grid,spherical,do_trunc)
!--------------------------------------------------------------------------------   
    real, intent(in) :: grid(:,:,:)
    complex, dimension(:,:,:), intent(out) :: spherical
    logical, intent(in), optional :: do_trunc
    complex :: four(size(grid,1)*size(grid,2), flen)
    complex, pointer :: four3(:,:,:)

    real, pointer :: grd(:,:)
    type(C_PTR) :: pgrd, pfour

    integer :: nk, nj, ni, howmany
    integer :: i, j, k
    logical :: do_trunc1
    
    do_trunc1 = .true.
    
    if(present(do_trunc)) do_trunc1 = do_trunc

    nk = size(grid,1); nj = size(grid,2); ni = size(grid,3)

    howmany = nk*nj

    pgrd = C_LOC(grid)
    call c_f_pointer(pgrd, grd, [howmany, ni])

    pfour = C_LOC(four)
    call c_f_pointer(pfour, four3, [nk, nj, flen])

    call grid_to_fourier(grd,four)

    call fourier_to_spherical(four3,spherical,do_trunc=do_trunc1)

    return
end subroutine grid_to_spherical3D 

!--------------------------------------------------------------------------------   
subroutine spherical_to_grid3D(spherical,grid,lat_deriv,lon_deriv)
!--------------------------------------------------------------------------------   
    complex,dimension(:,:,:), intent(in) :: spherical
    real, intent(out), optional :: grid(:,:,:)
    real, intent(out), optional :: lat_deriv(:,:,:)
    real, intent(out), optional :: lon_deriv(:,:,:)

    complex :: four(size(spherical,1)*jlen, isf:ief)
    complex, pointer :: four3(:,:,:)

    real, pointer :: grd(:,:)
    type(C_PTR) :: pgrd, pfour

    integer :: nk, nj, ni, howmany
    integer :: i, j, k
    real :: ma
    real, parameter :: rRADIUS = 1./RADIUS

    nk = size(spherical,1); nj = jlen; ni = iec-isc+1

    howmany = nk*nj

    pfour = C_LOC(four)
    call c_f_pointer(pfour, four3, [nk, nj, flen])

    if (present(lat_deriv)) then
        call spherical_to_fourier(spherical, four3, .true.)

        pgrd = C_LOC(lat_deriv)
        call c_f_pointer(pgrd, grd, [howmany, ni])

        call fourier_to_grid(four,grd)
    endif

    if (.not.present(lon_deriv).and. &
        .not.present(grid)) return

    call spherical_to_fourier(spherical, four3, .false.)

    if (present(grid)) then
        pgrd = C_LOC(grid)
        call c_f_pointer(pgrd, grd, [howmany, ni])

        call fourier_to_grid(four,grd)
    endif
     
    if (present(lon_deriv)) then
        pgrd = C_LOC(lon_deriv)
        call c_f_pointer(pgrd, grd, [howmany, ni])

        do m = isf, ief
            ma = Tshuff(m)*rRADIUS
            four(:,m) = ma*cmplx(-aimag(four(:,m)),real(four(:,m)))
        enddo
     
        call fourier_to_grid(four,grd)
    endif

    return
end subroutine spherical_to_grid3D

!--------------------------------------------------------------------------------   
subroutine grid_to_spherical2D(grid,spherical,do_trunc)
!--------------------------------------------------------------------------------   
    real, intent(in) :: grid(:,:)
    complex, dimension(:,:,:), intent(out) :: spherical
    logical, intent(in), optional :: do_trunc

    real :: buff(1,size(grid,1),size(grid,2))

    buff(1,:,:) = grid

    call grid_to_spherical3D(buff, spherical, do_trunc)

    return
end subroutine grid_to_spherical2D

!--------------------------------------------------------------------------------   
subroutine spherical_to_grid2D(spherical,grid,lat_deriv,lon_deriv)
!--------------------------------------------------------------------------------   
    complex,dimension(:,:), intent(in) :: spherical
    real, intent(out), optional :: grid(:,:)
    real, intent(out), optional :: lat_deriv(:,:)
    real, intent(out), optional :: lon_deriv(:,:)

    complex :: buff(1,size(spherical,1),size(spherical,2))
    real :: buff1(1,jec-jsc+1,iec-isc+1)
    real :: buff2(1,jec-jsc+1,iec-isc+1)

    buff(1,:,:) = spherical

    if (present(lat_deriv)) then
        call spherical_to_grid3D(buff,lat_deriv=buff1)
        lat_deriv=buff1(1,:,:)
    endif

    if (present(grid).and.present(lon_deriv)) then
        call spherical_to_grid3D(buff,grid=buff1,lon_deriv=buff2)
        grid = buff1(1,:,:)
        lon_deriv = buff2(1,:,:)
    elseif(present(grid).and..not.present(lon_deriv)) then
        call spherical_to_grid3D(buff,grid=buff1)
        grid = buff1(1,:,:)
    elseif(.not.present(grid).and.present(lon_deriv)) then
        call spherical_to_grid3D(buff,lon_deriv=buff2)
        lon_deriv = buff2(1,:,:)
    endif

    return
end subroutine spherical_to_grid2D


!--------------------------------------------------------------------------------
subroutine read_specdata(filename,fieldname,dat)
!--------------------------------------------------------------------------------   
    implicit none
    character (len=*), intent(in) :: filename, fieldname
    complex,dimension(:,:,:) :: dat

    real, allocatable :: rebuff(:,:), robuff(:,:)
    real, allocatable :: iebuff(:,:), iobuff(:,:)
    integer :: i, j, k, m
    character(len=len(fieldname)+3) :: iew, rew, iow, row 
    integer :: wdom(size(dat,2),2), neven_g, nodd_g

    iew = 'iew'//trim(fieldname)
    rew = 'rew'//trim(fieldname)
    iow = 'iow'//trim(fieldname)
    row = 'row'//trim(fieldname)

    call get_wdecomp(wdom,neven_g,nodd_g)

    allocate(rebuff(neven_g,size(dat,1)))
    allocate(iebuff(neven_g,size(dat,1)))
    allocate(robuff(nodd_g,size(dat,1)))
    allocate(iobuff(nodd_g,size(dat,1)))

    call read_data(filename,rew,rebuff)
    call read_data(filename,iew,iebuff)
    call read_data(filename,row,robuff)
    call read_data(filename,iow,iobuff)

    dat = cmplx(0.,0.)
    do k = 1, size(dat,1)
        do i = 1, size(dat,2)
            if (wdom(i,ev)>0) dat(k,i,ev) = cmplx(rebuff(wdom(i,ev),k),iebuff(wdom(i,ev),k))
            if (wdom(i,od)>0) dat(k,i,od) = cmplx(robuff(wdom(i,od),k),iobuff(wdom(i,od),k))
        enddo
    enddo

end subroutine read_specdata

end module transforms_mod

