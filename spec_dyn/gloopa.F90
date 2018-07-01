
#ifdef gloopa

program main

    use, intrinsic :: iso_c_binding

    use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
    use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
    use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
    use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
    use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

    use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

    use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init

    use fms_io_mod, only : fms_io_exit 

    use constants_mod, only : RADIUS

    use grid_fourier_mod, only : init_grid_fourier, fft_1dr2c_serial, fft_1dc2c_serial
    use grid_fourier_mod, only : end_grid_fourier, grid_to_fourier, fourier_to_grid

    use fourier_spherical_mod, only : init_fourier_spherical, fourier_to_spherical
    use fourier_spherical_mod, only : spherical_to_fourier

    use spherical_mod, only : compute_lon_deriv_cos, compute_lat_deriv_cos
    use spherical_mod, only : compute_vor_div, compute_ucos_vcos, cosm2_lat, Pnm
    use spherical_mod, only : triangle_mask, cosm_lat, sin_lat, ev, od
    use spherical_mod, only : do_truncation

    use vertical_levels_mod, only: init_vertical_levels, ak, bk

    use implicit_mod, only : init_implicit, do_implicit

    use gfidi_mod, only : gfidi_drv

    implicit none

    type satm_type
        complex, dimension(:,:,:),   allocatable :: vor
        complex, dimension(:,:,:),   allocatable :: div
        complex, dimension(:,:,:),   allocatable :: tem
        complex, dimension(:,:,:,:), allocatable :: tr
        complex, dimension(:,:,:),   allocatable :: prs
    end type satm_type

    type(satm_type) :: satm(3)

    complex, dimension(:,:,:),   allocatable :: sucos, svcos
    complex, dimension(:,:,:),   allocatable :: stopo

    real, allocatable, dimension(:,:,:) :: p, dpdphi, dpdlam, dpdt
    real, allocatable, dimension(:,:,:) :: u, dudphi, dudlam, dudt
    real, allocatable, dimension(:,:,:) :: v, dvdphi, dvdlam, dvdt
    real, allocatable, dimension(:,:,:) :: tem, dtemdphi, dtemdlam, dtemdt
    real, allocatable, dimension(:,:,:,:) :: tr, dtrdphi, dtrdlam, dtrdt
    real, allocatable, dimension(:,:,:) :: div, vor

    real, allocatable, dimension(:) :: spdmax

    type(domain2d) :: domainl
    type(domain2d) :: domainf

    logical :: fpe 
    integer :: nlon, nlat, nlev
      
    integer :: ilen, istart, olen, ostart, nlonb2, jlen
    integer :: nwaves_oe=0
        
    integer, allocatable :: pelist(:), extent(:)
    integer, allocatable :: fpelist(:), fextent(:)
    integer, allocatable :: Tshuff(:)
 
    character(len=32) :: routine='gloopa'

    integer :: comm, idfft3d, n, nt=1
    logical :: check=.false.
    real, allocatable :: tmp(:,:,:), tmp2d(:,:)
    complex, allocatable :: fldc1d(:,:), fldct(:,:,:)
    real :: ref_temp
    integer :: isc, iec, isg, ieg, m, l, t, i, ig, k, kstart=0, kend=0, kstep=1
    integer :: isf, ief, flen, jsc, jec, j, jsf, jef
    real :: scl, x, y, imgf=0.3, phi=0.15
    integer :: clck_grid2fourier, clck_fourier2grid, init, unit, cl, ck, num_fourier
    complex(kind=4) :: cpout(3)
    logical :: ideal_data=.false., debug=.false.
    real, parameter :: PI=4.D0*DATAN(1.D0)
    complex, parameter :: ui = cmplx(0.,1.), mui = -1.*ui
    integer :: pe, ntrac=3, ntr
    real :: deltim=600.
    character(len=16) :: rfile='gloopa', wfile='rgloopa'
    character(len=8) :: fldnm


    interface vor_div_to_uv_grid
        procedure vor_div_to_uv_grid3d
        procedure vor_div_to_uv_grid2d
    end interface vor_div_to_uv_grid

    interface vor_div_from_uv_grid
        procedure vor_div_from_uv_grid3D
        procedure vor_div_from_uv_grid2D
    end interface vor_div_from_uv_grid

    interface uv_grid_to_vor_div
        procedure uv_grid_to_vor_div3D
        procedure uv_grid_to_vor_div2D
    end interface uv_grid_to_vor_div

    interface read_griddataGFS
        procedure read_griddataGFS3D
        procedure read_griddataGFS2D
    end interface read_griddataGFS

    interface write_griddata
        procedure write_griddata3D
        procedure write_griddata2D
    end interface write_griddata

    interface grid_to_spherical
        procedure grid_to_spherical3D
        procedure grid_to_spherical2D
    end interface

    namelist/gloopa_nml/kstart, kend, kstep, ideal_data, imgf, ck, cl, &
                        check, num_fourier, nt, nlon, nlat, nlev, debug, deltim

    call mpp_init() 
    call fms_init()

    pe = mpp_pe()

    unit = open_namelist_file()
    read(unit,nml=gloopa_nml)
    call close_file(unit)

    clck_grid2fourier = mpp_clock_id('grid_to_fourier')
    clck_fourier2grid = mpp_clock_id('fourier_to_grid')

    allocate(pelist(mpp_npes()))
    allocate(extent(mpp_npes()))

    call mpp_define_domains( [1,nlat,1,nlon], [1,mpp_npes()], domainl, kxy=1)
    call mpp_get_compute_domain(domainl, jsc, jec, isc, iec)
    ilen = iec-isc+1
    jlen = jec-jsc+1

    call mpp_get_current_pelist(pelist,commid=comm)

    allocate(Tshuff(0:num_fourier))
    call init_grid_fourier (nlon, ilen, num_fourier, isf, flen, comm, Tshuff)

    call mpp_gather([flen], extent)
    call mpp_broadcast(extent,size(extent), mpp_root_pe())

    allocate(fextent(count(extent>0)))
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
        call mpp_define_domains( [1,nlat, 0,num_fourier], [1,mpp_npes()], domainf, yextent=fextent)
        call mpp_get_compute_domain(domainf, jsf, jef, isf, ief)
        if(flen /= ief-isf+1) call mpp_error('gloopa', 'flen /= ief-isf+1', FATAL)
        print *, 'pe, isf, flen, load=', mpp_pe(), isf, flen, sum(num_fourier-Tshuff(isf:ief)+1)

        call init_fourier_spherical(num_fourier, num_fourier, nlat, nwaves_oe, domainf, Tshuff) 
    endif
    call mpp_set_current_pelist()

    call init_vertical_levels(nlev)

    ref_temp = 300.
    if (nlev>100) ref_temp=1500.

    call init_implicit(ak,bk,ref_temp,deltim,num_fourier)

    allocate(sucos(nlev,nwaves_oe,2))
    allocate(svcos(nlev,nwaves_oe,2))

    allocate(satm(1)%vor(nlev,nwaves_oe,2))
    allocate(satm(1)%div(nlev,nwaves_oe,2))
    allocate(satm(1)%tem(nlev,nwaves_oe,2))
    allocate(satm(1)%tr(nlev,nwaves_oe,2,ntrac))
    allocate(satm(1)%prs(1,nwaves_oe,2))

    allocate(satm(2)%vor(nlev,nwaves_oe,2))
    allocate(satm(2)%div(nlev,nwaves_oe,2))
    allocate(satm(2)%tem(nlev,nwaves_oe,2))
    allocate(satm(2)%tr(nlev,nwaves_oe,2,ntrac))
    allocate(satm(2)%prs(1,nwaves_oe,2))

    allocate(satm(3)%vor(nlev,nwaves_oe,2))
    allocate(satm(3)%div(nlev,nwaves_oe,2))
    allocate(satm(3)%tem(nlev,nwaves_oe,2))
    allocate(satm(3)%tr(nlev,nwaves_oe,2,ntrac))
    allocate(satm(3)%prs(1,nwaves_oe,2))

    allocate(stopo(1,nwaves_oe,2))

    allocate(u(nlev,jsc:jec,isc:iec))
    allocate(dudlam(nlev,jsc:jec,isc:iec))
    allocate(dudphi(nlev,jsc:jec,isc:iec))
    allocate(dudt(nlev,jsc:jec,isc:iec))

    allocate(p(1,jsc:jec,isc:iec))
    allocate(dpdlam(1,jsc:jec,isc:iec))
    allocate(dpdphi(1,jsc:jec,isc:iec))
    allocate(dpdt(1,jsc:jec,isc:iec))

    allocate(v(nlev,jsc:jec,isc:iec))
    allocate(dvdlam(nlev,jsc:jec,isc:iec))
    allocate(dvdphi(nlev,jsc:jec,isc:iec))
    allocate(dvdt(nlev,jsc:jec,isc:iec))

    allocate(tem(nlev,jsc:jec,isc:iec))
    allocate(dtemdlam(nlev,jsc:jec,isc:iec))
    allocate(dtemdphi(nlev,jsc:jec,isc:iec))
    allocate(dtemdt(nlev,jsc:jec,isc:iec))

    allocate(tr(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dtrdlam(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dtrdphi(nlev,jsc:jec,isc:iec,ntrac))
    allocate(dtrdt(nlev,jsc:jec,isc:iec,ntrac))

    allocate(div(nlev,jsc:jec,isc:iec))
    allocate(vor(nlev,jsc:jec,isc:iec))
    allocate(tmp2d(jsc:jec,isc:iec))
    allocate(spdmax(nlev))
    spdmax = 0.

    call read_specdatagfs('specdata','topo',stopo)

    call read_specdatagfs('specdata','lnp_1',satm(1)%prs)
    call read_specdatagfs('specdata','lnp_2',satm(2)%prs)

    call read_specdatagfs('specdata','vor_1',satm(1)%vor)
    call read_specdatagfs('specdata','vor_2',satm(2)%vor)

    call read_specdatagfs('specdata','div_1',satm(1)%div)
    call read_specdatagfs('specdata','div_2',satm(2)%div)

    call read_specdatagfs('specdata','tem_1',satm(1)%tem)
    call read_specdatagfs('specdata','tem_2',satm(2)%tem)

    do ntr = 1, ntrac
        write(fldnm,'(A,I1)') 'tr',ntr
        call read_specdatagfs('specdata',trim(fldnm)//'_1',satm(1)%tr(:,:,:,ntr))
        call read_specdatagfs('specdata',trim(fldnm)//'_2',satm(2)%tr(:,:,:,ntr))
    enddo

    call compute_ucos_vcos(satm(2)%vor,satm(2)%div,sucos,svcos,do_trunc=.false.)

    call spherical_to_grid3d(satm(2)%div,grid=div)

    call spherical_to_grid3d(satm(2)%vor,grid=vor)

    call spherical_to_grid3d(sucos,grid=u,lon_deriv=dudlam)

    call spherical_to_grid3d(svcos,grid=v,lon_deriv=dvdlam)

    call spherical_to_grid3d(satm(2)%prs,grid=p,lat_deriv=dpdphi,lon_deriv=dpdlam)

    call spherical_to_grid3d(satm(2)%tem,grid=tem,lat_deriv=dtemdphi,lon_deriv=dtemdlam)

    do ntr = 1, ntrac
        call spherical_to_grid3d(satm(2)%tr(:,:,:,ntr),grid=tr(:,:,:,ntr), &
            lat_deriv=dtrdphi(:,:,:,ntr),lon_deriv=dtrdlam(:,:,:,ntr))
    enddo

    do j = jsc, jec
        dtemdphi(:,j,:) = dtemdphi(:,j,:) * cosm2_lat(j)
        dtrdphi(:,j,:,:) = dtrdphi(:,j,:,:) * cosm2_lat(j)

        dtemdlam(:,j,:) = dtemdlam(:,j,:) * cosm2_lat(j)
        dtrdlam(:,j,:,:) = dtrdlam(:,j,:,:) * cosm2_lat(j)

        dudlam(:,j,:) = dudlam(:,j,:) * cosm2_lat(j)
        dvdlam(:,j,:) = dvdlam(:,j,:) * cosm2_lat(j)
    enddo

    dudphi = dvdlam - vor
    dvdphi = div - dudlam

    call gfidi_drv(nlev, ntrac, ilen, jlen, deltim, sin_lat(jsc:jec), cosm2_lat(jsc:jec), &
            div, tem, u, v, tr, dpdphi, dpdlam, p, dtemdphi, dtemdlam, dtrdphi, &
            dtrdlam, dudlam, dvdlam, dudphi, dvdphi, dpdt, dtemdt, dtrdt, dudt, dvdt, spdmax)

    call write_griddata('rgloopa','div',div)
    call write_griddata('rgloopa','vor',vor)

    call write_griddata('rgloopa','dudt',dudt)
    call write_griddata('rgloopa','dudphi',dudphi)
    call write_griddata('rgloopa','dudlam',dudlam)
    call write_griddata('rgloopa','u',u)

    call write_griddata('rgloopa','dvdt',dvdt)
    call write_griddata('rgloopa','dvdphi',dvdphi)
    call write_griddata('rgloopa','dvdlam',dvdlam)
    call write_griddata('rgloopa','v',v)

    call write_griddata('rgloopa','dpdt',dpdt)
    call write_griddata('rgloopa','dpdphi',dpdphi)
    call write_griddata('rgloopa','dpdlam',dpdlam)
    call write_griddata('rgloopa','p',p)

    call write_griddata('rgloopa','dtemdt',dtemdt)
    call write_griddata('rgloopa','dtemdphi',dtemdphi)
    call write_griddata('rgloopa','dtemdlam',dtemdlam)
    call write_griddata('rgloopa','tem',tem)

    do ntr = 1, ntrac
        write(fldnm,'(A,I2.2)') 'tr',ntr
        print *, trim(fldnm)
        call write_griddata('rgloopa','d'//trim(fldnm)//'dt',dtrdt(:,:,:,ntr))
        call write_griddata('rgloopa','d'//trim(fldnm)//'dphi',dtrdphi(:,:,:,ntr))
        call write_griddata('rgloopa','d'//trim(fldnm)//'dlam',dtrdlam(:,:,:,ntr))
        call write_griddata('rgloopa',trim(fldnm),tr(:,:,:,ntr))
    enddo

    call grid_to_spherical(dpdt, satm(3)%prs, do_trunc=.true.)
    call grid_to_spherical(dtemdt, satm(3)%tem, do_trunc=.true.)
    call grid_to_spherical(dudt, sucos, do_trunc=.false.)
    call grid_to_spherical(dvdt, svcos, do_trunc=.false.)

    do ntr = 1, ntrac
        call grid_to_spherical(dtrdt(:,:,:,ntr),satm(3)%tr(:,:,:,ntr),do_trunc=.true.)
    enddo

    call compute_vor_div(sucos,svcos,satm(3)%vor,satm(3)%div)

    call spherical_to_grid3D(satm(3)%div,grid=div)

    call write_griddata('rgloopa','divdt2',div)

    satm(3)%vor = satm(1)%vor + 2.*deltim*satm(3)%vor
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr
    do k = 1, size(satm(3)%div,1)
        satm(3)%div(k,:,:) = satm(3)%div(k,:,:) + stopo(1,:,:)
    enddo
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr

    !call read_specdatagfs('specdata','lnp_1',satm(1)%prs)
    !call read_specdatagfs('specdata','lnp_2',satm(2)%prs)
    !call read_specdatagfs('specdata','lnp_3',satm(3)%prs)

    !call read_specdatagfs('specdata','div_1',satm(1)%div)
    !call read_specdatagfs('specdata','div_2',satm(2)%div)
    !call read_specdatagfs('specdata','div_3',satm(3)%div)

    !call read_specdatagfs('specdata','tem_1',satm(1)%tem)
    !call read_specdatagfs('specdata','tem_2',satm(2)%tem)
    !call read_specdatagfs('specdata','tem_3',satm(3)%tem)

    call spherical_to_grid3D(satm(3)%div,grid=div)
    call spherical_to_grid3D(satm(3)%tem,grid=tem)
    call spherical_to_grid3D(satm(3)%prs,grid=p)

    call write_griddata('rgloopa','divdt1',div)
    call write_griddata('rgloopa','temdt1',tem)
    call write_griddata('rgloopa','prsdt1',p)

    call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)

    call write_specdata('rgloopa','divdt',satm(3)%div)
    call write_specdata('rgloopa','temdt',satm(3)%tem)
    call write_specdata('rgloopa','prsdt',satm(3)%prs)

    call spherical_to_grid3D(satm(3)%div,grid=div)
    call spherical_to_grid3D(satm(3)%tem,grid=tem)
    call spherical_to_grid3D(satm(3)%prs,grid=p)

    call write_griddata('rgloopa','divdt',div)
    call write_griddata('rgloopa','temdt',tem)
    call write_griddata('rgloopa','prsdt',p)

    call fms_io_exit()
    call end_grid_fourier()
    call mpp_exit()

    contains

subroutine vor_div_to_uv_grid3d(vor,div,u,v,getcosuv)

    complex,dimension(:,:,:), intent(in) :: vor, div
    real, intent(out) :: u(:,:,:), v(:,:,:)
    logical, intent(in), optional :: getcosuv
    complex,dimension(size(vor,1),size(vor,2),2) :: sucos, svcos
    logical :: getcosuv1

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

    call vor_div_to_uv_grid3d(vor, div, u3d, v3d,getcosuv)

    u(:,:) = u3d(1,:,:)
    v(:,:) = v3d(1,:,:)

    return

end subroutine vor_div_to_uv_grid2d


subroutine vor_div_from_uv_grid3D(u,v,vor,div,uvcos_in)

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


subroutine uv_grid_to_vor_div3D(u,v,vor,div)
    real, intent(in) :: u(:,:,:), v(:,:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div

    complex,dimension(size(vor,1),size(vor,2),2) :: us, vs

    complex :: four(size(u,1)*size(u,2), flen)
    complex, pointer :: four3(:,:,:)

    real, pointer :: grd(:,:)
    type(C_PTR) :: pgrd, pfour

    integer :: nk, nj, ni, howmany
    integer :: i, j, k

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

subroutine vor_div_from_uv_grid2D(u,v,vor,div)
    real, intent(in) :: u(:,:), v(:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div
    real :: u3d(1,size(u,1),size(u,2))
    real :: v3d(1,size(u,1),size(u,2))

    u3d(1,:,:) = u(:,:)
    v3d(1,:,:) = v(:,:)
    call vor_div_from_uv_grid3D(u3d,v3d,vor,div)

end subroutine vor_div_from_uv_grid2D

subroutine uv_grid_to_vor_div2D(u,v,vor,div)
    real, intent(in) :: u(:,:), v(:,:)
    complex,dimension(:,:,:), intent(out) :: vor, div
    real :: u3d(1,size(u,1),size(u,2))
    real :: v3d(1,size(u,1),size(u,2))

    u3d(1,:,:) = u(:,:)
    v3d(1,:,:) = v(:,:)
    call uv_grid_to_vor_div3D(u3d,v3d,vor,div)

end subroutine uv_grid_to_vor_div2D


subroutine grid_to_spherical3D(grid,spherical,do_trunc)
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


subroutine grid_to_spherical2D(grid,spherical,do_trunc)
    real, intent(in) :: grid(:,:)
    complex, dimension(:,:,:), intent(out) :: spherical
    logical, intent(in), optional :: do_trunc

    real :: buff(1,size(grid,1),size(grid,2))

    buff(1,:,:) = grid

    call grid_to_spherical3D(buff, spherical, do_trunc)

    return
end subroutine grid_to_spherical2D


subroutine spherical_to_grid2D(spherical,grid,lat_deriv,lon_deriv)
    complex,dimension(:,:,:), intent(in) :: spherical
    real, intent(out), optional :: grid(:,:)
    real, intent(out), optional :: lat_deriv(:,:)
    real, intent(out), optional :: lon_deriv(:,:)

    real :: buff1(1,jec-jsc+1,iec-isc+1)
    real :: buff2(1,jec-jsc+1,iec-isc+1)

    if (present(lat_deriv)) then
        call spherical_to_grid3D(spherical,lat_deriv=buff1)
        lat_deriv=buff1(1,:,:)
    endif

    if (present(grid).and.present(lon_deriv)) then
        call spherical_to_grid3D(spherical,grid=buff1,lon_deriv=buff2)
        grid = buff1(1,:,:)
        lon_deriv = buff2(1,:,:)
    elseif(present(grid).and..not.present(lon_deriv)) then
        call spherical_to_grid3D(spherical,grid=buff1)
        grid = buff1(1,:,:)
    elseif(.not.present(grid).and.present(lon_deriv)) then
        call spherical_to_grid3D(spherical,lon_deriv=buff2)
        lon_deriv = buff2(1,:,:)
    endif

    return
end subroutine spherical_to_grid2D



subroutine write_specdata(filename,fieldname,dat)

    character (len=*), intent(in) :: filename, fieldname
    complex, dimension(:,:,:), intent(in) :: dat

    real :: robuff1(size(dat,2)-num_fourier/2,size(dat,1))
    real :: iobuff1(size(dat,2)-num_fourier/2,size(dat,1))
    real :: rebuff(size(dat,2),size(dat,1)), robuff(size(dat,2),size(dat,1))
    real :: iebuff(size(dat,2),size(dat,1)), iobuff(size(dat,2),size(dat,1))
    integer :: i, j, k, nlev, m, nodd
    character(len=len(fieldname)+3) :: iew, rew, iow, row 
    integer :: nlen1, nlen2, ns1, ns2, ne1, ne2
    logical :: gfs_type

    gfs_type=.true.

    iew = 'iew'//trim(fieldname)
    rew = 'rew'//trim(fieldname)
    iow = 'iow'//trim(fieldname)
    row = 'row'//trim(fieldname)

    nlev = size(dat,1)
    nodd = size(dat,2)-num_fourier/2

    robuff = 0.
    iobuff = 0.

    do k = 1, nlev
        rebuff(:,k) = real(dat(k,:,ev))
        iebuff(:,k) = aimag(dat(k,:,ev))
        robuff(:,k) = real(dat(k,:,od))
        iobuff(:,k) = aimag(dat(k,:,od))
        nlen1 = 32; nlen2 = 33
        ne1 = 0; ne2 = 0
        do m = isf, ief
            if (mod(m+1,2)==0) nlen1 = nlen1 - 1
            if (mod(m,2)==0) nlen2 = nlen2 - 1
            ns1 = ne1 + 1; ns2 = ne2 + 1
            ne1 = ns1 + nlen1 - 1; ne2 = ns2 + nlen2 -1
            robuff1(ns1:ne1,k) = real(dat(k,ns2:ns2+nlen1-1,od))
            iobuff1(ns1:ne1,k) = aimag(dat(k,ns2:ns2+nlen1-1,od))
        enddo
    enddo

    call write_data(filename,rew,rebuff)
    call write_data(filename,iew,iebuff)
    if (gfs_type) then
    call write_data(filename,row,robuff1)
    call write_data(filename,iow,iobuff1)
    else
    call write_data(filename,row,robuff)
    call write_data(filename,iow,iobuff)
    endif
end subroutine write_specdata

subroutine read_specdataGFS(filename,fieldname,dat)

    character (len=*), intent(in) :: filename, fieldname
    complex,dimension(:,:,:) :: dat

    real :: rebuff(size(dat,2),size(dat,1)), robuff(size(dat,2)-num_fourier/2,size(dat,1))
    real :: iebuff(size(dat,2),size(dat,1)), iobuff(size(dat,2)-num_fourier/2,size(dat,1))
    integer :: i, j, k, nlev, m, nodd
    character(len=len(fieldname)+3) :: iew, rew, iow, row 
    integer :: nlen1, nlen2, ns1, ns2, ne1, ne2

    iew = 'iew'//trim(fieldname)
    rew = 'rew'//trim(fieldname)
    iow = 'iow'//trim(fieldname)
    row = 'row'//trim(fieldname)

    nlev = size(dat,1)
    nodd = size(dat,2)-num_fourier/2
    call read_data(filename,rew,rebuff)
    call read_data(filename,iew,iebuff)
    call read_data(filename,row,robuff)
    call read_data(filename,iow,iobuff)

    dat = cmplx(0.,0.)
    do k = 1, nlev
        dat(k,:,ev) = cmplx(rebuff(:,k),iebuff(:,k))
        nlen1 = 32; nlen2 = 33
        ne1 = 0; ne2 = 0
        do m = isf, ief
            if (mod(m+1,2)==0) nlen1 = nlen1 - 1
            if (mod(m,2)==0) nlen2 = nlen2 - 1
            ns1 = ne1 + 1; ns2 = ne2 + 1
            ne1 = ns1 + nlen1 - 1; ne2 = ns2 + nlen2 -1
            dat(k,ns2:ns2+nlen1-1,od) = cmplx(robuff(ns1:ne1,k),iobuff(ns1:ne1,k))
        enddo
    enddo

end subroutine read_specdataGFS
    
subroutine read_griddataGFS3D(filename,fieldname,dat)
    character (len=*), intent(in) :: filename, fieldname
    real, intent(out) :: dat(:,:,:)
   
    real :: buff(size(dat,3),size(dat,2),size(dat,1))
    integer :: i, j, k, j_s, j_n, j_sb, j_nb

    call read_data(filename,fieldname,buff)

    do k = 1, size(dat,1)
        do j = 1, size(dat,2)/2
            j_s = 2*(j-1)+1; j_n = 2*j
            j_nb = j; j_sb = size(dat,2) - j + 1
            do i = 1, size(dat,3)
                dat(k,j_s,i) = buff(i,j_sb,k)
                dat(k,j_n,i) = buff(i,j_nb,k)
            enddo
        enddo
    enddo
    return
end subroutine read_griddataGFS3D

        
subroutine read_griddataGFS2D(filename,fieldname,dat)
    character (len=*), intent(in) :: filename, fieldname
    real, intent(out) :: dat(:,:)
   
    real :: buff(1,size(dat,1),size(dat,2))
    integer :: i, j, k, j_s, j_n, j_sb, j_nb

    call read_griddataGFS3D(filename,fieldname,buff)

    dat(:,:) = buff(1,:,:)

    return
end subroutine read_griddataGFS2D

subroutine write_griddata3D(filename,fieldname,dat)
    character (len=*), intent(in) :: filename, fieldname
    real, intent(in) :: dat(:,:,:)
   
    real :: buff(size(dat,3),size(dat,2),size(dat,1))
    integer :: i, j, k, j_s, j_n, j_sb, j_nb

    do k = 1, size(dat,1)
        do j = 1, size(dat,2)/2
            j_s = 2*(j-1)+1; j_n = 2*j
            j_nb = j; j_sb = size(dat,2) - j + 1
            do i = 1, size(dat,3)
                buff(i,j_sb,k) = dat(k,j_s,i)
                buff(i,j_nb,k) = dat(k,j_n,i)
            enddo
        enddo
    enddo

    call write_data(filename,fieldname,buff)

    return
end subroutine write_griddata3D


subroutine write_griddata2D(filename,fieldname,dat)
    character (len=*), intent(in) :: filename, fieldname
    real, intent(in) :: dat(:,:)
   
    real :: buff(1,size(dat,1),size(dat,2))
    integer :: i, j, k, j_s, j_n, j_sb, j_nb


    buff(1,:,:) = dat(:,:)
    call write_griddata3D(filename,fieldname,buff)

    return
end subroutine write_griddata2D

end program main

#endif
