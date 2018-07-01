
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

    use transforms_mod, only : read_specdatagfs, write_specdata, write_griddata
    use transforms_mod, only : compute_ucos_vcos, compute_vor_div
    use transforms_mod, only : spherical_to_grid, grid_to_spherical
    use transforms_mod, only : cosm2_lat, sin_lat, init_transforms

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

    integer :: nlon, nlat, nlev
      
    integer :: nwaves_oe=0
        
    integer :: isc, iec, ilen, i, k
    integer :: jsc, jec, jlen, j
    integer :: unit, trunc
    integer :: pe, ntrac=3, ntr
    real :: deltim=600., ref_temp
    character(len=16) :: rfile='gloopa', wfile='rgloopa'
    character(len=8) :: fldnm

    namelist/gloopa_nml/ trunc, nlon, nlat, nlev, deltim

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=gloopa_nml)
    call close_file(unit)

    call init_transforms(nlon,nlat,trunc,isc,iec,jsc,jec,nwaves_oe)
    jlen = jec-jsc+1
    ilen = iec-isc+1

    call init_vertical_levels(nlev)

    ref_temp = 300.
    if (nlev>100) ref_temp=1500.

    call init_implicit(ak,bk,ref_temp,deltim,trunc)

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

    call spherical_to_grid(satm(2)%div,grid=div)

    call spherical_to_grid(satm(2)%vor,grid=vor)

    call spherical_to_grid(sucos,grid=u,lon_deriv=dudlam)

    call spherical_to_grid(svcos,grid=v,lon_deriv=dvdlam)

    call spherical_to_grid(satm(2)%prs,grid=p,lat_deriv=dpdphi,lon_deriv=dpdlam)

    call spherical_to_grid(satm(2)%tem,grid=tem,lat_deriv=dtemdphi,lon_deriv=dtemdlam)

    do ntr = 1, ntrac
        call spherical_to_grid(satm(2)%tr(:,:,:,ntr),grid=tr(:,:,:,ntr), &
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

    call spherical_to_grid(satm(3)%div,grid=div)

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

    call spherical_to_grid(satm(3)%div,grid=div)
    call spherical_to_grid(satm(3)%tem,grid=tem)
    call spherical_to_grid(satm(3)%prs,grid=p)

    call write_griddata('rgloopa','divdt1',div)
    call write_griddata('rgloopa','temdt1',tem)
    call write_griddata('rgloopa','prsdt1',p)

    call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)

    call write_specdata('rgloopa','divdt',satm(3)%div)
    call write_specdata('rgloopa','temdt',satm(3)%tem)
    call write_specdata('rgloopa','prsdt',satm(3)%prs)

    call spherical_to_grid(satm(3)%div,grid=div)
    call spherical_to_grid(satm(3)%tem,grid=tem)
    call spherical_to_grid(satm(3)%prs,grid=p)

    call write_griddata('rgloopa','divdt',div)
    call write_griddata('rgloopa','temdt',tem)
    call write_griddata('rgloopa','prsdt',p)

    call fms_io_exit()
    call mpp_exit()


end program main
