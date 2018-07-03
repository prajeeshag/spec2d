
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

    use transforms_mod, only : read_specdata, get_spherical_wave
    use transforms_mod, only : compute_ucos_vcos, compute_vor_div
    use transforms_mod, only : spherical_to_grid, grid_to_spherical
    use transforms_mod, only : cosm2_lat, sin_lat, init_transforms

    use vertical_levels_mod, only: init_vertical_levels, get_ak_bk

    use implicit_mod, only : init_implicit, do_implicit

    use gfidi_mod, only : gfidi_drv
    
    use horiz_diffusion_mod, only : init_horiz_diffusion, horiz_diffusion

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
    type(domain2d) :: domain_g

    real, allocatable :: ak(:), bk(:), sl(:)
    integer, allocatable :: sph_wave(:,:)

    namelist/gloopa_nml/ trunc, nlon, nlat, nlev, deltim

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=gloopa_nml)
    call close_file(unit)

    call mpp_define_domains( [1,nlat,1,nlon], [1,mpp_npes()], domain_g, kxy=1, ishuff=1)
    call mpp_get_compute_domain(domain_g, jsc, jec, isc, iec)
    ilen = iec-isc+1
    jlen = jec-jsc+1

    call init_transforms(domain_g,trunc,nwaves_oe)
    jlen = jec-jsc+1
    ilen = iec-isc+1

    call init_vertical_levels(nlev)

    allocate(ak(nlev+1),bk(nlev+1),sl(nlev))

    allocate(sph_wave(nwaves_oe,2))

    call get_ak_bk(ak_out=ak,bk_out=bk,sl_out=sl)
    call get_spherical_wave(sph_wave)

    ref_temp = 300.
    if (nlev>100) ref_temp=1500.

    call init_implicit(ak,bk,ref_temp,deltim,trunc)

    call init_horiz_diffusion(trunc,deltim,sl,sph_wave,bk)

    allocate(sucos(nlev,nwaves_oe,2))
    allocate(svcos(nlev,nwaves_oe,2))

    do i = 1, 3
        allocate(satm(i)%vor(nlev,nwaves_oe,2))
        allocate(satm(i)%div(nlev,nwaves_oe,2))
        allocate(satm(i)%tem(nlev,nwaves_oe,2))
        allocate(satm(i)%tr(nlev,nwaves_oe,2,ntrac))
        allocate(satm(i)%prs(1,nwaves_oe,2))
    enddo

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

    call read_specdata('specdata','topo',stopo)

    call read_specdata('specdata','lnp_1',satm(1)%prs)
    call read_specdata('specdata','lnp_2',satm(2)%prs)

    call read_specdata('specdata','vor_1',satm(1)%vor)
    call read_specdata('specdata','vor_2',satm(2)%vor)

    call read_specdata('specdata','div_1',satm(1)%div)
    call read_specdata('specdata','div_2',satm(2)%div)

    call read_specdata('specdata','tem_1',satm(1)%tem)
    call read_specdata('specdata','tem_2',satm(2)%tem)

    do ntr = 1, ntrac
        write(fldnm,'(A,I1)') 'tr',ntr
        call read_specdata('specdata',trim(fldnm)//'_1',satm(1)%tr(:,:,:,ntr))
        call read_specdata('specdata',trim(fldnm)//'_2',satm(2)%tr(:,:,:,ntr))
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

    call write_data('rgloopa','div',div,domain_g)
    call write_data('rgloopa','vor',vor,domain_g)

    call write_data('rgloopa','dudt',dudt,domain_g)
    call write_data('rgloopa','dudphi',dudphi,domain_g)
    call write_data('rgloopa','dudlam',dudlam,domain_g)
    call write_data('rgloopa','u',u,domain_g)

    call write_data('rgloopa','dvdt',dvdt,domain_g)
    call write_data('rgloopa','dvdphi',dvdphi,domain_g)
    call write_data('rgloopa','dvdlam',dvdlam,domain_g)
    call write_data('rgloopa','v',v,domain_g)

    call write_data('rgloopa','dpdt',dpdt,domain_g)
    call write_data('rgloopa','dpdphi',dpdphi,domain_g)
    call write_data('rgloopa','dpdlam',dpdlam,domain_g)
    call write_data('rgloopa','p',p,domain_g)

    call write_data('rgloopa','dtemdt',dtemdt,domain_g)
    call write_data('rgloopa','dtemdphi',dtemdphi,domain_g)
    call write_data('rgloopa','dtemdlam',dtemdlam,domain_g)
    call write_data('rgloopa','tem',tem,domain_g)

    do ntr = 1, ntrac
        write(fldnm,'(A,I2.2)') 'tr',ntr
        call write_data('rgloopa','d'//trim(fldnm)//'dt',dtrdt(:,:,:,ntr),domain_g)
        call write_data('rgloopa','d'//trim(fldnm)//'dphi',dtrdphi(:,:,:,ntr),domain_g)
        call write_data('rgloopa','d'//trim(fldnm)//'dlam',dtrdlam(:,:,:,ntr),domain_g)
        call write_data('rgloopa',trim(fldnm),tr(:,:,:,ntr),domain_g)
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

    call write_data('rgloopa','divdt2',div,domain_g)

    satm(3)%vor = satm(1)%vor + 2.*deltim*satm(3)%vor
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr
    do k = 1, size(satm(3)%div,1)
        satm(3)%div(k,:,:) = satm(3)%div(k,:,:) + stopo(1,:,:)
    enddo
    satm(3)%tr = satm(1)%tr + 2.*deltim*satm(3)%tr

    call spherical_to_grid(satm(3)%div,grid=div)
    call spherical_to_grid(satm(3)%tem,grid=tem)
    call spherical_to_grid(satm(3)%prs,grid=p)

    call write_data('rgloopa','divdt1',div,domain_g)
    call write_data('rgloopa','temdt1',tem,domain_g)
    call write_data('rgloopa','prsdt1',p,domain_g)

    call do_implicit(satm(1)%div, satm(1)%tem, satm(1)%prs, &
                     satm(2)%div, satm(2)%tem, satm(2)%prs, &
                     satm(3)%div, satm(3)%tem, satm(3)%prs, deltim)

    call horiz_diffusion(satm(3)%tr,satm(3)%vor,satm(3)%div,satm(3)%tem,satm(1)%prs(1,:,:))

    call spherical_to_grid(satm(3)%div,grid=div)
    call spherical_to_grid(satm(3)%tem,grid=tem)
    call spherical_to_grid(satm(3)%prs,grid=p)

    call write_data('rgloopa','divdt',div,domain_g)
    call write_data('rgloopa','temdt',tem,domain_g)
    call write_data('rgloopa','prsdt',p,domain_g)

    call fms_io_exit()
    call mpp_exit()

end program main
