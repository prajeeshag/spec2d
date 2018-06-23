
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

    use grid_fourier_mod, only : init_grid_fourier, fft_1dr2c_serial, fft_1dc2c_serial
    use grid_fourier_mod, only : end_grid_fourier, grid_to_fourier, fourier_to_grid

    use fourier_spherical_mod, only : init_fourier_spherical, fourier_to_spherical, spherical_to_fourier
    use spherical_mod, only : specVar, compute_lon_deriv_cos, compute_lat_deriv_cos
    use spherical_mod, only : compute_vor_div, compute_ucos_vcos, cos_lat, legendre

    implicit none

    type(domain2d) :: domainl
    type(domain2d) :: domainf

    logical :: fpe 
    integer :: nlon, nlat, nlev
      
    integer :: ilen, istart, olen, ostart, nlonb2
    integer :: nwaves_oe=0
        
    integer, allocatable :: pelist(:), extent(:)
    integer, allocatable :: fpelist(:), fextent(:)
    integer, allocatable :: Tshuff(:)
 
    character(len=32) :: routine='gloopa'

    integer :: comm, idfft3d, n, nt=1
    logical :: check=.false.
    real, allocatable :: lnp(:,:), dpdphi(:,:), dpdlam(:,:)
    real, allocatable :: tmp(:,:,:), tmp2d(:,:)
    complex, allocatable :: fldc1d(:,:), fldct(:,:,:)
    integer :: isc, iec, isg, ieg, m, l, t, i, ig, k, kstart=0, kend=0, kstep=1
    integer :: isf, ief, flen, jsc, jec, j, jsf, jef
    real :: scl, x, y, imgf=0.3, phi=0.15
    integer :: clck_grid2fourier, clck_fourier2grid, init, unit, cl, ck, num_fourier
    complex(kind=4) :: cpout(3)
    logical :: ideal_data=.false., debug=.false.
    real, parameter :: PI=4.D0*DATAN(1.D0)
    complex, parameter :: ui = cmplx(0.,1.), mui = -1.*ui
    integer :: pe

    type(specVar(n=:,nlev=:)), allocatable :: slnp, slnpdlam, slnpdphi
    type(specVar(n=:,nlev=:)), allocatable :: ucos, vcos, vor, div

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

    interface spherical_to_grid
        procedure spherical_to_grid3D
        procedure spherical_to_grid2D
    end interface

    namelist/gloopa_nml/kstart, kend, kstep, ideal_data, imgf, ck, cl, &
                                      check, num_fourier, nt, nlon, nlat, nlev, debug

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

    allocate(lnp(jsc:jec,isc:iec))
    allocate(tmp(nlev,jsc:jec,isc:iec))
    allocate(tmp2d(jsc:jec,isc:iec))
    allocate(dpdphi(jsc:jec,isc:iec))
    allocate(dpdlam(jsc:jec,isc:iec))

    allocate(specVar(n=nwaves_oe,nlev=1) :: slnp, slnpdlam, slnpdphi)
    allocate(specVar(n=nwaves_oe,nlev=nlev) :: vor, div, ucos, vcos)

    call read_specdataGFS('gloopa', 'vor', vor)

    call spherical_to_grid(vor,tmp)
    call write_griddata('rgloopa', 'vor', tmp)

    vor%ev = cmplx(0.,0.)
    vor%od = cmplx(0.,0.)
    call grid_to_spherical(tmp,vor)

    tmp = 0.
    call spherical_to_grid(vor,tmp)
    call write_griddata('rgloopa', 'vor2', tmp)

    
    call write_data('rgloopa', 'legev', legendre%ev)
    call write_data('rgloopa', 'legod', legendre%od)

    call fms_io_exit()
    call end_grid_fourier()
    call mpp_exit()


    contains


    subroutine grid_to_spherical3D(grid,spherical)
        real, intent(in) :: grid(:,:,:)
        type(specVar(n=*,nlev=*)), intent(out) :: spherical
    
        complex :: four(size(grid,1)*size(grid,2), flen)
        complex, pointer :: four3(:,:,:)

        real, pointer :: grd(:,:)
        type(C_PTR) :: pgrd, pfour

        integer :: nk, nj, ni, howmany
        integer :: i, j, k

        nk = size(grid,1); nj = size(grid,2); ni = size(grid,3)

        howmany = nk*nj

        pgrd = C_LOC(grid)
        call c_f_pointer(pgrd, grd, [howmany, ni])

        pfour = C_LOC(four)
        call c_f_pointer(pfour, four3, [nk, nj, flen])

        call grid_to_fourier(grd,four)

        call fourier_to_spherical(four3,spherical)

        return
    end subroutine grid_to_spherical3D 


    subroutine spherical_to_grid3D(spherical,grid,lat_deriv)
        real, intent(out) :: grid(:,:,:)
        type(specVar(n=*,nlev=*)), intent(in) :: spherical
        logical, intent(in), optional :: lat_deriv
    
        complex :: four(size(grid,1)*size(grid,2), flen)
        complex, pointer :: four3(:,:,:)

        real, pointer :: grd(:,:)
        type(C_PTR) :: pgrd, pfour
        logical :: lat_deriv1

        integer :: nk, nj, ni, howmany
        integer :: i, j, k
    
        lat_deriv1 = .false.

        if (present(lat_deriv)) lat_deriv1 = lat_deriv

        nk = size(grid,1); nj = size(grid,2); ni = size(grid,3)

        howmany = nk*nj

        pgrd = C_LOC(grid)
        call c_f_pointer(pgrd, grd, [howmany, ni])

        pfour = C_LOC(four)
        call c_f_pointer(pfour, four3, [nk, nj, flen])

        call spherical_to_fourier(spherical, four3, lat_deriv1)

        call fourier_to_grid(four,grd)

        return
    end subroutine spherical_to_grid3D


    subroutine grid_to_spherical2D(grid,spherical)
        real, intent(in) :: grid(:,:)
        type(specVar(n=*,nlev=*)), intent(out) :: spherical

        real :: buff(1,size(grid,1),size(grid,2))

        buff(1,:,:) = grid

        call grid_to_spherical3D(buff, spherical)

        return
    end subroutine grid_to_spherical2D


    subroutine spherical_to_grid2D(spherical,grid,lat_deriv)
        real, intent(out) :: grid(:,:)
        type(specVar(n=*,nlev=*)), intent(in) :: spherical
        logical, intent(in), optional :: lat_deriv

        real :: buff(1,size(grid,1),size(grid,2))

        call spherical_to_grid3D(spherical,buff,lat_deriv)

        grid(:,:) = buff(1,:,:)

        return
    end subroutine spherical_to_grid2D

    subroutine read_specdataGFS(filename,fieldname,dat)

        character (len=*), intent(in) :: filename, fieldname
        type(specVar(nlev=*,n=*)) :: dat

        real :: rebuff(dat%n,dat%nlev), robuff(dat%n-num_fourier/2,dat%nlev)
        real :: iebuff(dat%n,dat%nlev), iobuff(dat%n-num_fourier/2,dat%nlev)
        integer :: i, j, k, nlev, m, nodd
        character(len=len(fieldname)+3) :: iew, rew, iow, row 
        integer :: nlen1, nlen2, ns1, ns2, ne1, ne2

        iew = 'iew'//trim(fieldname)
        rew = 'rew'//trim(fieldname)
        iow = 'iow'//trim(fieldname)
        row = 'row'//trim(fieldname)

        nlev = dat%nlev
        nodd = dat%n-num_fourier/2
        call read_data(filename,rew,rebuff)
        call read_data(filename,iew,iebuff)
        call read_data(filename,row,robuff)
        call read_data(filename,iow,iobuff)

        dat%od = cmplx(0.,0.)
        do k = 1, nlev
            dat%ev(k,:) = cmplx(rebuff(:,k),iebuff(:,k))
            nlen1 = 32; nlen2 = 33
            ne1 = 0; ne2 = 0
            do m = isf, ief
                if (mod(m+1,2)==0) nlen1 = nlen1 - 1
                if (mod(m,2)==0) nlen2 = nlen2 - 1
                ns1 = ne1 + 1; ns2 = ne2 + 1
                ne1 = ns1 + nlen1 - 1; ne2 = ns2 + nlen2 -1
                print *, ns2, nlen2, ns1, nlen1
                dat%od(k,ns2:ns2+nlen1-1) = cmplx(robuff(ns1:ne1,k),iobuff(ns1:ne1,k))
            enddo
            !dat%od(k,1:nodd) = cmplx(robuff(:,k),iobuff(:,k))
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
