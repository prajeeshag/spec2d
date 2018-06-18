
#ifdef test_grid_to_fourier

program main

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

    use fourier_spherical_mod, only : init_fourier_spherical, fourier_to_spherical

    implicit none

    type(domain2d) :: domainl
    type(domain2d) :: domainf

    logical :: fpe 
    integer :: nlon, nlat, nlev
      
    integer :: ilen, istart, olen, ostart, nlonb2
        
    integer, allocatable :: pelist(:), extent(:)
    integer, allocatable :: fpelist(:), fextent(:)
    integer, allocatable :: Tshuff(:)
 
    character(len=32) :: routine='test_grid_to_fourier'

    integer :: comm, idfft3d, n, nt=1
    logical :: check=.false.
    real, allocatable :: fld(:,:,:), fld1d(:), fld1dout(:), fldout(:,:,:)
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

    namelist/test_grid_to_fourier_nml/kstart, kend, kstep, ideal_data, imgf, ck, cl, &
                                      check, num_fourier, nt, nlon, nlat, nlev, debug

    call mpp_init() 
    call fms_init()

    pe = mpp_pe()

    unit = open_namelist_file()
    read(unit,nml=test_grid_to_fourier_nml)
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
    call init_grid_fourier(nlon, ilen, num_fourier, isf, flen, comm, nlev, nlat, Tshuff=Tshuff)

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
        if(flen /= ief-isf+1) call mpp_error('test_grid_to_fourier', 'flen /= ief-isf+1', FATAL)
        print *, 'pe, isf, flen, load=', mpp_pe(), isf, flen, sum(num_fourier-Tshuff(isf:ief)+1)

        call init_fourier_spherical(num_fourier, num_fourier+1, nlat, domainf, Tshuff) 
    endif
    call mpp_set_current_pelist()

    allocate(fld(nlev,jsc:jec,isc:iec))
    allocate(fldout(nlev,jsc:jec,isc:iec))
    allocate(fldct(jsf:jef,nlev,isf:ief))

    allocate(fld1d(1:nlon))
    allocate(fld1dout(1:nlon))
    allocate(fldc1d(0:nlon-1,2))
    
    if(.not.ideal_data) call read_data('test_data.nc', 'tas', fld1d)

    scl=1./nlon

    if (ideal_data) then
        kend = num_fourier-1
        fld1d=0.
         do i = 0, nlon-1
            if (ideal_data) then
                do k = kstart, kend, kstep
                    fld1d(i+1) = fld1d(i+1) + 2*k*cos(2.*PI*(i-phi)*real(k)/nlon)
                enddo
            endif
         enddo
     endif

     k = 0
     do l = 1, nlev
        do m = 1, nlat/2
            k = k + 1
           do i = isc, iec
               fld(l,m,i) = k + k*fld1d(i)
               fld(l,nlat/2+m,i) = k - k*fld1d(i)
               !fld(l,m,i) = (m+(l-1)*nlat)
           enddo
        enddo
    enddo
   
    call write_data('test_grid2four', 'fld', fld(:,:,:), domain=domainl)
    call write_data('test_grid2four', 'fld1', fld(1,:,:), domain=domainl)

    !idrestart = register_restart_field_r2d(fileObj, filename, fieldname, data, domain)
 
    call mpp_sync()


    !call mpp_clock_begin(clck_grid2fourier)
    do t = 1, nt
        call grid_to_fourier(fld, fldct)
    enddo
    !call mpp_clock_end(clck_grid2fourier)

    if (fpe) then
        call mpp_set_current_pelist(fpelist)
        call write_data('test_grid2four', 'fldc_r', real(fldct(:,1,:)), domain=domainf)
    endif
    call mpp_set_current_pelist()

    !call mpp_clock_begin(clck_fourier2grid)
    do t = 1, nt
        call fourier_to_grid(fldct, fldout)
    enddo
    !call mpp_clock_end(clck_fourier2grid)


    if (ideal_data) then
        k = 0
        if (mpp_pe()==mpp_root_pe()) then
            print *, ''
            print *, ''
            print *, 'with 1dr2c'
            print *, ''
            print *, 'printing for lev and lat index :',ck, cl
        endif
        do l = 1, nlev
        do m = 1, nlat/2
            k = k+1
                fld1dout(:) = k + k*fld1d(:)
            call fft_1dr2c_serial(fld1dout(:)*scl,fldc1d(:,1))
            if (debug) then 
                call mpp_sync()  
                if (mpp_pe()==mpp_root_pe()) then
                   print *, 'full= ', real(fldc1d(1:num_fourier,1))
                   print *, 'half= ', Tshuff(isf:ief), real(fldct(m,l,:))
                endif
            endif
            call mpp_sync()  
            do i = isf, ief
                j = Tshuff(i)
                cpout(1) = fldc1d(j,1)
                cpout(2) = fldct(m,l,i)
                cpout(3) = fldct(nlat/2+m,l,i)
                if(abs(cpout(1)-cpout(2))>1.e-10) then
                    print *,'forward check:', k, j, i, cpout(1),cpout(1)
                    call mpp_error('test_grid_to_fourier','forward check error', WARNING)
                endif
            enddo
        enddo
        enddo
    endif

    if (ideal_data) then
        k = 0
        if (mpp_pe()==mpp_root_pe()) then
            print *, ''
            print *, ''
            print *, 'Backward'
            print *, ''
        endif
        do m = 1, nlat
            do l = 1, nlev
                if (debug) then
                    call mpp_sync()
                    print '(A,1x,2(I3,1x),100(F13.6,1x))', 'backward check1:', l, m, fld(l,m,isc:iec)
                    print '(A,1x,2(I3,1x),100(F13.6,1x))', 'backward check2:', l, m, fldout(l,m,isc:iec)
                endif
                call mpp_sync()
                do i = isc, iec
                    if(abs(fldout(l,m,i)-fld(l,m,i))>1.e-6) then
                        print *,'backward check:', l, m, i, fldout(l,m,i), fld(l,m,i)
                        call mpp_error('test_grid_to_fourier','backward check error', FATAL)
                    endif
                enddo
            enddo
        enddo
    endif

    call fms_io_exit()
    call end_grid_fourier()
    call mpp_exit()
            
end program main

#endif
