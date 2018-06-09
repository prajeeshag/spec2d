module grid_to_fourier_mod

    use fft_guru, only : init_fft_guru, fft, fft_1dr2c_serial, fft_1dc2c_serial, end_fft_guru

    use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
    use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
    use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
    use mpp_mod, only : mpp_sync, mpp_root_pe
    use mpp_domains_mod, only : mpp_define_domains, domain1d, mpp_get_compute_domain
    use fms_mod, only : read_data, write_data, open_namelist_file, close_file
    use fms_io_mod, only : fms_io_exit 

    implicit none
    integer, parameter :: sp = kind(1e0), dp = kind(1d0)    
    
    type(domain1d) :: domainl
    type(domain1d) :: domainf
 
    contains


    subroutine init_grid_to_fourier(nlon, nlat, nlev)

        integer, intent(in) :: nlon, nlat, nlev
      
        integer :: ilen, istart, olen, ostart, nlonb2
        
        integer, allocatable :: pelist(:) 
 
        character(len=32) :: routine='init_grid_to_fourier'

        integer :: comm, idfft3d, n, nt=1
        logical :: check=.false.
        real, allocatable :: fld(:,:,:), fld1d(:), fld1dout(:)
        complex, allocatable :: fldc1d(:,:), fldc(:,:,:,:)
        integer :: isc, iec, isg, ieg, m, l, t, i, ig, k, kstart=0, kend=0, kstep=1
        integer :: isf, ief, flen
        real :: scl, x, y, imgf=0.3, phi=0.15
        integer :: clck_fftw3, clck_drcft, init, unit, cl, ck, num_fourier
        complex :: wgt, a, b, c, d, e
        complex(kind=4) :: cpout(3)
        logical :: oddeven=.true., ideal_data=.false.
        real, parameter :: PI=4.D0*DATAN(1.D0)
        complex, parameter :: ui = cmplx(0.,1.), mui = -1.*ui
        integer :: pe

        namelist/grid_to_fourier_nml/kstart, kend, kstep, ideal_data, imgf, ck, cl, check, num_fourier, nt

        call mpp_init() 
        pe = mpp_pe()
        num_fourier = nlon
 
        unit = open_namelist_file()
        read(unit,nml=grid_to_fourier_nml)
        call close_file(unit)

        clck_fftw3 = mpp_clock_id('fftw3')

        allocate(pelist(mpp_npes()))

        call mpp_define_domains( (/1,nlon/), mpp_npes(), domainl, halo=0)
        call mpp_get_compute_domain(domainl, isc, iec)
        ilen = iec-isc+1

        call mpp_define_domains( (/1,num_fourier/), mpp_npes(), domainf, halo=0)
        call mpp_get_compute_domain(domainf, isf, ief)
        flen = ief-isf+1

        call mpp_get_current_pelist(pelist,commid=comm)

        !print *, 'comm=', comm
        call init_fft_guru(nlon, ilen, num_fourier, isf, flen, comm, nlev, nlat)
        ief = isf+flen-1
        print *, 'pe, isf, ief, flen=', mpp_pe(), isf, ief, flen

!        call mpp_sync()
!        call mpp_error('grid_to_fourier','TESTING...',FATAL)

        allocate(fld(nlev,nlat,isc:iec))
        allocate(fldc(isf:ief,2,nlat,nlev))

        allocate(fld1d(1:nlon))
        allocate(fld1dout(1:nlon))
        allocate(fldc1d(1:nlon,2))
        
        if(.not.ideal_data) call read_data('test_data.nc', 'tas', fld1d)

        scl=1./nlon

        if (ideal_data) fld1d=0.
        do i = 0, nlon-1
            if (ideal_data) then
                do k = kstart, kend, kstep
                    fld1d(i+1) = fld1d(i+1) + 2*k*cos(2.*PI*(i-phi)*real(k)/nlon)
                enddo
            endif
         enddo

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
        
        fldc=cmplx(0.,0.)
      
        call mpp_sync()

        call mpp_clock_begin(clck_fftw3)
        do t = 1, nt
            call fft(fld,fldc)
        enddo
        call mpp_clock_end(clck_fftw3)

        if(cl>nlat/2) cl = nlat/2
        if(cl<1) cl = 1
        if(ck>nlev) ck = nlev
        if(ck<1) ck = 1

        if (check.or.mpp_npes()==1) then
        if(mpp_pe()==mpp_root_pe()) print *, 'printing for lev and lat index :',ck, cl
        do i = isf, ief
            cpout(1) = fldc(i,1,cl,ck) 
            cpout(2) = fldc(i,2,cl,ck) 
            !cpout(1:2) = fldc(1,:,i+1)
            if (ideal_data) print *, i, cpout(1:2)
        enddo
        endif

        if (mpp_pe()==mpp_root_pe().and.check) then
            k = (cl+(ck-1)*nlat/2)
            fld1dout(:) = k + k*fld1d(:)
            call fft_1dr2c_serial(fld1dout(:)*scl,fldc1d(:,1))
            fld1dout(:) = k - k*fld1d(:)
            call fft_1dr2c_serial(fld1dout(:)*scl,fldc1d(:,2))
            if (ideal_data) then
            print *, ''
            print *, ''
            print *, 'with 1dr2c'
            print *, ''
            do i = 0, nlon/2
                cpout(1) = fldc1d(i+1,1) + fldc1d(i+1,2)
                cpout(2) = fldc1d(i+1,1) - fldc1d(i+1,2)
                print *, i, cpout(1:2)
            enddo
            print *, '1dr2c'
            print *, ''
            print *, ''
            endif 
        endif

        call fms_io_exit()
        call end_fft_guru()
        call mpp_exit()
            
    end subroutine init_grid_to_fourier

end module grid_to_fourier_mod

