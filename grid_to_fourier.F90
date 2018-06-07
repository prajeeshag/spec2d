module grid_to_fourier_mod

    use fft_guru, only : init_fft_guru, fft, fft_1dr2c_serial, fft_1dc2c_serial, end_fft_guru

    use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
    use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
    use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
    use mpp_mod, only : mpp_sync
    use mpp_domains_mod, only : mpp_define_domains, domain1d, mpp_get_compute_domain
    use fms_mod, only : read_data, write_data, open_namelist_file, close_file
    use fms_io_mod, only : fms_io_exit 

    implicit none
    integer, parameter :: sp = kind(1e0), dp = kind(1d0)    
    
    type(domain1d) :: domain
 
    contains


    subroutine init_grid_to_fourier(nlon, nlat, nlev)

        integer, intent(in) :: nlon, nlat, nlev
      
        integer :: ilen, istart, olen, ostart, nlonb2
        
        integer, allocatable :: pelist(:) 
 
        character(len=32) :: routine='init_grid_to_fourier'

        integer :: comm, idfft3d, n

        real, allocatable :: fld(:,:,:), fld1d(:), fld1dout(:)
        complex, allocatable :: fldc1d(:,:), fldc(:,:,:), fldc1(:,:,:)
        real :: AUX1CRS(42002)
        complex, allocatable :: four_fld(:,:,:)
        integer :: isc, iec, isg, ieg, m, l, t, i, ig, k, kstart=0, kend=0, kstep=1
        real :: scl, x, y, imgf=0.3, phi=0.15
        integer :: clck_fftw3, clck_drcft, init, unit
        complex :: wgt, a, b, c, d, e
        complex(kind=4) :: cpout(3)
        logical :: oddeven=.true., ideal_data=.false.
        real, parameter :: PI=4.D0*DATAN(1.D0)
        complex, parameter :: ui = cmplx(0.,1.), mui = -1.*ui
        complex, parameter :: cone = cmplx(1.,0.), w1b2=exp(mui*2.*pi*0.5)
        complex, parameter :: one_w1b2 = cone - w1b2
        complex :: wlbN

        namelist/grid_to_fourier_nml/kstart, kend, kstep, oddeven, ideal_data, imgf

        call mpp_init() 
 
        unit = open_namelist_file()
        read(unit,nml=grid_to_fourier_nml)
        call close_file(unit)

        clck_fftw3 = mpp_clock_id('fftw3')

        allocate(pelist(mpp_npes()))


        call mpp_define_domains( (/1,nlon/), mpp_npes(), domain, halo=0)
        
        call mpp_get_compute_domain(domain, isc, iec)
       
        ilen = iec-isc+1

        call mpp_get_current_pelist(pelist,commid=comm)

        !print *, 'comm=', comm
        call init_fft_guru(nlon, ilen, nlon, ilen, comm, nlev, nlat)

        allocate(fld(nlev,nlat,isc:iec))
        allocate(fldc(nlev,nlat,isc:iec))
        allocate(fld1d(1:nlon))
        allocate(fldc1d(2,1:nlon))
        
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

         do l = 1, nlev
            do m = 1, nlat
               do i = isc, iec
                   fld(l,m,i) = m*fld1d(i)
               enddo
            end do
        end do
      
        
        fldc=cmplx(0.,0.)
      
        call mpp_sync()

        call mpp_clock_begin(clck_fftw3)
        do t = 1, 100
            call fft(fld, fldc)
        enddo
        call mpp_clock_end(clck_fftw3)

        do k = isc, iec
            i = k - 1
            cpout(1) = (fldc(1,1,i+1) + fldc(1,2,i+1))*0.5
            cpout(2) = (fldc(1,1,i+1) - fldc(1,2,i+1))*0.5
            !cpout(1:2) = fldc(1,:,i+1)
            if (ideal_data) print *, i, cpout(1:2)
        enddo

        if (mpp_npes()==1) then

            call fft_1dr2c_serial(fld(1,1,:)*scl,fldc1d(1,:))
            call fft_1dr2c_serial(fld(1,2,:)*scl,fldc1d(2,:))

            if (ideal_data) then
            print *, '1dr2c'
            do i = 0, nlon
                cpout(1) = fldc1d(1,i+1)
                cpout(2) = fldc1d(2,i+1)
                print *, i, cpout(1:2)
            enddo
            print *, '1dr2c'
            print *, ''
            print *, ''
            endif 
        endif

        call end_fft_guru()
        call fms_io_exit()
        call mpp_exit()
            
    end subroutine init_grid_to_fourier

end module grid_to_fourier_mod

