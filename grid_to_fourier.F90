module grid_to_fourier_mod

    use FFTW3, only : init_fftw3, fft, method, fft_1dr2c_serial, fft_1dc2c_serial

    use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
    use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
    use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
    use mpp_mod, only : mpp_sync
    use fms_mod, only : read_data, write_data, open_namelist_file, close_file
    use fms_io_mod, only : fms_io_exit 

    implicit none
    integer, parameter :: sp = kind(1e0), dp = kind(1d0)    
     
    contains


    subroutine init_grid_to_fourier(nlon, nlat, nlev)

        integer, intent(in) :: nlon, nlat, nlev
      
        integer :: ilen, istart, olen, ostart, nlonb2
        
        integer, allocatable :: pelist(:) 
 
        character(len=32) :: routine='init_grid_to_fourier'

        integer :: comm

        real, allocatable :: fld(:,:,:), fld1d(:), fld1dout(:)
        complex, allocatable :: fldc1d(:), fldc(:,:,:), fldc1(:,:,:)
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
 
        if (mod(nlon,2)/=0) call mpp_error(routine, 'nlon is not an ever number', FATAL)
        
        nlonb2 = nlon/2

        if (mod(nlonb2,mpp_npes())/=0)  &
         call mpp_error(routine, 'No: procs in x-dir is not a factor of nlon/2', FATAL)

        unit = open_namelist_file()
        read(unit,nml=grid_to_fourier_nml)
        call close_file(unit)

        clck_fftw3 = mpp_clock_id('fftw3')

        allocate(pelist(mpp_npes()))

        call mpp_get_current_pelist(pelist,commid=comm)
        
        call init_fftw3 (nlonb2, nlat, nlev, comm, ilen, istart, olen, ostart)
        
        isc = istart + 1; iec = isc + ilen - 1

        allocate(fld(nlat*nlev,2,isc:iec))
        allocate(fldc(nlat*nlev,2,isc:iec))
        allocate(fldc1(nlat*nlev,2,isc:iec))
        fld = 0.
        allocate(four_fld(nlat*nlev,2,isc:iec))
        allocate(fld1d(nlon))
        allocate(fld1dout(nlon))
        allocate(fldc1d(nlon))
        
        if(.not.ideal_data) call read_data('test_data.nc', 'tas', fld1d)

        scl=1./nlon

         if (ideal_data) fld1d=0.
         do i = 0, nlon-1
            if (ideal_data) then
                do k = kstart, kend, kstep
                    fld1d(i+1) = fld1d(i+1) + 2*k*cos(2.*PI*(i-phi)*real(k)/nlon)
                enddo
            endif
            fldc1d(i+1) = cmplx(fld1d(i+1),fld1d(i+1)*imgf)*scl
         enddo

        do l = 1, nlev*nlat
            do m = 1, 2
               do i = isc, iec
                   ig = 2*(i-1)+m
                   fld(l,m,i) = fld1d(ig)*scl
                   fldc1(l,m,i) = cmplx(fld1d(ig),fld1d(ig)*imgf)*scl
               enddo
            end do
        end do
      
        call mpp_sync()
 
        
       ! print *, 'testing 2D method'
       ! print *, '' 
       ! print *, 'X1...'
       ! fldc = fldc1
       ! do m = 1, size(fldc,2)
       !     call fft_1dc2c_serial(fldc(1,m,:))
       ! enddo

       ! do i = 0, nlon/2-1
       !     cpout(1) = (fldc(1,1,i+1)) !* exp(-1.*ui*2.*pi*real(i)/nlon) 
       !     cpout(2) = (fldc(1,2,i+1)) !* exp(-1.*ui*2.*pi*real(i)/nlon)
       !     print *, i, cpout(1:2)
       ! enddo

       ! print *, '' 
       ! print *, 'X2...'
       ! print *, 'ui=', ui
       ! do l = 1, size(fldc,3)
       !     print *, 'wgt=', wgt
       !     do i = 1, size(fldc,2)
       !         wgt = exp((-1.*ui*2.*PI/nlon)*(l-1)*(i-1))
       !         fldc(1,i,l) = fldc(1,i,l) * wgt
       !     enddo
       !     call fft_1dc2c_serial(fldc(1,:,l))
       ! enddo

       ! do i = 0, nlon/2
       !     cpout(1) = (fldc(1,1,i+1)) !* exp(-1.*ui*2.*pi*real(i)/nlon) 
       !     cpout(2) = (fldc(1,2,i+1)) !* exp(-1.*ui*2.*pi*real(i)/nlon)
       !     print *, i, cpout(1:2)
       ! enddo

       ! print *, 'end testing 2d method'
       ! print *, '' 

       
        call mpp_clock_begin(clck_fftw3)
        if (trim(method)=='2dc2c') then
            fldc = fldc1
            call fft(fldc)
            do i = 0, nlonb2-1
                wlbN = exp(mui*2.*pi*real(i)/nlon)
                a = wlbN - w1b2
                b = cone - wlbN
                c = cone - w1b2
                d = cone - w1b2*wlbN
                e = w1b2*b
                cpout(1) = (a*fldc(1,1,i+1) + b*fldc(1,2,i+1))/c
                cpout(2) = (d*fldc(1,2,i+1) - e*fldc(1,1,i+1))/c
                if (ideal_data) print *, i, cpout(1:2)
            !    cpout(1:2) = fldc(1,:,i+1)
            !    if (ideal_data) print *, i, cpout(1:2)
                fldc1(1,:,i+1) = cpout(1:2)
            enddo
            call write_data('fftw_out.nc','realc', real(fldc1(1,:,:)))
            call write_data('fftw_out.nc','imagc', dimag(fldc1(1,:,:)))
            call write_data('fftw_out.nc','magc', abs(fldc1(1,:,:)))
           ! print *, '' 
           ! print *, '2dc2c'
           ! print *, '2dc2c'
           ! print *, ''
           ! print *, ''
        endif
        
        fldc=cmplx(0.,0.)
        if (trim(method)=='2dr2c') then
            call fft(fld, fldc)
            do i = 0, nlonb2-1
                wlbN = exp(mui*2.*pi*real(i)/nlon)
                a = wlbN - w1b2
                b = cone - wlbN
                c = cone - w1b2
                d = cone - w1b2*wlbN
                e = w1b2*b
                cpout(1) = (a*fldc(1,1,i+1) + b*fldc(1,2,i+1))/c
                cpout(2) = (d*fldc(1,2,i+1) - e*fldc(1,1,i+1))/c
                if (ideal_data) print *, i, cpout(1:2)
                cpout(1:2) = fldc(1,:,i+1)
            !    if (ideal_data) print *, i, cpout(1:2)
            !    fldc1(1,:,i+1) = cpout(1:2)
            enddo
            call write_data('fftw_out.nc','realc', real(fldc1(1,:,:)))
            call write_data('fftw_out.nc','imagc', dimag(fldc1(1,:,:)))
            call write_data('fftw_out.nc','magc', abs(fldc1(1,:,:)))
           ! print *, '' 
           ! print *, '2dc2c'
           ! print *, '2dc2c'
           ! print *, ''
           ! print *, ''
        endif
        call mpp_clock_end(clck_fftw3)

        if (mpp_npes()==1) then
           ! clck_drcft = mpp_clock_id('drcft')
           ! INIT=1
           ! CALL DCRFT(INIT, fld1d, 1, fld1dout, 1, nlon, 1, 1, 1./nlon, & 
           !             AUX1CRS, 22000, AUX1CRS(22001), 20000)

           ! call mpp_clock_begin(clck_drcft)
           ! INIT=0
           ! CALL DCRFT(INIT, fld1d, 1, fld1dout, 1, nlon, 1, 1, 1./nlon, & 
           !             AUX1CRS, 22000, AUX1CRS(22001), 20000)
           ! call mpp_clock_end(clck_drcft)
           ! print *, 'write data dcrft_out'
           ! call write_data('dcrft_out.nc','fld', fld1dout(:))
                
            !call fft_1dr2c_serial(fld1d*scl, fldc1d) 
            call fft_1dc2c_serial(fldc1d(:)) 
            call write_data('fftw_out.nc','real', real(fldc1d(:)))
            call write_data('fftw_out.nc','imag', dimag(fldc1d(:)))
            call write_data('fftw_out.nc','mag', abs(fldc1d(:)))

            if (ideal_data) then
            print *, '1dr2c'
            do i = 0, nlon
                cpout(1) = fldc1d(i+1)
                print *, i, cpout(1)
            enddo
            print *, '1dr2c'
            print *, ''
            print *, ''
            !print *, ''
            endif 
        endif

        call fms_io_exit()
        call mpp_exit()
            
    end subroutine init_grid_to_fourier

end module grid_to_fourier_mod

