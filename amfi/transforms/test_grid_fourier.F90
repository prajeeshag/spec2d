
#ifdef test_grid_to_fourier

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

    use constants_mod, only : PI

    use grid_fourier_mod, only : init_grid_fourier, fft_1dr2c_serial, fft_1dc2c_serial
    use grid_fourier_mod, only : end_grid_fourier, grid_to_fourier, fourier_to_grid, &
                                 fftw_iodim, fft_1dr2c_serial_guru

    implicit none

    type(domain2d) :: domainl
    type(domain2d) :: domainf

    integer, parameter :: reg_nlon = 204
    integer, parameter :: oc_mlon = 224
    integer :: is(2), ie(2), ilen(2), l, i, N, ii
    
    real :: dat1d(oc_mlon) = 0., dat1d1(reg_nlon)
    complex :: fdat1d(oc_mlon)

    call mpp_init() 
    call fms_init()

    is(1) = 1; ie(1) = 204
    is(2) = 205; ie(2) = oc_mlon

    ilen = ie - is + 1

    do l = 1, 2
        N=ie(l)-is(l)+1
        do i = is(l), ie(l)
            ii = i-is(l)
            dat1d(i) = l*cos(2*PI*ii/N)
        enddo
    enddo

    do l = 1, 2
        dat1d1 = 0.
        fdat1d = cmplx(0.,0.)
        dat1d1(1:ilen(l)-1) = dat1d(is(l):ie(l))
        call fft_1dr2c_serial(dat1d1, fdat1d(1:reg_nlon/2+1))
        call write_datac('test_fftw','f1d',fdat1d)    
    enddo

    call fms_io_exit()
    call end_grid_fourier()
    call mpp_exit()

    contains

subroutine write_datac(fnm,fldnm,cdat)
    character(len=*) :: fnm, fldnm
    complex, intent(in) :: cdat(:)
    type(C_PTR) :: cptr
    real, pointer :: rdat(:,:) 

    cptr = C_LOC(cdat)

    call c_f_pointer(cptr, rdat, [2, size(cdat)])

    call write_data(fnm,fldnm,rdat)

    return

end subroutine write_datac
            
end program main

#endif
