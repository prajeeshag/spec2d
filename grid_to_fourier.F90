module grid_to_fourier_mod

    use FFTW3

    use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
    use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
    use mpp_mod, only : mpp_exit
 
    implicit none
    
     
    contains


    subroutine init_grid_to_fourier(nlon, nlat, nlev)

        integer, intent(in) :: nlon, nlat, nlev
      
        integer :: ilen, istart, olen, ostart, nlonb2
        
        integer, allocatable :: pelist(:) 
 
        character(len=32) :: routine='init_grid_to_fourier'
    
        integer :: comm
 
        call mpp_init() 
       
        if (mod(nlon,2)/=0) call mpp_error(routine, 'nlon is not an ever number', FATAL)
        
        nlonb2 = nlon/2

        allocate(pelist(mpp_npes()))

        call mpp_get_current_pelist(pelist,commid=comm)
        
        call init_fftw3 (nlonb2, nlat, nlev, comm, ilen, istart, olen, ostart)
            
        print *, mpp_pe(), ilen, istart, olen, ostart

        call mpp_exit()
 
    end subroutine init_grid_to_fourier

end module grid_to_fourier_mod

