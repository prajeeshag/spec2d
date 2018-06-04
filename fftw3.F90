module FFTW3

    use, intrinsic :: iso_c_binding

    use fms_mod, only : open_namelist_file, close_file, error_mesg, FATAL, WARNING, NOTE

    include 'fftw3-mpi.f03'

    type(C_PTR) :: plan3d, c_cout, c_rin
    real(C_DOUBLE), pointer :: rin(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:)
    integer(C_INTPTR_T), parameter :: rank=2, M=2, M2=M*2

    integer :: howmany3d

    logical :: initialized = .false.
    integer :: plan_level=-1

    namelist/fftw3_nml/plan_level

    contains

    subroutine init_fftw3 (n0_in, nlat, nlev, comm, ni_local, i_start_local, no_local, o_start_local)

        implicit none
        integer, intent(in) :: n0_in, nlat, nlev, comm
        integer, intent(out) :: ni_local, i_start_local, no_local, o_start_local
        integer(C_INTPTR_T) :: n0, iblock, oblock, n(rank)
        integer(C_INTPTR_T) :: i, j, alloc_local, local_ni, local_i_start, local_no, local_o_start
        integer :: unit, flags
        character (len=32) :: routine = 'init_fftw3'

        unit = open_namelist_file()
        read(unit, nml=fftw3_nml)
        call close_file(unit)

        select case (plan_level)
        case(0)
            flags = FFTW_ESTIMATE
        case(1)
            flags = FFTW_MEASURE
        case(2)
            flags = FFTW_PATIENT
        case(3)
            flags = FFTW_EXHAUSTIVE
        case default
            call error_mesg(routine,'Wrong option for plan_level, &
            & set plan_level (accepted values are 0-3) in fftw3_nml', FATAL)
        end select

        iblock = FFTW_MPI_DEFAULT_BLOCK
        oblock = FFTW_MPI_DEFAULT_BLOCK

        n(1) = n0_in
        n(2) = M

        howmany3d = nlat*nlev
        
        call fftw_mpi_init()

        alloc_local = fftw_mpi_local_size_many(rank, n, howmany3d, iblock, comm, &
                                                        local_ni, local_i_start)
        local_no = local_ni
        local_o_start = local_i_start

        c_rin = fftw_alloc_real(2*alloc_local)
        c_cout = fftw_alloc_complex(alloc_local)

        print *, 'alloc_local=', alloc_local

        call c_f_pointer(c_rin, rin, [M2, local_ni])
        call c_f_pointer(c_cout, cout, [M, local_no])

        plan3d = fftw_mpi_plan_many_dft_r2c (rank, n, howmany, iblock, oblock, &
                                           rin, cout, comm, flags)

        ni_local = local_ni
        i_start_local = local_i_start
        no_local = local_no
        o_start_local = local_o_start
        
        initialized = .true.

        call error_mesg(routine, 'FFTW3 initialized', NOTE)

    end subroutine init_fftw3


   ! subroutine fft_3d(rinp,coutp)

   !     real, dimension(2,:) ::  
   !      
end module
