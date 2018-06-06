module FFTW3
    
    use, intrinsic :: iso_c_binding

    use fms_mod, only : open_namelist_file, close_file, error_mesg, FATAL, WARNING, NOTE

    include 'fftw3-mpi.f03'

    private

    integer(C_INTPTR_T), parameter :: M=2, M2=M*2
    integer :: rank=2
    integer(C_INTPTR_T) :: howmany3d, nlat, nlev

    type(C_PTR) :: plan3d, c_cout, c_rin
    real(C_DOUBLE), pointer :: rin(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:,:)

    logical :: initialized = .false.
    integer :: plan_level=-1
    character(len=8) :: method='2dc2c'

    public :: init_fftw3, fft, method, fft_1dr2c_serial, fft_1dc2c_serial

    interface fft
        module procedure fft_2dr2c
        module procedure fft_2dc2c
    end interface fft

    namelist/fftw3_nml/plan_level, method

    contains

    subroutine init_fftw3 (n0_in, nlat_in, nlev_in, comm, ni_local, i_start_local, no_local, o_start_local)

        implicit none
        integer, intent(in) :: n0_in, nlat_in, nlev_in, comm
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
        n0 = n0_in

        nlat = nlat_in
        nlev = nlev_in

        howmany3d = nlat*nlev
        
        call fftw_mpi_init()

        if (trim(method)=='2dr2c') then
            rank = 2
            alloc_local = fftw_mpi_local_size_many(rank, n, howmany3d, iblock, comm, &
                                                            local_ni, local_i_start)
            local_no = local_ni
            local_o_start = local_i_start

            c_rin = fftw_alloc_real(2*alloc_local)
            c_cout = fftw_alloc_complex(alloc_local)

            call c_f_pointer(c_rin, rin, [howmany3d, M2, local_ni])
            call c_f_pointer(c_cout, cout, [howmany3d, M, local_no])

            plan3d = fftw_mpi_plan_many_dft_r2c (rank, n, howmany3d, iblock, oblock, &
                                               rin, cout, comm, flags)

        elseif (trim(method)=='2dc2c') then

            rank = 2
            alloc_local = fftw_mpi_local_size_many(rank, n, howmany3d, iblock, comm, &
                                                            local_ni, local_i_start)
            local_no = local_ni
            local_o_start = local_i_start

            c_cout = fftw_alloc_complex(alloc_local)

            call c_f_pointer(c_cout, cout, [howmany3d, M, local_no])
            
            plan3d = fftw_mpi_plan_many_dft(rank, n, howmany3d, iblock, oblock, & 
                                            cout, cout, comm, FFTW_FORWARD, flags)
        endif

        ni_local = local_ni
        i_start_local = local_i_start
        no_local = local_no
        o_start_local = local_o_start
        
        initialized = .true.

        call error_mesg(routine, 'FFTW3 initialized with method :'//trim(method), NOTE)

    end subroutine init_fftw3

    subroutine fft_2dr2c(rinp, coutp)

        real, dimension (:, :, :) :: rinp
        complex, dimension (:, :, :) :: coutp
        
        if (trim(method)/='2dr2c') call error_mesg('fft_2dr2c', 'Method is not 2dr2c', FATAL) 

        rin(:,1:M,:) = rinp(:,1:M,:)

        call fftw_mpi_execute_dft_r2c(plan3d, rin, cout) 

        coutp(:,1:M,:) = cout(:,1:M,:)

    end subroutine fft_2dr2c

    subroutine fft_2dc2c(coutp)
        complex, dimension (:, :, :) :: coutp
        
        if (trim(method)/='2dc2c') call error_mesg('fft_2dc2c', 'Method is not 2dc2c', FATAL) 

        cout(:,1:M,:) = coutp(:,1:M,:)

        call fftw_mpi_execute_dft(plan3d, cout, cout) 

        coutp(:,1:M,:) = cout(:,1:M,:)
    end subroutine fft_2dc2c

    subroutine fft_1dr2c_serial(rinp, coutp)

        real :: rinp(:)
        complex :: coutp(:)
        real(C_DOUBLE), pointer :: rinpl(:)
        complex(C_DOUBLE_COMPLEX), pointer :: coutpl(:)
        type(C_PTR) :: plan1d, data1d
        integer :: L
        L = size(rinp)

        data1d = fftw_alloc_complex(int(L/2+1, C_SIZE_T))
        call c_f_pointer(data1d, rinpl, [2*(L/2+1)])
        call c_f_pointer(data1d, coutpl, [L/2+1])
   
        plan1d = fftw_plan_dft_r2c_1d(L, rinpl, coutpl, FFTW_ESTIMATE) 

        rinpl = rinp
        call fftw_execute_dft_r2c(plan1d, rinpl, coutpl)
        coutp = coutpl
        
        call fftw_destroy_plan(plan1d)
        call fftw_free(data1d)
    end subroutine fft_1dr2c_serial

    subroutine fft_1dc2c_serial(coutp)
        complex :: coutp(:)
        complex(C_DOUBLE_COMPLEX), pointer :: coutpl(:)
        type(C_PTR) :: plan1d, data1d
        integer :: L
        L = size(coutp)

        data1d = fftw_alloc_complex(int(L, C_SIZE_T))
        call c_f_pointer(data1d, coutpl, [L])
   
        plan1d = fftw_plan_dft_1d(L, coutpl, coutpl, FFTW_FORWARD, FFTW_ESTIMATE) 

        coutpl = coutp
        call fftw_execute_dft(plan1d, coutpl, coutpl)
        coutp = coutpl
        
        call fftw_destroy_plan(plan1d)
        call fftw_free(data1d)
    end subroutine fft_1dc2c_serial
 
end module
