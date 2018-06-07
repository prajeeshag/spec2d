module fft_guru
    
    use, intrinsic :: iso_c_binding

    use fms_mod, only : open_namelist_file, close_file, error_mesg, FATAL, WARNING, NOTE

    include 'fftw3-mpi.f03'

    private

    integer(C_INTPTR_T), parameter :: M=2, M2=M*2
    integer :: nplan=0, max_plans=10
    integer, parameter :: rank=2

    type plan_type
        integer(C_INTPTR_T) :: howmany
        type(C_PTR) :: plan, cdat
        type(C_PTR) :: tplan, tcdat
        real(C_DOUBLE), pointer :: rin(:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:,:)
        integer(C_INTPTR_T) :: n0, iblock, oblock
        integer(C_INTPTR_T) :: iblockt, oblockt
        character(len=8) :: method
    endtype plan_type

    integer(C_INTPTR_T) :: NTOTAL, FTOTAL, FTRUNC, NLOCAL, FLOCAL
    real :: RSCALE
    integer :: COMM
    integer :: plan_level = 1, plan_flags
    logical :: transpos=.true.
    logical :: initialized=.false.

    character (len=16), parameter :: modul = 'fft_guru'

    type(plan_type), allocatable :: myplans(:)

    public :: init_fft_guru, fft, fft_1dr2c_serial, fft_1dc2c_serial

    namelist/fft_guru_nml/plan_level, transpos, max_plans

    contains

    subroutine init_fft_guru (nlon, nlon_local, nfourier, nfourier_local, comm_in)
        implicit none

        integer, intent(in) :: nlon, comm_in, nlon_local, nfourier, nfourier_local
        integer(C_INTPTR_T) :: i, j, alloc_local, local_ni, local_i_start, local_no, local_o_start
        character (len=32) :: routine = 'init_fft_guru'
        integer :: unit, flags, stat

        unit = open_namelist_file()
        read(unit, nml=fft_guru_nml,iostat=stat)
        call close_file(unit)

        select case (plan_level)
        case(0)
            plan_flags = FFTW_ESTIMATE
        case(1)
            plan_flags = FFTW_MEASURE
        case(2)
            plan_flags = FFTW_PATIENT
        case(3)
            plan_flags = FFTW_EXHAUSTIVE
        case default
            call error_mesg(routine,'Wrong option for plan_level, &
            & set plan_level (accepted values are 0-3) in fft_guru_nml', FATAL)
        end select

        if(mod(nlon,2)) call error_mesg(routine,'NLON must be an integer',FATAL)

        NTOTAL = nlon
        NLOCAL = nlon_local
        FTOTAL = nlon/2
        FTRUNC = nfourier 
        FLOCAL = nfourier_local 
        COMM = comm_in
        RSCALE = 1./(nlon*2)

        call fftw_mpi_init()

        allocate(myplans(max_plans))
        
        initialized = .true.

        call error_mesg(routine, 'FFT_GURU initialized with method !', NOTE)

    end subroutine init_fft_guru

    subroutine fft(rinp, coutp, id, execute)
        real, intent(in) :: rinp(:,:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:,:)
        integer, intent(inout), optional :: id
        logical, optional :: execute
        integer(C_INTPTR_T) :: howmany, nx, ny
        integer :: idl

        nx = size(rinp,1); ny=size(rinp,2)
        howmany=nx*ny/2

        idl = 0
        if(present(id).and.id>0) idl=id

        if (idl<1) then
            do i = 1, nplan
                if (myplans(i)%howmany==howmany) then
                    idl = i
                    exit
                endif
            enddo
        endif
        
        if (idl<1) then
            idl = register_plan(howmany)
        endif

        if(present(id)) id = idl

        if(present(execute).and..not.(execute)) return

        myplans(idl)%rin(:,1:2,:) = reshape(rinp(:,:,:)*RSCALE,shape=[howmany,M,NLOCAL])

        call fftw_mpi_execute_dft_r2c(myplans(idl)%plan, myplans(idl)%rin, myplans(idl)%cout) 

        coutp(:,:,:) = reshape(myplans(idl)%cout(:,:,:),shape=[nx,ny,NLOCAL])

    end subroutine fft

    function register_plan(howmany)
        implicit none
        integer(C_INTPTR_T), intent(in) :: howmany
        integer :: register_plan, n, flags
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
        integer(C_INTPTR_T) :: alloc_local

        if (howmany<1) call error_mesg('register_plan', 'howmany cannot be Zero', FATAL)

        nplan = nplan + 1
        
        if(nplan>max_plans) call error_mesg(modul,'No: plans > Max_plans', FATAL)

        n = nplan

        myplans(n)%howmany = howmany
        
        !alloc_local = fftw_mpi_local_size_many_transposed(rank, (/NTOTAL,M/), howmany, &
        !                   NLOCAL, FFTW_MPI_DEFAULT_BLOCK, comm, local_n0, local_0_start, &
        !                   local_n1, local_1_start)

        alloc_local = fftw_mpi_local_size_many(rank, (/NTOTAL,M/), howmany, &
                           NLOCAL, COMM, local_n0, local_0_start)
  
        myplans(n)%cdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(myplans(n)%cdat, myplans(n)%rin, [howmany, M2, local_n0])
        call c_f_pointer(myplans(n)%cdat, myplans(n)%cout, [howmany, M, local_n0])

        flags = plan_flags
        myplans(n)%plan = fftw_mpi_plan_many_dft_r2c (rank, (/NTOTAL,M/), howmany, NLOCAL, &
                            NLOCAL, myplans(n)%rin, myplans(n)%cout, COMM, flags)
         
        register_plan = n
    end function register_plan

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
