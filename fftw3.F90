module fft_guru
    
    use, intrinsic :: iso_c_binding

    use fms_mod, only : open_namelist_file, close_file, error_mesg, FATAL, WARNING, NOTE

    implicit none

    include 'fftw3-mpi.f03'

    private

    integer :: nplan=0
    integer, parameter :: max_plans=3
    integer, parameter :: rank=2

    type plan_type
        integer(C_INTPTR_T) :: howmany
        type(C_PTR) :: plan, cdat, rdat
        type(C_PTR) :: tplan, tcdat
        real(C_DOUBLE), pointer :: rin(:,:,:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:,:,:,:)
        integer(C_INTPTR_T) :: n0, iblock, oblock
        integer(C_INTPTR_T) :: iblockt, oblockt
        character(len=8) :: method
    endtype plan_type

    integer(C_INTPTR_T) :: NLON, FTOTAL, FTRUNC, NLON_LOCAL, FLOCAL
    integer(C_INTPTR_T) :: NLEV, NLAT, NVAR
    integer :: COMM_FFT

    real :: RSCALE
    integer :: plan_level = 1, plan_flags, id2d, id3d, id3dext
    logical :: transpos=.true.
    logical :: initialized=.false.

    character (len=16), parameter :: modul = 'fft_guru'

    type(plan_type) :: myplans(max_plans)

    interface fft
        module procedure fft3d
        module procedure fft2d
    end interface

    public :: init_fft_guru, fft, fft_1dr2c_serial, fft_1dc2c_serial, end_fft_guru

    namelist/fft_guru_nml/plan_level, transpos

    contains

    subroutine init_fft_guru (nlons, nlons_local, nfourier, nfourier_local, comm_in, nlevs, nlats, nvars)
        implicit none

        integer, intent(in) :: nlons, comm_in, nlons_local, nfourier, nfourier_local, nlats, nlevs
        integer, intent(in), optional :: nvars
        integer(C_INTPTR_T) :: alloc_local, local_ni, local_i_start, local_no, local_o_start, nvar1
        integer(C_INTPTR_T) :: howmany
        character (len=32) :: routine = 'init_fft_guru'
        integer :: unit, flags, stat

        unit = open_namelist_file()
        read(unit, nml=fft_guru_nml,iostat=stat)
        call close_file(unit)

        NLON = nlons
        NLAT = nlats
        NLEV = nlevs
        NVAR = 1 
        if (present(nvars))NVAR=nvars
        COMM_FFT = comm_in
        NLON_LOCAL = nlons_local
        FTOTAL = nlons/2
        FTRUNC = nfourier 
        FLOCAL = nfourier_local 
        RSCALE = 1./nlons

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

        if(mod(NLON,2)) call error_mesg(routine,'NLON must be an integer',FATAL)

        call fftw_mpi_init()

        !plan 2d
        id2d = register_plan(1,1,NLAT,NLON_LOCAL,comm_in)

        !plan 3d
        id3d = register_plan(1,NLEV,NLAT,NLON_LOCAL,comm_in)
        
        !plan 3d
        id3dext = register_plan(NVAR,NLEV,NLAT,NLON_LOCAL,comm_in)
        
        initialized = .true.

        call error_mesg(routine, 'FFT_GURU initialized !!!', NOTE)

    end subroutine init_fft_guru

    function register_plan(nvars, nlevs, nlats, nlons_local, comm_in)
        implicit none
        integer(C_INTPTR_T), intent(in) :: nvars, nlevs, nlats, nlons_local
        integer, intent(in) :: comm_in
        integer(C_INTPTR_T) :: howmany
        integer :: register_plan, n, flags
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
        integer(C_INTPTR_T) :: alloc_local, n0(2)

        howmany = nvars*nlevs*nlats/2

        if (howmany<1) call error_mesg('register_plan', 'howmany cannot be Zero', FATAL)

        nplan = nplan + 1
        
        if(nplan>max_plans) call error_mesg(modul,'No: plans > Max_plans', FATAL)

        n = nplan

        myplans(n)%howmany = howmany
       
        n0=(/NLON,2/) 
        !alloc_local = fftw_mpi_local_size_many_transposed(rank, (/NLON,M/), howmany, &
        !                   NLON_LOCAL, FFTW_MPI_DEFAULT_BLOCK, comm, local_n0, local_0_start, &
        !                   local_n1, local_1_start)

        
        alloc_local = fftw_mpi_local_size_many(rank, n0, howmany, &
                           nlons_local, comm_in, local_n0, local_0_start)
  
        myplans(n)%cdat = fftw_alloc_complex(alloc_local)
        !myplans(n)%rdat = fftw_alloc_real(alloc_local*2)

        call c_f_pointer(myplans(n)%cdat, myplans(n)%rin, [nvars, nlevs, nlats/2, 4, local_n0])
        call c_f_pointer(myplans(n)%cdat, myplans(n)%cout, [nvars, nlevs, nlats/2, 2, local_n0])

        flags = plan_flags
        myplans(n)%plan = fftw_mpi_plan_many_dft_r2c (rank, n0, howmany, nlons_local, &
                            nlons_local, myplans(n)%rin, myplans(n)%cout, comm_in, flags)
         
        register_plan = n
    end function register_plan

    subroutine fft3d(rinp, coutp)
        implicit none
        real, intent(in) :: rinp(:,:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:,:,:)
        integer :: id

        id = id3d

        myplans(id)%rin(1,:,:,1:2,:) = reshape(rinp,[NLEV,NLAT/2,2,NLON_LOCAL])*RSCALE
       
        call fftw_mpi_execute_dft_r2c(myplans(id)%plan, myplans(id)%rin, myplans(id)%cout) 

        coutp(:,:,:,:) = reshape(myplans(id)%cout(1,:,:,:,:),shape=[NLON_LOCAL,2,NLAT/2,NLEV],order=[4,3,2,1])

    end subroutine fft3d

    subroutine fft2d(rinp, coutp)
        implicit none
        real, intent(in) :: rinp(:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:,:)
        integer :: id

        id = id2d

        myplans(id)%rin(1,1,:,1:2,:) = reshape(rinp,[NLAT/2,2,NLON_LOCAL])*RSCALE
       
        call fftw_mpi_execute_dft_r2c(myplans(id)%plan, myplans(id)%rin, myplans(id)%cout) 

        coutp(:,:,:) = reshape(myplans(id)%cout(1,1,:,:,:),shape=[NLON_LOCAL,2,NLAT/2],order=[3,2,1])

    end subroutine fft2d

    subroutine end_fft_guru()
        implicit none
        integer :: i

        do i = 1, nplan
    
            call fftw_destroy_plan(myplans(i)%plan)
            call fftw_free(myplans(i)%cdat)
            call fftw_free(myplans(i)%rdat)
        enddo 
    end subroutine end_fft_guru

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
