module grid_to_fourier_mod
    
    use, intrinsic :: iso_c_binding

    use mpp_mod, only : mpp_pe, mpp_npes, mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
    use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, WARNING, NOTE

    implicit none

    include 'fftw3-mpi.f03'

    private

    integer :: nplan=0
    integer, parameter :: max_plans=3
    integer, parameter :: rank=2
    integer(C_INTPTR_T), parameter :: TWO=2, ONE=1

    type plan_type
        integer(C_INTPTR_T) :: howmany
        type(C_PTR) :: plan, tplan, splan, tcdat, cdat, r2c
        real(C_DOUBLE), pointer :: rin(:,:), trin(:,:), srin(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: scout(:,:), scouttr(:,:)
        real(C_DOUBLE), pointer :: srout(:,:,:), tsrout(:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:)
    endtype plan_type

    integer(C_INTPTR_T) :: NLON, FTOTAL, FTRUNC, NLON_LOCAL, FLOCAL
    integer(C_INTPTR_T) :: NLEV, NLAT, NVAR
    integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
    integer :: COMM_FFT


    real :: RSCALE
    integer :: plan_level = 3, plan_flags, id_grid2four, id_four2grid
    logical :: initialized=.false.


    character (len=16), parameter :: modul = 'grid_to_fourier_mod'


    type(plan_type) :: myplans(max_plans)


    public :: init_grid_to_fourier, fft_1dr2c_serial, fft_1dc2c_serial, end_grid_to_fourier, grid_to_fourier


    namelist/grid_to_fourier_nml/plan_level


    contains


    subroutine init_grid_to_fourier (nlons, ilen, nfourier, isf, flen, comm_in, nlevs, nlats, nvars, tshuffle)

        implicit none

        integer, intent(in) :: nlons, comm_in, ilen, nfourier, nlats, nlevs
        integer, intent(inout) :: flen, isf
        integer, intent(in), optional :: nvars
        integer, intent(out), optional :: tshuffle(:) !if present shuffle the fourier for 
                                                      !load balance for triangular truncation
        character (len=32) :: routine = 'init_grid_to_fourier'
        integer :: unit

        unit = open_namelist_file()

        read(unit, nml=grid_to_fourier_nml)

        call close_file(unit)

        NLON = nlons
        NLAT = nlats
        NLEV = nlevs
        NVAR = 1
        if (present(nvars))NVAR=nvars

        COMM_FFT = comm_in
        NLON_LOCAL = ilen
        FTOTAL = nlons/2
        FTRUNC = nfourier 
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
            call mpp_error(routine,'Wrong option for plan_level, &
            & set plan_level (accepted values are 0-3) in grid_to_fourier_nml', FATAL)
        end select

        if(mod(NLON,mpp_npes())/=0) &
            call mpp_error(modul,'No: of Pes in x-direction should be a factor of NLONS', FATAL)

        call fftw_mpi_init()

        flen = 0

        !grid_to_fourier
        id_grid2four = plan_grid_to_fourier(1,NLEV,NLAT,NLON_LOCAL,comm_in, isf, flen)

        !fourier_to_grid -> should be called after plan_grid_to_fourier
        id_four2grid = plan_fourier_to_grid(1,NLEV,NLAT,NLON_LOCAL,comm_in, flen)

        initialized = .true.
        call mpp_error(routine, 'grid_to_fourier initialized !!!', NOTE)

    end subroutine init_grid_to_fourier


    function plan_grid_to_fourier(nvars, nlevs, nlats, ilen, comm_in, isf, flen)

        implicit none
        integer(C_INTPTR_T), intent(in) :: nvars, nlevs, nlats, ilen
        integer, intent(inout) :: isf, flen
        integer, intent(in) :: comm_in
        integer(C_INTPTR_T) :: howmany
        integer :: plan_grid_to_fourier, n, flags, t, clck_transpose
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
        integer(C_INTPTR_T) :: alloc_local, n0(2), oblock
        integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

        howmany = nvars*nlevs*nlats

        if (howmany<1) call mpp_error('plan_grid_to_fourier', 'howmany cannot be Zero', FATAL)

        nplan = nplan + 1
        
        if(nplan>max_plans) call mpp_error(modul,'No: plans > Max_plans', FATAL)

        n = nplan
        plan_grid_to_fourier = n

        myplans(n)%howmany = howmany
       
        !Transpose
        n0=[NLON,howmany]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (ilen/=local_n0) & 
            call mpp_error('plan_grid_to_fourier', 'ilen/=local_n0', FATAL)

        flags = plan_flags

        myplans(n)%tcdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(myplans(n)%tcdat, myplans(n)%rin, [howmany, local_n0])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%trin, [NLON, local_n1])
    
        myplans(n)%plan = fftw_mpi_plan_many_transpose(NLON, howmany, 1, &
                                block0, block0, myplans(n)%rin, myplans(n)%trin, &
                                comm_in, flags)

        !multi-threaded shared memory fft
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%srin, [TWO*(NLON/2+1),local_n1])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%scout, [(NLON/2+1),local_n1])
      
        nn(1) = NLON 
        idist = NLON; odist= NLON/2+1
        istride = 1; ostride = 1
        inembed = [NLON]; onembed = [NLON/2+1]
        flags = plan_flags

        myplans(n)%splan = fftw_plan_many_dft_r2c(ONE, nn, int(local_n1), &
                                myplans(n)%srin, inembed, istride, idist, &
                                myplans(n)%scout, onembed, ostride, odist, flags) 
       
        local_n1_prev = local_n1
        !Transpose back

        n0=[howmany,FTRUNC]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (local_n1_prev/=local_n0) &
            call mpp_error('plan_grid_to_fourier', 'local_n1_prev/=local_n0', FATAL)

        FLOCAL = local_n1
        flen = local_n1
        isf = local_1_start

        myplans(n)%cdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(myplans(n)%cdat, myplans(n)%scouttr, [FTRUNC,local_n0])
        call c_f_pointer(myplans(n)%cdat, myplans(n)%srout, [TWO,FTRUNC,local_n0])
        call c_f_pointer(myplans(n)%cdat, myplans(n)%tsrout, [TWO,howmany,local_n1])

        myplans(n)%tplan = fftw_mpi_plan_many_transpose(howmany, FTRUNC, 2, &
                                block0, block0, myplans(n)%srout, myplans(n)%tsrout, &
                                comm_in, flags) 

        myplans(n)%r2c = c_loc(myplans(n)%tsrout)
        call c_f_pointer(myplans(n)%r2c, myplans(n)%cout, [howmany, local_n1])

    end function plan_grid_to_fourier

    subroutine grid_to_fourier(rinp, coutp)

        implicit none
        real, intent(in) :: rinp(:,:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:,:) ! fourier, lat, lev
        integer :: id, i, j, ci=1, cj=2, ct
        integer(C_INTPTR_T) :: howmany

        id = id_grid2four
        
        howmany = myplans(id)%howmany

        myplans(id)%rin(1:howmany,:) = reshape(rinp(:,:,:), shape=[howmany,NLON_LOCAL])*RSCALE

        !Transpose
        call fftw_mpi_execute_r2r(myplans(id)%plan, myplans(id)%rin, myplans(id)%trin)
        
        !Serial FFT
        call fftw_execute_dft_r2c(myplans(id)%splan, myplans(id)%srin, myplans(id)%scout) 

        !Truncation        
        myplans(id)%scouttr(1:FTRUNC,:) = myplans(id)%scout(1:FTRUNC,:)

        !Transpose Back
        call fftw_mpi_execute_r2r(myplans(id)%tplan, myplans(id)%srout, myplans(id)%tsrout)

        coutp = reshape(myplans(id)%cout(1:howmany,1:FLOCAL),shape=[FLOCAL,NLAT,NLEV],order=[3,2,1])
                
    end subroutine grid_to_fourier



    function plan_fourier_to_grid(nvars, nlevs, nlats, ilen, comm_in, flen)

        implicit none
        integer(C_INTPTR_T), intent(in) :: nvars, nlevs, nlats, ilen
        integer, intent(in) :: flen
        integer, intent(in) :: comm_in
        integer(C_INTPTR_T) :: howmany
        integer :: plan_fourier_to_grid, n, flags, t, clck_transpose
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
        integer(C_INTPTR_T) :: local_n0_prev
        integer(C_INTPTR_T) :: alloc_local, n0(2), oblock
        integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

        howmany = nvars*nlevs*nlats

        if (howmany<1) call mpp_error('plan_fourier_to_grid', 'howmany cannot be Zero', FATAL)

        if (flen<0) &
            call mpp_error('plan_fourier_to_grid', 'plan_grid_to_fourier should be called first', FATAL)

        nplan = nplan + 1
        
        if(nplan>max_plans) call mpp_error(modul,'No: plans > Max_plans', FATAL)

        n = nplan
        plan_fourier_to_grid = n

        myplans(n)%howmany = howmany
       
        !Transpose
        n0=[howmany,NLON]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (ilen/=local_n1) & 
            call mpp_error('plan_fourier_to_grid', 'ilen/=local_n1', FATAL)

        myplans(n)%tcdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(myplans(n)%tcdat, myplans(n)%trin, [NLON, local_n0])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%rin, [howmany, local_n1])
    
        flags = plan_flags

        myplans(n)%plan = fftw_mpi_plan_many_transpose(howmany, NLON, 1, &
                                block0, block0, myplans(n)%trin, myplans(n)%rin, &
                                comm_in, flags)

        !multi-threaded shared memory fft
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%srin, [TWO*(NLON/2+1),local_n0])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%scout, [(NLON/2+1),local_n0])
      
        nn(1) = NLON 
        idist = NLON; odist= NLON/2+1
        istride = 1; ostride = 1
        inembed = [NLON]; onembed = [NLON/2+1]
        flags = plan_flags

        myplans(n)%splan = fftw_plan_many_dft_c2r(ONE, nn, int(local_n0), &
                                myplans(n)%scout, onembed, ostride, odist, & 
                                myplans(n)%srin, inembed, istride, idist, flags)
       
        local_n0_prev = local_n0 
        !Transpose back

        n0=[FTRUNC,howmany]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (local_n0_prev /= local_n1) &
            call mpp_error('plan_fourier_to_grid', 'local_n0_prev /= local_n1', FATAL)

        if(flen/=local_n0) call mpp_error('plan_fourier_to_grid', 'flen/=local_n0', FATAL)

        call c_f_pointer(myplans(n)%tcdat, myplans(n)%srout, [TWO,NLON/2+1,local_n1])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%tsrout, [TWO,howmany,local_n0])

        myplans(n)%tplan = fftw_mpi_plan_many_transpose(FTRUNC, howmany, 2, &
                                block0, block0, myplans(n)%tsrout, myplans(n)%srout, &
                                comm_in, flags) 

        myplans(n)%r2c = c_loc(myplans(n)%tsrout)
        call c_f_pointer(myplans(n)%r2c, myplans(n)%cout, [howmany, local_n0])

    end function plan_fourier_to_grid

    subroutine fourier_to_grid(coutp, rinp)

        implicit none
        real, intent(out) :: rinp(:,:,:) ! lev, lat, lon
        complex, intent(in) :: coutp(:,:,:) ! fourier, lat, lev
        integer :: id, i, j, ci=1, cj=2, ct
        integer(C_INTPTR_T) :: howmany

        id = id_four2grid
        
        howmany = myplans(id)%howmany

        myplans(id)%cout = reshape(coutp(:,:,:),shape=[howmany,FLOCAL],order=[2,1])

        !Transpose Back
        call fftw_mpi_execute_r2r(myplans(id)%tplan, myplans(id)%tsrout, myplans(id)%srout)

        !Serial FFT
        myplans(id)%scout(FTRUNC+1:NLON/2+1,:) = 0. !Truncation
        call fftw_execute_dft_c2r(myplans(id)%splan, myplans(id)%scout, myplans(id)%srin) 

        !Transpose
        call fftw_mpi_execute_r2r(myplans(id)%plan, myplans(id)%trin, myplans(id)%rin)
        
        rinp = reshape(myplans(id)%rin(1:howmany,:), shape=[NLEV, NLAT, NLON_LOCAL])
                
    end subroutine fourier_to_grid



    subroutine end_grid_to_fourier()
        implicit none

        call fftw_mpi_cleanup()

    end subroutine end_grid_to_fourier



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

