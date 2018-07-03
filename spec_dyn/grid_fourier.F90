module grid_fourier_mod
    
    use, intrinsic :: iso_c_binding

    use mpp_mod, only : mpp_pe, mpp_npes, mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
    use mpp_mod, only : mpp_root_pe
    use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, WARNING, NOTE

    implicit none

    include 'fftw3-mpi.f03'

    private

    integer :: nplang2f=0, nplanf2g=0
    integer, parameter :: max_plans=10
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

    integer(C_INTPTR_T) :: NLON, FTOTAL, FTRUNC, NLON_LOCAL, FLOCAL=-1
    integer(C_INTPTR_T) :: NLEV, NLAT, NVAR
    integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
    integer :: COMM_FFT
    integer, allocatable :: Tshuffle(:)
    logical :: shuffle=.false.

    logical :: debug=.false.
    integer :: clck_grid_to_fourier, clck_fourier_to_grid
    integer :: clck_plan_g2f, clck_plan_f2g

    real :: RSCALE
    integer :: plan_level = 3, plan_flags, id_grid2four, id_four2grid
    logical :: initialized=.false.


    character (len=16), parameter :: modul = 'grid_fourier_mod'


    type(plan_type) :: g2fplans(max_plans)
    type(plan_type) :: f2gplans(max_plans)

    public :: init_grid_fourier, end_grid_fourier
    public :: grid_to_fourier, fourier_to_grid
    public :: fft_1dr2c_serial, fft_1dc2c_serial


    namelist/grid_fourier_nml/plan_level, debug


    contains


    subroutine init_grid_fourier (nlons, ilen, nfourier, isf, flen, comm_in, Tshuff)

        implicit none

        integer, intent(in) :: nlons ! Total no: of longitudes
        integer, intent(in) :: comm_in !MPI Communicator
        integer, intent(in) :: ilen ! No: of longitudes in this proc
        integer, intent(in) :: nfourier ! fourier truncation
        integer, intent(inout) :: flen ! No: fouriers in this proc
        integer, intent(inout) :: isf ! Starting of the fourier in this proc
        integer, intent(out), optional :: Tshuff(nfourier+1) !if present shuffle the fourier for 
                                                             !load balance for triangular truncation
                                                             !and give the order of shuffled fouriers
                                                             !in Tshuff

                                                                    
        character (len=32) :: routine = 'init_grid_fourier'
        integer :: unit, i, k

        unit = open_namelist_file()

        read(unit, nml=grid_fourier_nml)

        call close_file(unit)

        NLON = nlons
        COMM_FFT = comm_in
        NLON_LOCAL = ilen
        FTOTAL = nlons/2
        FTRUNC = nfourier + 1 !because here everything starts from 1 not from Zero
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
            & set plan_level (accepted values are 0-3) in grid_fourier_nml', FATAL)
        end select

        if(mod(NLON,mpp_npes())/=0) &
            call mpp_error(modul,'No: of Pes in x-direction should be a factor of NLONS', FATAL)

        call fftw_mpi_init()

        clck_grid_to_fourier = mpp_clock_id('grid_to_fourier')
        clck_fourier_to_grid = mpp_clock_id('fourier_to_grid')
        clck_plan_g2f = mpp_clock_id('plan_grid_to_fourier')
        clck_plan_f2g = mpp_clock_id('plan_fourier_to_grid')

        flen = 0

        if (present(Tshuff)) then
            allocate(Tshuffle(FTRUNC))
            
            if (mpp_npes() > 1) then
                if (mod(FTRUNC,2)==0) then
                    k = FTRUNC + 1
                    do i = 1, FTRUNC/2
                        k = k - 1 
                        Tshuffle(2*(i-1)+1) = i
                        Tshuffle(2*i) = k
                    enddo
                else
                    Tshuffle(FTRUNC) = 1
                    k = FTRUNC + 1
                    do i = 1, (FTRUNC-1)/2
                        k = k - 1 
                        Tshuffle(2*(i-1)+1) = i + 1
                        Tshuffle(2*i) = k
                    enddo 
                endif
                shuffle=.true.
            else 
                shuffle=.false.
                forall(i=1:FTRUNC) Tshuffle(i) = i
            endif

            Tshuff = Tshuffle-1
            
            if (debug.and.mpp_pe()==mpp_root_pe()) then
                print *, 'Tshuffle=', Tshuff
            endif
        endif

        !grid_to_fourier
        id_grid2four = plan_grid_to_fourier(1,NLON_LOCAL,comm_in, isf, flen)

        !fourier_to_grid -> should be called after plan_grid_to_fourier
        id_four2grid = plan_fourier_to_grid(1,NLON_LOCAL,comm_in)

        initialized = .true.
        call mpp_error(routine, 'grid_to_fourier initialized !!!', NOTE)

    end subroutine init_grid_fourier


    function plan_grid_to_fourier(howmany, ilen, comm_in, isf, flen)

        implicit none
        integer(C_INTPTR_T), intent(in) :: howmany, ilen
        integer, intent(inout), optional :: isf, flen
        integer, intent(in) :: comm_in
        integer :: plan_grid_to_fourier, n, flags, t, clck_transpose
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
        integer(C_INTPTR_T) :: alloc_local, n0(2)
        integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

        call mpp_clock_begin(clck_plan_g2f)
        if (howmany<1) call mpp_error('plan_grid_to_fourier', 'howmany cannot be Zero', FATAL)

        nplang2f = nplang2f + 1
        
        if(nplang2f>max_plans) call mpp_error(modul,'No: g2f plans > Max_plans', FATAL)

        n = nplang2f
        plan_grid_to_fourier = n

        g2fplans(n)%howmany = howmany
       
        !Transpose
        n0=[NLON,howmany]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (ilen/=local_n0) & 
            call mpp_error('plan_grid_to_fourier', 'ilen/=local_n0', FATAL)

        flags = plan_flags

        g2fplans(n)%tcdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%rin, [howmany, local_n0])
        call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%trin, [NLON, local_n1])
    
        g2fplans(n)%plan = fftw_mpi_plan_many_transpose(NLON, howmany, 1, &
                                block0, block0, g2fplans(n)%rin, g2fplans(n)%trin, &
                                comm_in, flags)

        !multi-threaded shared memory fft
        call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%srin, [TWO*(NLON/2+1),local_n1])
        call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%scout, [(NLON/2+1),local_n1])
      
        nn(1) = NLON 
        idist = NLON; odist= NLON/2+1
        istride = 1; ostride = 1
        inembed = [NLON]; onembed = [NLON/2+1]
        flags = plan_flags

        g2fplans(n)%splan = fftw_plan_many_dft_r2c(ONE, nn, int(local_n1), &
                                g2fplans(n)%srin, inembed, istride, idist, &
                                g2fplans(n)%scout, onembed, ostride, odist, flags) 
       
        local_n1_prev = local_n1
        !Transpose back

        n0=[howmany,FTRUNC]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (local_n1_prev/=local_n0) &
            call mpp_error('plan_grid_to_fourier', 'local_n1_prev/=local_n0', FATAL)

        if(present(flen)) then
            flen = local_n1
            FLOCAL = local_n1
        else
            !FLOCAL should already have assinged and should be equal to local_n1
            if (FLOCAL/=local_n1) call mpp_error('plan_grid_to_fourier', 'FLOCAL/=local_n1', FATAL)
        endif 
 
        if(present(isf))isf = local_1_start

        g2fplans(n)%cdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%scouttr, [FTRUNC,local_n0])
        call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%srout, [TWO,FTRUNC,local_n0])
        call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%tsrout, [TWO,howmany,local_n1])

        g2fplans(n)%tplan = fftw_mpi_plan_many_transpose(howmany, FTRUNC, 2, &
                                block0, block0, g2fplans(n)%srout, g2fplans(n)%tsrout, &
                                comm_in, flags) 

        g2fplans(n)%r2c = c_loc(g2fplans(n)%tsrout)
        call c_f_pointer(g2fplans(n)%r2c, g2fplans(n)%cout, [howmany, local_n1])
        call mpp_clock_end(clck_plan_g2f)

    end function plan_grid_to_fourier

    subroutine grid_to_fourier(rinp, coutp, id_in)

        implicit none
        real, intent(in) :: rinp(:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:) ! fourier, lat, lev
        integer, intent(inout), optional :: id_in
        integer :: id, i, j, ci=1, cj=2, ct
        integer(C_INTPTR_T) :: howmany


        howmany = size(rinp,1)
        id = 0
        if (present(id_in)) id = id_in

        if (id<1) then
            do i = 1, nplang2f
                if (howmany==g2fplans(i)%howmany) then
                    id = i
                    exit
                endif
            enddo
        endif
        
        if (id<1) then
            id = plan_grid_to_fourier(howmany,NLON_LOCAL,COMM_FFT)
        endif
                     
        call mpp_clock_begin(clck_grid_to_fourier)

        if (present(id_in)) id_in = id
        
        howmany = g2fplans(id)%howmany

        !g2fplans(id)%rin(1:howmany,:) = reshape(rinp(:,:,:), shape=[howmany,NLON_LOCAL])*RSCALE
        g2fplans(id)%rin(1:howmany,:) = rinp(:,:)*RSCALE

        !Transpose
        call fftw_mpi_execute_r2r(g2fplans(id)%plan, g2fplans(id)%rin, g2fplans(id)%trin)
        
        !Serial FFT
        call fftw_execute_dft_r2c(g2fplans(id)%splan, g2fplans(id)%srin, g2fplans(id)%scout) 

        !Truncation
        if (shuffle) then
            do i = 1, FTRUNC
                j = Tshuffle(i)
                g2fplans(id)%scouttr(i,:) = g2fplans(id)%scout(j,:)
            enddo 
        else 
            g2fplans(id)%scouttr(1:FTRUNC,:) = g2fplans(id)%scout(1:FTRUNC,:)
        endif

        !Transpose Back
        call fftw_mpi_execute_r2r(g2fplans(id)%tplan, g2fplans(id)%srout, g2fplans(id)%tsrout)

        !coutp = reshape(g2fplans(id)%cout(1:howmany,1:FLOCAL),shape=[FLOCAL,NLAT,NLEV],order=[3,2,1])
        !coutp = reshape(g2fplans(id)%cout(1:howmany,1:FLOCAL),shape=[NLAT,NLEV,FLOCAL],order=[2,1,3])
        coutp = g2fplans(id)%cout(1:howmany,1:FLOCAL)

        call mpp_clock_end(clck_grid_to_fourier)
                
    end subroutine grid_to_fourier



    function plan_fourier_to_grid(howmany, ilen, comm_in)

        implicit none
        integer(C_INTPTR_T), intent(in) :: howmany, ilen
        integer, intent(in) :: comm_in
        integer :: plan_fourier_to_grid, n, flags, t, clck_transpose
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
        integer(C_INTPTR_T) :: local_n0_prev
        integer(C_INTPTR_T) :: alloc_local, n0(2)
        integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

        call mpp_clock_begin(clck_plan_f2g)

        if (howmany<1) call mpp_error('plan_fourier_to_grid', 'howmany cannot be Zero', FATAL)

        if (FLOCAL<0) &
            call mpp_error('plan_fourier_to_grid', 'plan_grid_to_fourier should be called first', FATAL)

        nplanf2g = nplanf2g + 1
        
        if(nplanf2g>max_plans) call mpp_error(modul,'No: f2g plans > Max_plans', FATAL)

        n = nplanf2g
        plan_fourier_to_grid = n

        f2gplans(n)%howmany = howmany
       
        !Transpose
        n0=[howmany,NLON]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (ilen/=local_n1) & 
            call mpp_error('plan_fourier_to_grid', 'ilen/=local_n1', FATAL)

        f2gplans(n)%tcdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%trin, [NLON, local_n0])
        call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%rin, [howmany, local_n1])
    
        flags = plan_flags

        f2gplans(n)%plan = fftw_mpi_plan_many_transpose(howmany, NLON, 1, &
                                block0, block0, f2gplans(n)%trin, f2gplans(n)%rin, &
                                comm_in, flags)

        !multi-threaded shared memory fft
        call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%srin, [TWO*(NLON/2+1),local_n0])
        call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%scout, [(NLON/2+1),local_n0])
      
        nn(1) = NLON 
        idist = NLON; odist= NLON/2+1
        istride = 1; ostride = 1
        inembed = [NLON]; onembed = [NLON/2+1]
        flags = plan_flags

        f2gplans(n)%splan = fftw_plan_many_dft_c2r(ONE, nn, int(local_n0), &
                                f2gplans(n)%scout, onembed, ostride, odist, & 
                                f2gplans(n)%srin, inembed, istride, idist, flags)
       
        local_n0_prev = local_n0 

        !Transpose back
        n0=[FTRUNC,howmany]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (local_n0_prev /= local_n1) &
            call mpp_error('plan_fourier_to_grid', 'local_n0_prev /= local_n1', FATAL)

        if(FLOCAL/=local_n0) call mpp_error('plan_fourier_to_grid', 'FLOCAL/=local_n0', FATAL)

        f2gplans(n)%cdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%scouttr, [FTRUNC,local_n1])
        call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%srout, [TWO,FTRUNC,local_n1])
        call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%tsrout, [TWO,howmany,local_n0])

        f2gplans(n)%tplan = fftw_mpi_plan_many_transpose(FTRUNC, howmany, 2, &
                                block0, block0, f2gplans(n)%tsrout, f2gplans(n)%srout, &
                                comm_in, flags) 

        f2gplans(n)%r2c = c_loc(f2gplans(n)%tsrout)
        call c_f_pointer(f2gplans(n)%r2c, f2gplans(n)%cout, [howmany, local_n0])

        call mpp_clock_end(clck_plan_f2g)

    end function plan_fourier_to_grid


    subroutine fourier_to_grid(coutp, rinp, id_in)

        implicit none
        real, intent(out) :: rinp(:,:) ! howmany, lon
        complex, intent(in) :: coutp(:,:) ! howmany, fourier
        integer, intent(inout), optional :: id_in
        integer :: id, i, j, ci=1, cj=2, ct
        integer(C_INTPTR_T) :: howmany


        howmany = size(coutp,1)

        id = 0
        if (present(id_in)) id = id_in

        if (id<1) then
            do i = 1, nplanf2g
                if (howmany==f2gplans(i)%howmany) then
                    id = i
                    exit
                endif
            enddo
        endif
        
        if (id<1) then
            id = plan_fourier_to_grid(howmany,NLON_LOCAL,COMM_FFT)
        endif
                     
        call mpp_clock_begin(clck_fourier_to_grid)

        if (present(id_in)) id_in = id
        
        howmany = f2gplans(id)%howmany
        
        !f2gplans(id)%cout = reshape(coutp,shape=[howmany,FLOCAL])
        f2gplans(id)%cout = coutp

        !Transpose Back
        call fftw_mpi_execute_r2r(f2gplans(id)%tplan, f2gplans(id)%tsrout, f2gplans(id)%srout)

        !Serial FFT
        f2gplans(id)%scout(:,:) = 0.
        if (shuffle) then
            do i = 1, FTRUNC
                j = Tshuffle(i)
                f2gplans(id)%scout(j,:) = f2gplans(id)%scouttr(i,:) !Truncation & Shuffle
            enddo
        else
            f2gplans(id)%scout(1:FTRUNC,:) = f2gplans(id)%scouttr(1:FTRUNC,:) !Truncation
        endif
        call fftw_execute_dft_c2r(f2gplans(id)%splan, f2gplans(id)%scout, f2gplans(id)%srin) 

        !Transpose
        call fftw_mpi_execute_r2r(f2gplans(id)%plan, f2gplans(id)%trin, f2gplans(id)%rin)
        
        !rinp = reshape(f2gplans(id)%rin(1:howmany,:), shape=[NLEV, NLAT, NLON_LOCAL], order=[2,1,3])
        rinp = f2gplans(id)%rin(1:howmany,:)
                
        call mpp_clock_end(clck_fourier_to_grid)
    end subroutine fourier_to_grid


    subroutine end_grid_fourier()
        implicit none

        call fftw_mpi_cleanup()

    end subroutine end_grid_fourier



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

