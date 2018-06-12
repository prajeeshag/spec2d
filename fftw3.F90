module fft_guru
    
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
        type(C_PTR) :: plan, tplan, splan, tcdat, r2c
        real(C_DOUBLE), pointer :: rin(:,:), trin(:,:), srin(:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: scout(:,:)
        real(C_DOUBLE), pointer :: srout(:,:,:), tsrout(:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:)
    endtype plan_type

    integer(C_INTPTR_T) :: NLON, FTOTAL, FTRUNC, NLON_LOCAL, FLOCAL
    integer(C_INTPTR_T) :: NLEV, NLAT, NVAR
    integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
    integer :: COMM_FFT

    real :: RSCALE
    integer :: plan_level = 1, plan_flags, id2d, id3d, id3dext, idtrans
    logical :: transpos=.true.
    logical :: initialized=.false.

    character (len=16), parameter :: modul = 'fft_guru'

    type(plan_type) :: myplans(max_plans)

    public :: init_fft_guru, fft_1dr2c_serial, fft_1dc2c_serial, end_fft_guru, fft_trans

    namelist/fft_guru_nml/plan_level, transpos

    contains

    subroutine init_fft_guru (nlons, nlons_local, nfourier, fourier_start_local, &
                                nfourier_local, comm_in, nlevs, nlats, nvars)
        implicit none

        integer, intent(in) :: nlons, comm_in, nlons_local, nfourier, nlats, nlevs
        integer, intent(inout) :: nfourier_local, fourier_start_local
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
            call mpp_error(routine,'Wrong option for plan_level, &
            & set plan_level (accepted values are 0-3) in fft_guru_nml', FATAL)
        end select

        if(mod(NLON,2)) call mpp_error(routine,'NLON must be an integer',FATAL)

        call fftw_mpi_init()


        !transpose plan 
        idtrans = register_plan(1,NLEV,NLAT,NLON_LOCAL,comm_in, fourier_start_local, nfourier_local)

        nfourier_local = FLOCAL 
        initialized = .true.

        call mpp_error(routine, 'FFT_GURU initialized !!!', NOTE)

    end subroutine init_fft_guru

    function register_plan(nvars, nlevs, nlats, nlons_local, comm_in, fourier_start_local, nfourier_local)

        implicit none
        integer(C_INTPTR_T), intent(in) :: nvars, nlevs, nlats, nlons_local
        integer, intent(inout) :: fourier_start_local, nfourier_local
        integer, intent(in) :: comm_in
        integer(C_INTPTR_T) :: howmany
        integer :: register_plan, n, flags, t, clck_transpose
        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
        integer(C_INTPTR_T) :: alloc_local, n0(2), oblock
        integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

        howmany = nvars*nlevs*nlats

        if (howmany<1) call mpp_error('register_plan', 'howmany cannot be Zero', FATAL)

        nplan = nplan + 1
        
        if(nplan>max_plans) call mpp_error(modul,'No: plans > Max_plans', FATAL)

        n = nplan
        register_plan = n

        myplans(n)%howmany = howmany
       
        !Transpose
        n0=[NLON,howmany]

        if(mod(NLON,mpp_npes())/=0) &
            call mpp_error(modul,'No: of Pes in x-direction should be a factor of NLONS', FATAL)

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                       nlons_local, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        if (nlons_local/=local_n0) & 
            call mpp_error('register_plan', 'nlons_local/=local_n0', FATAL)

        flags = plan_flags

        myplans(n)%tcdat = fftw_alloc_complex(alloc_local)

        call c_f_pointer(myplans(n)%tcdat, myplans(n)%rin, [howmany, local_n0])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%trin, [NLON, local_n1])
    
        myplans(n)%plan = fftw_mpi_plan_many_transpose(NLON, howmany, 1, &
                                nlons_local, block0, myplans(n)%rin, myplans(n)%trin, &
                                comm_in, flags)

        !multi-threaded shared memory fft
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%srin, [TWO*(NLON/2+1),local_n1])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%scout, [(NLON/2+1),local_n1])
      
        nn(1) = NLON 
        idist = NLON; odist= FTRUNC
        istride = 1; ostride = 1
        inembed = [NLON]; onembed = [NLON/2+1]
        flags = plan_flags

        myplans(n)%splan = fftw_plan_many_dft_r2c(ONE, nn, int(local_n1), &
                                myplans(n)%srin, inembed, istride, idist, &
                                myplans(n)%scout, onembed, ostride, odist, flags) 
        
        !Transpose back

        n0=[howmany,FTRUNC]

        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                       block0, block0, comm_in, local_n0, &
                       local_0_start, local_n1, local_1_start)

        FLOCAL = local_n1
        nfourier_local = local_n1
        fourier_start_local = local_1_start

        call c_f_pointer(myplans(n)%tcdat, myplans(n)%srout, [TWO,FTRUNC,local_n0])
        call c_f_pointer(myplans(n)%tcdat, myplans(n)%tsrout, [TWO,howmany,FLOCAL])

        myplans(n)%tplan = fftw_mpi_plan_many_transpose(howmany, FTRUNC, 2, &
                                block0, block0, myplans(n)%srout, myplans(n)%tsrout, &
                                comm_in, flags) 

        myplans(n)%r2c = c_loc(myplans(n)%tsrout)
        call c_f_pointer(myplans(n)%r2c, myplans(n)%cout, [howmany, FLOCAL])

    end function register_plan


    subroutine fft_trans(rinp, coutp)
        implicit none
        real, intent(in) :: rinp(:,:,:) ! lev, lat, lon
        complex, intent(out) :: coutp(:,:,:) ! fourier, lat, lev
        integer :: id, i, j, ci=1, cj=2, ct
        integer(C_INTPTR_T) :: howmany

        id = idtrans
        
        howmany = myplans(id)%howmany

        myplans(id)%rin(1:howmany,:) = reshape(rinp(:,:,:), shape=[howmany,NLON_LOCAL])*RSCALE

        !Transpose
        call fftw_mpi_execute_r2r(myplans(id)%plan, myplans(id)%rin, myplans(id)%trin)
        
        !Serial FFT
        call fftw_execute_dft_r2c(myplans(id)%splan, myplans(id)%srin, myplans(id)%scout) 

        !Transpose Back
        call fftw_mpi_execute_r2r(myplans(id)%tplan, myplans(id)%srout, myplans(id)%tsrout)

        coutp = reshape(myplans(id)%cout(1:howmany,1:FLOCAL),shape=[NLEV,NLAT,FLOCAL])

    end subroutine fft_trans


    subroutine end_fft_guru()
        implicit none

        call fftw_mpi_cleanup()

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

!    function register_plan(nvars, nlevs, nlats, nlons_local, comm_in, fourier_start_local, nfourier_local)
!        implicit none
!        integer(C_INTPTR_T), intent(in) :: nvars, nlevs, nlats, nlons_local
!        integer, intent(in) :: comm_in
!        integer, intent(out) :: fourier_start_local, nfourier_local
!        integer(C_INTPTR_T) :: howmany
!        integer :: register_plan, n, flags
!        integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
!        integer(C_INTPTR_T) :: alloc_local, n0(2), oblock
!
!        howmany = nvars*nlevs*nlats/2
!
!        if (howmany<1) call mpp_error('register_plan', 'howmany cannot be Zero', FATAL)
!
!        nplan = nplan + 1
!        
!        if(nplan>max_plans) call mpp_error(modul,'No: plans > Max_plans', FATAL)
!
!        n = nplan
!        register_plan = n
!
!        myplans(n)%howmany = howmany
!       
!        n0=[NLON,TWO]
!
!        if(mpp_npes()==1) transpos = .false.
!
!        if (.not.transpos) then
!            alloc_local = fftw_mpi_local_size_many(rank, n0, howmany, &
!                           nlons_local, comm_in, local_n0, local_0_start)
!            local_n1 = local_n0
!            local_1_start = local_0_start
!            flags = plan_flags
!            oblock = nlons_local
!            myplans(n)%tn0 = local_n1
!
!            myplans(n)%cdat = fftw_alloc_complex(alloc_local)
!
!            call c_f_pointer(myplans(n)%cdat, myplans(n)%rin, [nvars, nlevs, nlats/2, 4, local_n0])
!            call c_f_pointer(myplans(n)%cdat, myplans(n)%cout, [nvars, nlevs, nlats/2, 2, local_n1])
!
!            myplans(n)%plan = fftw_mpi_plan_many_dft_r2c (rank, n0, howmany, nlons_local, &
!                                oblock, myplans(n)%rin, myplans(n)%cout, comm_in, flags)
!
!            return
!        endif
!
!        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, howmany, &
!                       nlons_local, block0, comm_in, local_n0, &
!                       local_0_start, local_n1, local_1_start)
!
!        flags = ior(plan_flags,FFTW_MPI_TRANSPOSED_OUT)
!        oblock = block0
!        myplans(n)%tn0 = local_n1
!
!        myplans(n)%cdat = fftw_alloc_complex(alloc_local*2)
!
!        call c_f_pointer(myplans(n)%cdat, myplans(n)%rin, [nvars, nlevs, nlats/2, TWO*2, local_n0])
!        call c_f_pointer(myplans(n)%cdat, myplans(n)%tcout, [nvars, nlevs, nlats/2, NLON, local_n1])
!
!        myplans(n)%plan = fftw_mpi_plan_many_dft_r2c (rank, n0, howmany, nlons_local, &
!                                oblock, myplans(n)%rin, myplans(n)%tcout, comm_in, flags)
!
!        !Transpose
!        n0=[TWO,FTRUNC]
!
!        alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2*howmany, &
!                       myplans(n)%tn0, FLOCAL, comm_in, local_n0, &
!                       local_0_start, local_n1, local_1_start)
!
!!        if (FLOCAL/=local_n1) &
!!            call mpp_error('register_plan', 'Asked ', WARNING)
!
!        FLOCAL = local_n1
!        nfourier_local = local_n1
!        fourier_start_local = local_1_start
!
!        if (myplans(n)%tn0/=local_n0) &
!            call mpp_error('register_plan', 'myplans(n)%tn0/=local_n0, try a different no: &
!                                               pes on logitudinal direction',FATAL)
!
!        flags = plan_flags
!
!        myplans(n)%tcdat = fftw_alloc_complex(alloc_local*2)
!        !myplans(n)%tcdato = fftw_alloc_complex(alloc_local*2)
!        if (myplans(n)%tn0>1) call mpp_error('register_plan', 'tn0 problem...', fatal)
!        call c_f_pointer(myplans(n)%tcdat, myplans(n)%tcin, [nvars, nlevs, nlats/2, FTRUNC, myplans(n)%tn0])
!        call c_f_pointer(myplans(n)%tcdat, myplans(n)%trin, [TWO, nvars, nlevs, nlats/2, FTRUNC, myplans(n)%tn0])
!        call c_f_pointer(myplans(n)%tcdat, myplans(n)%trout, [TWO, nvars, nlevs, nlats/2, TWO, FLOCAL])
!    
!        myplans(n)%tplan = fftw_mpi_plan_many_transpose(TWO, FTRUNC, howmany*2, &
!                                myplans(n)%tn0, FLOCAL, myplans(n)%trin, myplans(n)%trout, &
!                                comm_in, flags) 
!
!        !!trout -> cout
!        !allocate(myplans(n)%cout(nvars, nlevs, nlats/2, 2, FLOCAL))
!        myplans(n)%r2c = c_loc(myplans(n)%trout)
!        call c_f_pointer(myplans(n)%r2c, myplans(n)%cout, [nvars, nlevs, nlats/2, TWO, FLOCAL])
!         
!    end function register_plan

!    subroutine fft3d(rinp, coutp)
!        implicit none
!        real, intent(in) :: rinp(:,:,:) ! lev, lat, lon
!        complex, intent(out) :: coutp(:,:,:,:)
!        integer :: id, j
!
!        id = id3d
!
!        myplans(id)%rin(1,:,:,1:2,:) = reshape(rinp,[NLEV,NLAT/2,2,NLON_LOCAL])*RSCALE
!      
!        if (.not.transpos) then 
!            call fftw_mpi_execute_dft_r2c(myplans(id)%plan, myplans(id)%rin, myplans(id)%cout) 
!        else
!            call fftw_mpi_execute_dft_r2c(myplans(id)%plan, myplans(id)%rin, myplans(id)%tcout)
!            if(myplans(id)%tn0>0) then
!                myplans(id)%tcin(:, :, :, 1:FTRUNC, :) = myplans(id)%tcout(:, :, :, 1:FTRUNC, :)
!            !    do j = 1, FTRUNC
!            !        print *, 'trin=',mpp_pe(), j, myplans(id)%trin(1,1,1,1,j,:), myplans(id)%trin(2,1,1,1,j,:)
!            !    enddo
!            endif
!            call fftw_mpi_execute_r2r(myplans(id)%tplan, myplans(id)%trin, myplans(id)%trout)
!            !print *, 'trout=',mpp_pe(), myplans(id)%trout(1,1,1,1,:,1:FLOCAL)
!        endif
!
!        coutp = reshape(myplans(id)%cout(1,:,:,:,1:FLOCAL), &
!                        shape=[FLOCAL,2,NLAT/2,NLEV],order=[4,3,2,1])
!
!    end subroutine fft3d

