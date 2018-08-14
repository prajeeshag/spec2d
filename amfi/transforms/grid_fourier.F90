module grid_fourier_mod

!--------------------------------------------------------------------------------   
! Module for grid to fourier and fourier to grid transformation. 
!

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_pe, mpp_npes, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                    mpp_sync, mpp_root_pe, CLOCK_MODULE, CLOCK_ROUTINE
use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, WARNING, NOTE, &
                    open_file, close_file, file_exist

implicit none

include 'fftw3-mpi.f03'

private

integer :: nplang2f=0, nplanf2g=0
integer, parameter :: max_plans=10
integer, parameter :: rank=2
integer(C_INTPTR_T), parameter :: TWO=2, ONE=1

type ocplan_type
    type(C_PTR) :: plan
    integer :: is, ie, ilen, rlat
    real :: rscale
    integer(C_INTPTR_T) :: howmany
    type(C_PTR) :: iptr, optr
    real(C_DOUBLE), pointer :: din(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: dout(:,:)
end type ocplan_type

type plan_type
    integer(C_INTPTR_T) :: howmany
    type(C_PTR) :: tplan1, tplan2, fplan
    type(C_PTR) :: t1dat, t2dat, r2c
    real(C_DOUBLE), pointer :: G(:,:), GT(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: FT(:,:), sFT(:,:)
    real(C_DOUBLE), pointer :: rsFT(:,:,:), rsF(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: sF(:,:)
endtype plan_type

integer(C_INTPTR_T) :: NLON, FTRUNC, NLON_LOCAL, FLOCAL=-1
integer(C_INTPTR_T) :: NLEV, NLAT, NVAR
integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
integer :: COMM_FFT
integer, allocatable :: Tshuffle(:)
logical :: shuffle=.false.

logical :: debug=.false.
integer :: clck_grid_to_fourier, clck_fourier_to_grid, clck_g2f_tran, clck_g2f_dft, &
           clck_f2g_tran, clck_f2g_dft, clck_plan_g2f, clck_plan_f2g

real :: RSCALE
integer :: plan_level = 3, plan_flags, id_grid2four, id_four2grid
logical :: initialized=.false.

character (len=256) :: wsdmfnm='fftw.wisdom'
character (len=256) :: null_plan_msg

type(plan_type) :: g2fp(max_plans), f2gp(max_plans)

public :: init_grid_fourier, end_grid_fourier, grid_to_fourier, fourier_to_grid

namelist/grid_fourier_nml/plan_level, debug

contains

!--------------------------------------------------------------------------------   
subroutine init_grid_fourier (nlons, ilen, nfourier, isf, flen, comm_in, Tshuff)
!--------------------------------------------------------------------------------   
    implicit none
    integer, intent(in) :: nlons ! Total no: of longitudes
    integer, intent(in) :: comm_in !MPI Communicator
    integer, intent(in) :: ilen ! No: of longitudes in this proc
    integer, intent(in) :: nfourier ! fourier truncation
    integer, intent(out) :: flen ! No: fouriers in this proc
    integer, intent(out) :: isf ! Starting of the fourier in this proc
    integer, intent(out), optional :: Tshuff(nfourier+1) !if present shuffle the fourier for 
                                                         !load balance for triangular truncation
                                                         !and give the order of shuffled fouriers
                                                         !in Tshuff
    integer :: unit, i, k, stat, n, flags, t
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2)
    integer(C_INTPTR_T) :: howmany

    unit = open_namelist_file()

    read(unit, nml=grid_fourier_nml,iostat=stat)

    call close_file(unit)

    NLON = nlons
    COMM_FFT = comm_in
    NLON_LOCAL = ilen
    FTRUNC = nfourier + 1 !because here everything starts from 1 not from Zero
    RSCALE = 1./nlons

    flen = 0
    isf = 0

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
        call mpp_error('init_grid_fourier','Wrong option for plan_level, &
        & set plan_level (accepted values are 0-3) in grid_fourier_nml', FATAL)
    end select

    if(mod(NLON,mpp_npes())/=0) &
        call mpp_error('grid_fourier_mod','No: of Pes in x-direction should be a factor of NLONS', FATAL)

    call fftw_mpi_init()

    clck_grid_to_fourier = mpp_clock_id('grid_to_fourier')
    clck_fourier_to_grid = mpp_clock_id('fourier_to_grid')
    clck_plan_g2f = mpp_clock_id('plan_grid_to_fourier')
    clck_plan_f2g = mpp_clock_id('plan_fourier_to_grid')
    clck_g2f_tran = mpp_clock_id('g2f_tran')
    clck_g2f_dft = mpp_clock_id('g2f_dft')
    clck_f2g_tran = mpp_clock_id('f2g_tran')
    clck_f2g_dft = mpp_clock_id('f2g_dft')

    if (import_wisdom(nlons)==1) then
        plan_flags = ior(plan_flags,FFTW_WISDOM_ONLY)
    endif

    null_plan_msg = 'NULL PLAN: try rerunning after removing the wisdom file '//trim(wsdmfnm)

    if (present(Tshuff)) then
        allocate(Tshuffle(FTRUNC))
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

        Tshuff = Tshuffle-1
        
        if (debug.and.mpp_pe()==mpp_root_pe()) then
            print *, 'Tshuffle=', Tshuff
        endif
    endif

    howmany = 1
    n0=[NLON,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    local_n1_prev = local_n1

    if (ilen/=local_n0) & 
        call mpp_error('init_grid_fourier', 'ilen/=local_n0', FATAL)

    n0=[howmany,FTRUNC]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n1_prev/=local_n0) &
        call mpp_error('init_grid_fourier', 'local_n1_prev/=local_n0', FATAL)

    isf = local_1_start
    flen = local_n1
    FLOCAL = local_n1

    initialized = .true.
    call mpp_error('init_grid_fourier', '----Initialized----', NOTE)

    return
end subroutine init_grid_fourier

!--------------------------------------------------------------------------------   
function plan_grid_to_fourier(howmany)
!--------------------------------------------------------------------------------   
    implicit none
    integer(C_INTPTR_T), intent(in) :: howmany
    integer :: plan_grid_to_fourier, n, flags, t, clck_transpose
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2)
    integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)

    call mpp_clock_begin(clck_plan_g2f)
    if (howmany<1) call mpp_error('plan_grid_to_fourier', 'howmany cannot be Zero', FATAL)

    nplang2f = nplang2f + 1
    
    if(nplang2f>max_plans) call mpp_error('grid_fourier_mod','No: g2f plans > Max_plans', FATAL)

    n = nplang2f
    plan_grid_to_fourier = n

    g2fp(n)%howmany = howmany
   
    !Transpose
    n0=[NLON,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (NLON_LOCAL/=local_n0) & 
        call mpp_error('plan_grid_to_fourier', 'ilen/=local_n0', FATAL)

    flags = plan_flags

    g2fp(n)%t1dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%G, [howmany, local_n0])
    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%GT, [NLON, local_n1])

    g2fp(n)%tplan1 = fftw_mpi_plan_many_transpose(NLON, howmany, 1, &
                            block0, block0, g2fp(n)%G, g2fp(n)%GT, &
                            COMM_FFT, flags)
    
    if (.not.c_associated(g2fp(n)%tplan1)) &
        call mpp_error('plan_grid_to_fourier: transpose1:',trim(null_plan_msg),FATAL)

    !multi-threaded shared memory fft
    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%FT, [(NLON/2+1),local_n1])
  
    nn(1) = NLON 
    idist = NLON; odist= NLON/2+1
    istride = 1; ostride = 1
    inembed = [NLON]; onembed = [NLON/2+1]
    flags = plan_flags

    g2fp(n)%fplan = fftw_plan_many_dft_r2c(ONE, nn, int(local_n1), &
                            g2fp(n)%GT, inembed, istride, idist, &
                            g2fp(n)%FT, onembed, ostride, odist, flags) 
    if (.not.c_associated(g2fp(n)%fplan)) &
        call mpp_error('plan_grid_to_fourier: dft_r2c:',trim(null_plan_msg),FATAL)
   
    local_n1_prev = local_n1
    !Transpose back

    n0=[howmany,FTRUNC]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n1_prev/=local_n0) &
        call mpp_error('plan_grid_to_fourier', 'local_n1_prev/=local_n0', FATAL)

    if (FLOCAL/=local_n1) call mpp_error('plan_grid_to_fourier', 'FLOCAL/=local_n1', FATAL)

    g2fp(n)%t2dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(g2fp(n)%t2dat, g2fp(n)%sFT, [FTRUNC,local_n0])
    call c_f_pointer(g2fp(n)%t2dat, g2fp(n)%rsFT, [TWO,FTRUNC,local_n0])
    call c_f_pointer(g2fp(n)%t2dat, g2fp(n)%rsF, [TWO,howmany,local_n1])

    g2fp(n)%tplan2 = fftw_mpi_plan_many_transpose(howmany, FTRUNC, 2, &
                            block0, block0, g2fp(n)%rsFT, g2fp(n)%rsF, &
                            COMM_FFT, flags) 

    if (.not.c_associated(g2fp(n)%tplan2)) &
        call mpp_error('plan_grid_to_fourier: transpose2:',trim(null_plan_msg),FATAL)

    g2fp(n)%r2c = c_loc(g2fp(n)%rsF)
    call c_f_pointer(g2fp(n)%r2c, g2fp(n)%sF, [howmany, local_n1])

    call save_wisdom()
    call mpp_clock_end(clck_plan_g2f)

    return

end function plan_grid_to_fourier

!--------------------------------------------------------------------------------   
function import_wisdom(nl)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: nl
    integer :: unitw
    integer :: isuccess, import_wisdom
    character(len=256) :: tmpc

    isuccess = 0

    write(tmpc,*) nl
    wsdmfnm = trim(wsdmfnm)//'_'//trim(adjustl(tmpc)) 

    write(tmpc,*) mpp_npes()
    wsdmfnm = trim(wsdmfnm)//'_'//trim(adjustl(tmpc)) 

    if (file_exist(trim(wsdmfnm))) then
        unitw = open_file(file=trim(wsdmfnm),action='read',form='ascii') 
        call import_wisdom_from_file(isuccess, unitw)
        if(isuccess==1) then
            call mpp_error(NOTE, 'grid_fourier_mod: wisdom import successful')
        else
            call mpp_error(NOTE, 'grid_fourier_mod: wisdom import failed')
        endif
        call close_file(unitw)
    else
        call mpp_error(NOTE, 'grid_fourier_mod: no wisdom file found')
    endif

    call mpp_sync()

    import_wisdom = isuccess
    return
end function import_wisdom 

!--------------------------------------------------------------------------------   
subroutine save_wisdom()
!--------------------------------------------------------------------------------   
    integer :: unitw

    call fftw_mpi_gather_wisdom(COMM_FFT)

    if (mpp_pe()==mpp_root_pe()) then
        unitw = open_file(file=trim(wsdmfnm),action='write',form='ascii') 
        call export_wisdom_to_file(unitw)
        call close_file(unitw)
    endif

    call mpp_sync()

    return
end subroutine save_wisdom 

!--------------------------------------------------------------------------------     
subroutine grid_to_fourier(Gp, sFp, id_in)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(in) :: Gp(:,:) ! lev, lat, lon
    complex, intent(out) :: sFp(:,:) ! fourier, lat, lev
    integer, intent(inout), optional :: id_in
    integer :: id, i, j, ci=1, cj=2, ct
    integer(C_INTPTR_T) :: howmany


    howmany = size(Gp,1)
    id = 0
    if (present(id_in)) id = id_in

    if (id<1) then
        do i = 1, nplang2f
            if (howmany==g2fp(i)%howmany) then
                id = i
                exit
            endif
        enddo
    endif
    
    if (id<1) then
        id = plan_grid_to_fourier(howmany)
    end if
                 
    call mpp_clock_begin(clck_grid_to_fourier)

    if (present(id_in)) id_in = id
    
    howmany = g2fp(id)%howmany

    g2fp(id)%G(1:howmany,:) = Gp(:,:)*RSCALE

    !Transpose
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fp(id)%tplan1, g2fp(id)%G, g2fp(id)%GT)
    call mpp_clock_end(clck_g2f_tran)

    !Serial FFT
    call mpp_clock_begin(clck_g2f_dft)
    call fftw_execute_dft_r2c(g2fp(id)%fplan, g2fp(id)%GT, g2fp(id)%FT) 
    call mpp_clock_end(clck_g2f_dft)

    !Truncation
    if (shuffle) then
        do i = 1, FTRUNC
            j = Tshuffle(i)
            g2fp(id)%sFT(i,:) = g2fp(id)%FT(j,:)
        enddo 
    else 
        g2fp(id)%sFT(1:FTRUNC,:) = g2fp(id)%FT(1:FTRUNC,:)
    endif

    !Transpose Back
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fp(id)%tplan2, g2fp(id)%rsFT, g2fp(id)%rsF)
    call mpp_clock_end(clck_g2f_tran)

    sFp = g2fp(id)%sF(1:howmany,1:FLOCAL)

    call mpp_clock_end(clck_grid_to_fourier)
            
end subroutine grid_to_fourier


!--------------------------------------------------------------------------------   
function plan_fourier_to_grid(howmany)
!--------------------------------------------------------------------------------   
    implicit none
    integer(C_INTPTR_T), intent(in) :: howmany
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
    
    if(nplanf2g>max_plans) call mpp_error('grid_fourier_mod','No: f2g plans > Max_plans', FATAL)

    n = nplanf2g
    plan_fourier_to_grid = n

    f2gp(n)%howmany = howmany
   
    !Transpose
    n0=[howmany,NLON]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (NLON_LOCAL/=local_n1) & 
        call mpp_error('plan_fourier_to_grid', 'ilen/=local_n1', FATAL)

    f2gp(n)%t1dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%GT, [NLON, local_n0])
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%G, [howmany, local_n1])

    flags = plan_flags

    f2gp(n)%tplan1 = fftw_mpi_plan_many_transpose(howmany, NLON, 1, &
                            block0, block0, f2gp(n)%GT, f2gp(n)%G, &
                            COMM_FFT, flags)
    if (.not.c_associated(f2gp(n)%tplan1)) &
        call mpp_error('plan_fourier_to_grid: transpose1:',trim(null_plan_msg),FATAL)

    !multi-threaded shared memory fft
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%GT, [NLON,local_n0])
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%FT, [(NLON/2+1),local_n0])
  
    nn(1) = NLON 
    idist = NLON; odist= NLON/2+1
    istride = 1; ostride = 1
    inembed = [NLON]; onembed = [NLON/2+1]
    flags = plan_flags

    f2gp(n)%fplan = fftw_plan_many_dft_c2r(ONE, nn, int(local_n0), &
                            f2gp(n)%FT, onembed, ostride, odist, & 
                            f2gp(n)%GT, inembed, istride, idist, flags)
    if (.not.c_associated(f2gp(n)%fplan)) &
        call mpp_error('plan_fourier_to_grid: dft_c2r:',trim(null_plan_msg),FATAL)
   
    local_n0_prev = local_n0 

    !Transpose back
    n0=[FTRUNC,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n0_prev /= local_n1) &
        call mpp_error('plan_fourier_to_grid', 'local_n0_prev /= local_n1', FATAL)

    if(FLOCAL/=local_n0) call mpp_error('plan_fourier_to_grid', 'FLOCAL/=local_n0', FATAL)

    f2gp(n)%t2dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(f2gp(n)%t2dat, f2gp(n)%sFT, [FTRUNC,local_n1])
    call c_f_pointer(f2gp(n)%t2dat, f2gp(n)%rsFT, [TWO,FTRUNC,local_n1])
    call c_f_pointer(f2gp(n)%t2dat, f2gp(n)%rsF, [TWO,howmany,local_n0])

    f2gp(n)%tplan2 = fftw_mpi_plan_many_transpose(FTRUNC, howmany, 2, &
                            block0, block0, f2gp(n)%rsF, f2gp(n)%rsFT, &
                            COMM_FFT, flags) 
    if (.not.c_associated(f2gp(n)%tplan2)) &
        call mpp_error('plan_fourier_to_grid: transpose2:',trim(null_plan_msg),FATAL)

    f2gp(n)%r2c = c_loc(f2gp(n)%rsF)
    call c_f_pointer(f2gp(n)%r2c, f2gp(n)%sF, [howmany, local_n0])

    call save_wisdom()

    call mpp_clock_end(clck_plan_f2g)

    return

end function plan_fourier_to_grid

!--------------------------------------------------------------------------------   
subroutine fourier_to_grid(sFp, Gp, id_in)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(out) :: Gp(:,:) ! howmany, lon
    complex, intent(in) :: sFp(:,:) ! howmany, fourier
    integer, intent(inout), optional :: id_in
    integer :: id, i, j, ci=1, cj=2, ct
    integer(C_INTPTR_T) :: howmany


    howmany = size(sFp,1)

    id = 0
    if (present(id_in)) id = id_in

    if (id<1) then
        do i = 1, nplanf2g
            if (howmany==f2gp(i)%howmany) then
                id = i
                exit
            endif
        enddo
    endif
    
    if (id<1) then
        id = plan_fourier_to_grid(howmany)
    endif
                 
    call mpp_clock_begin(clck_fourier_to_grid)

    if (present(id_in)) id_in = id
    
    howmany = f2gp(id)%howmany
    
    f2gp(id)%sF = sFp

    !Transpose Back
    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gp(id)%tplan2, f2gp(id)%rsF, f2gp(id)%rsFT)
    call mpp_clock_end(clck_f2g_tran)

    !Serial FFT
    f2gp(id)%FT(:,:) = 0.
    if (shuffle) then
        do i = 1, FTRUNC
            j = Tshuffle(i)
            f2gp(id)%FT(j,:) = f2gp(id)%sFT(i,:) !Truncation & Shuffle
        enddo
    else
        f2gp(id)%FT(1:FTRUNC,:) = f2gp(id)%sFT(1:FTRUNC,:) !Truncation
    endif

    call mpp_clock_begin(clck_f2g_dft)
    call fftw_execute_dft_c2r(f2gp(id)%fplan, f2gp(id)%FT, f2gp(id)%GT) 
    call mpp_clock_end(clck_f2g_dft)

    !Transpose
    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gp(id)%tplan1, f2gp(id)%GT, f2gp(id)%G)
    call mpp_clock_end(clck_f2g_tran)
    
    Gp = f2gp(id)%G(1:howmany,:)
            
    call mpp_clock_end(clck_fourier_to_grid)
end subroutine fourier_to_grid

!--------------------------------------------------------------------------------   
subroutine end_grid_fourier()
!--------------------------------------------------------------------------------   
    implicit none

    call fftw_mpi_cleanup()

end subroutine end_grid_fourier


!--------------------------------------------------------------------------------   
subroutine fft_1dr2c_serial(Gp, sFp)
!--------------------------------------------------------------------------------   
    real :: Gp(:)
    complex :: sFp(:)
    real(C_DOUBLE), pointer :: Gpl(:)
    complex(C_DOUBLE_COMPLEX), pointer :: sFpl(:)
    type(C_PTR) :: plan1d, data1d
    integer :: L
    L = size(Gp)

    data1d = fftw_alloc_complex(int(L/2+1, C_SIZE_T))
    call c_f_pointer(data1d, Gpl, [2*(L/2+1)])
    call c_f_pointer(data1d, sFpl, [L/2+1])

    plan1d = fftw_plan_dft_r2c_1d(L, Gpl, sFpl, FFTW_ESTIMATE) 

    if (.not.c_associated(plan1d)) &
        call mpp_error('fft_1dr2c_serial:','NULL PLAN',FATAL)

    Gpl = Gp
    call fftw_execute_dft_r2c(plan1d, Gpl, sFpl)
    sFp = sFpl
    
    call fftw_destroy_plan(plan1d)
    call fftw_free(data1d)

    return
end subroutine fft_1dr2c_serial


!--------------------------------------------------------------------------------   
subroutine fft_1dc2c_serial(sFp)
!--------------------------------------------------------------------------------   

    complex :: sFp(:)
    complex(C_DOUBLE_COMPLEX), pointer :: sFpl(:)
    type(C_PTR) :: plan1d, data1d
    integer :: L
    L = size(sFp)

    data1d = fftw_alloc_complex(int(L, C_SIZE_T))
    call c_f_pointer(data1d, sFpl, [L])

    plan1d = fftw_plan_dft_1d(L, sFpl, sFpl, FFTW_FORWARD, FFTW_ESTIMATE) 
    if (.not.c_associated(plan1d)) &
        call mpp_error('fft_1dc2c_serial:','NULL PLAN',FATAL)

    sFpl = sFp
    call fftw_execute_dft(plan1d, sFpl, sFpl)
    sFp = sFpl
    
    call fftw_destroy_plan(plan1d)
    call fftw_free(data1d)

end subroutine fft_1dc2c_serial

end module grid_fourier_mod 


#ifdef test_grid_fourier

program main

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_sum
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use grid_fourier_mod

implicit none

integer, allocatable :: pelist(:)
integer :: comm, nlon, ilen, trunc, isf, flen, nlat, howmany
integer :: j
real, allocatable :: grd(:,:,:)
real, allocatable :: grd1(:,:)
complex, allocatable :: four(:,:)

call mpp_init()

nlon = 10
ilen = nlon/mpp_npes()
trunc = nlon/2
nlat = 3
howmany = 6


allocate(pelist(mpp_npes()))
call mpp_get_current_pelist(pelist,commid=comm)
call init_grid_fourier (nlon, ilen, trunc, isf, flen, comm)

allocate(grd(howmany,nlat,ilen), four(howmany*nlat,flen), &
grd1(howmany*nlat,ilen))

do j = 1, nlat
    grd(:,j,:) = j
end do

grd1 = reshape(grd,[howmany*nlat,ilen])

call grid_to_fourier(grd1,four)

call mpp_exit()

end program main

#endif

