module grid_fourier_mod
    
use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_pe, mpp_npes, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                    mpp_sync, mpp_root_pe, CLOCK_MODULE, CLOCK_ROUTINE
use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, WARNING, NOTE, &
                    open_file, close_file, file_exist
use ocpack_mod, only : ocpack_type 
implicit none

include 'fftw3-mpi.f03'

private

integer :: nplang2f=0, nplanf2g=0
integer, parameter :: max_plans=10
integer, parameter :: rank=2
integer(C_INTPTR_T), parameter :: TWO=2, ONE=1

type ocplan_type
    type(C_PTR) :: plan
    integer(C_INTPTR_T) :: howmany
    integer :: is, ie, ilen
end type ocplan_type

type plan_type
    integer(C_INTPTR_T) :: howmany
    type(C_PTR) :: plan, tplan, splan, tcdat, cdat, r2c
    real(C_DOUBLE), pointer :: rin(:,:), trin(:,:), srin(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: scout(:,:), scouttr(:,:)
    real(C_DOUBLE), pointer :: srout(:,:,:), tsrout(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: cout(:,:)
    type(ocplan_type), allocatable :: oc(:,:)
endtype plan_type

integer(C_INTPTR_T) :: NX, FTRUNC, NX_LOCAL, FLOCAL=-1
integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
integer :: COMM_FFT
integer :: js, je, jlen
type(ocpack_type), allocatable :: ocP(:,:)

integer, allocatable :: Tshuffle(:)
logical :: shuffle=.false.

logical :: debug=.false.
integer :: clck_grid_to_fourier, clck_fourier_to_grid, clck_g2f_tran, clck_g2f_dft, &
           clck_f2g_tran, clck_f2g_dft, clck_plan_g2f, clck_plan_f2g

integer :: plan_level = 3, plan_flags
logical :: initialized=.false.

character (len=256) :: wsdmfnm='fftw.wisdom'
character (len=256) :: null_plan_msg

type(plan_type) :: g2fplans(max_plans)
type(plan_type) :: f2gplans(max_plans)

public :: init_grid_fourier, end_grid_fourier, grid_to_fourier, fourier_to_grid, &
          plan_grid_to_fourier, plan_fourier_to_grid
!public :: fft_1dr2c_serial, fft_1dc2c_serial

namelist/grid_fourier_nml/plan_level, debug

contains

!--------------------------------------------------------------------------------   
subroutine init_grid_fourier (oc, jsp, jep, ilen, nfourier, isf, flen, comm_in, Tshuff)
!--------------------------------------------------------------------------------   
    implicit none
    type(ocpack_type), intent(in) :: oc(:,:)
    integer, intent(in) :: jsp, jep
    integer, intent(in) :: ilen ! No: of longitudes in this proc
    integer, intent(in) :: nfourier ! fourier truncation
    integer, intent(out) :: isf ! Starting of the fourier in this proc
    integer, intent(out) :: flen ! No: fouriers in this proc
    integer, intent(in) :: comm_in !MPI Communicator
    integer, intent(out), optional :: Tshuff(nfourier+1) !if present shuffle the fourier for 
    !load balance for triangular truncation and give the order of shuffled fouriers in Tshuff.

    integer :: unit, i, k, stat
    integer :: n, flags, t
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2)
    integer(C_INTPTR_T) :: howmany

    unit = open_namelist_file()
    read(unit, nml=grid_fourier_nml,iostat=stat)
    call close_file(unit)

    NX = oc(npack(),1)%ie
    COMM_FFT = comm_in
    NX_LOCAL = ilen
    FTRUNC = nfourier + 1 !because here everything starts from 1 not from Zero
    js = jsp
    je = jep
    jlen = jep - jsp + 1

    allocate(ocP(size(oc,1),size(oc,2)))
    ocP = oc
 
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

    if(mod(NX,mpp_npes())/=0) &
        call mpp_error(grid_fourier_mod,'No: of Pes is not a factor of NX', FATAL)

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
    n0=[NX,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, comm_in, local_n0, &
                   local_0_start, local_n1, local_1_start)

    local_n1_prev = local_n1

    if (ilen/=local_n0) & 
        call mpp_error('init_grid_fourier', 'ilen/=local_n0', FATAL)

    n0=[howmany,FTRUNC]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 2, &
                   block0, block0, comm_in, local_n0, &
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

    if (mod(howmany,jlen)/=0) &
        call mpp_error('plan_grid_to_fourier', 'howmany should be a factor of jlen', FATAL)

    nplang2f = nplang2f + 1
    if(nplang2f>max_plans) call mpp_error(grid_fourier_mod,'No: g2f plans > Max_plans', FATAL)

    n = nplang2f

    plan_grid_to_fourier = n
    g2fplans(n)%howmany = howmany
   
    !Transpose
    n0=[NX,howmany]
    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (NX_LOCAL/=local_n0) & 
        call mpp_error('plan_grid_to_fourier', 'NX_LOCAL/=local_n0', FATAL)

    flags = plan_flags
    g2fplans(n)%tcdat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%rin, [howmany, local_n0])
    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%trin, [NX, local_n1])

    g2fplans(n)%plan = fftw_mpi_plan_many_transpose(NX, howmany, 1, &
                            block0, block0, g2fplans(n)%rin, g2fplans(n)%trin, &
                            COMM_FFT, flags)
    
    if (.not.c_associated(g2fplans(n)%plan)) &
        call mpp_error('plan_grid_to_fourier: transpose1:',trim(null_plan_msg),FATAL)

    !FFT
    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%srin, [TWO*(NX/2+1),local_n1])
    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%scout, [(NX/2+1),local_n1])
   
    n_oc_plans = num_plans_for(local_n1, local_1_start, howmany) 
    allocate(g2fplans(n)%oc(size(ocP,1),n_oc_plans))
    call set_ocplan_parm(local_n1, local_1_start, howmany, g2fplans(n)%oc)

    nn(1) = g2fplans(n)%oc(i,j)%ilen
    idist = NX; odist= NX/2+1
    istride = 1; ostride = 1
    inembed = [g2fplans(n)%oc(i,j)%ilen]; 
    onembed = [g2fplans(n)%oc(i,j)%ilen/2+1]
    flags = plan_flags
    howmany1 = int(local_n1)

    g2fplans(n)%splan = fftw_plan_many_dft_r2c(ONE, nn, howmany1, &
                            g2fplans(n)%srin, inembed, istride, idist, &
                            g2fplans(n)%scout, onembed, ostride, odist, flags) 
    if (.not.c_associated(g2fplans(n)%splan)) &
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

    g2fplans(n)%cdat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%scouttr, [FTRUNC,local_n0])
    call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%srout, [TWO,FTRUNC,local_n0])
    call c_f_pointer(g2fplans(n)%cdat, g2fplans(n)%tsrout, [TWO,howmany,local_n1])

    g2fplans(n)%tplan = fftw_mpi_plan_many_transpose(howmany, FTRUNC, 2, &
                            block0, block0, g2fplans(n)%srout, g2fplans(n)%tsrout, &
                            COMM_FFT, flags) 

    if (.not.c_associated(g2fplans(n)%tplan)) &
        call mpp_error('plan_grid_to_fourier: transpose2:',trim(null_plan_msg),FATAL)

    g2fplans(n)%r2c = c_loc(g2fplans(n)%tsrout)
    call c_f_pointer(g2fplans(n)%r2c, g2fplans(n)%cout, [howmany, local_n1])

    call save_wisdom()
    call mpp_clock_end(clck_plan_g2f)

    return

end function plan_grid_to_fourier

!--------------------------------------------------------------------------------   
function num_plans_for(hmlen, shm, howmany)
    integer, intent(in) :: hmlen, shm, howmany
    integer :: num_plans_for
    integer :: hms, hme, numj, j, h

    hms = shm + 1 ! shmn is start of homany from fftw, which starts from 0.
    hme = hms + hmlen - 1
    numj = howmany/jlen

    num_plans_for = 0
    h = hms
    do while(h <= hme)
        rem = numj - mod(j,numj)
        h = h + rem + 1
        num_plans_for = num_plans_for+1
    end do
    
    return
end function num_plans_for
          

!--------------------------------------------------------------------------------   
subroutine set_ocpack_parm(hmlen, shmn, howmany, oc)
    integer, intent(in) :: hmlen, shmn, howmany
    type(ocplan_type), intent(out) :: oc(:,:)

    integer :: hms, hme, numj, j, h, n, howmany1

    if (size(oc,1)/=size(ocP,1)) call mpp_error(FATAL,'set_ocpack_for: size(oc,1)/=size(ocP,1)')

    hms = shmn + 1 ! shmn is start of homany from fftw, which starts from 0.
    hme = hms + hmlen - 1
    numj = howmany/jlen

    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%srin, [TWO*(NX/2+1),local_n1])
    call c_f_pointer(g2fplans(n)%tcdat, g2fplans(n)%scout, [(NX/2+1),local_n1])
   
    n = 0
    h = hms
    do while(h<=hme)
        rem = numj - mod(j,numj)
        howmany1 = rem + 1
        n = n + 1
        j = js + (h-1)/numj
        do i = 1, size(ocP,1)

            oc(i,n)%ilen = ocP(i,j)%ilen
            oc(i,n)%is = ocP(i,j)%is
            oc(i,n)%ie = ocP(i,j)%ie
            oc(i,n)%rlat = ocP(i,j)%rlat
            oc(i,n)%howmany = howmany1

            nn(1) = ocP(i,j)%ilen
            idist = NX; odist= NX/2+1
            istride = 1; ostride = 1
            inembed = [ocP(i,j)%ilen]; 
            onembed = [ocP(i,j)%ilen/2+1]
            flags = plan_flags

            oc(i,n)%iptr = C_LOC(g2fplans(n)%srin(1,h))
            call c_f_pointer(iptr,g2fplans(n)%oc(i,n)%srin,[TWO*(NX/2+1),

            g2fplans(n)%splan = fftw_plan_many_dft_r2c(ONE, nn, howmany1, &
                                    g2fplans(n)%oc(i,n)%srin, inembed, istride, idist, &
                                    g2fplans(n)%oc(i,n)%scout, onembed, ostride, odist, flags)

            if (.not.c_associated(g2fplans(n)%splan)) &
                call mpp_error('plan_grid_to_fourier: dft_r2c:',trim(null_plan_msg),FATAL)

        end do
        h = h + rem + 1
    end do

    return 
end subroutine set_ocpack_for
        

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
subroutine grid_to_fourier(rinp, coutp, id_in)
!--------------------------------------------------------------------------------   
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
        id = plan_grid_to_fourier(howmany)
    endif
                 
    call mpp_clock_begin(clck_grid_to_fourier)

    if (present(id_in)) id_in = id
    
    howmany = g2fplans(id)%howmany

    g2fplans(id)%rin(1:howmany,:) = rinp(:,:)*RSCALE

    !Transpose
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fplans(id)%plan, g2fplans(id)%rin, g2fplans(id)%trin)
    call mpp_clock_end(clck_g2f_tran)
    
    !Serial FFT
    call mpp_clock_begin(clck_g2f_dft)
    call fftw_execute_dft_r2c(g2fplans(id)%splan, g2fplans(id)%srin, g2fplans(id)%scout) 
    call mpp_clock_end(clck_g2f_dft)

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
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fplans(id)%tplan, g2fplans(id)%srout, g2fplans(id)%tsrout)
    call mpp_clock_end(clck_g2f_tran)

    coutp = g2fplans(id)%cout(1:howmany,1:FLOCAL)

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
    
    if(nplanf2g>max_plans) call mpp_error(grid_fourier_mod,'No: f2g plans > Max_plans', FATAL)

    n = nplanf2g
    plan_fourier_to_grid = n

    f2gplans(n)%howmany = howmany
   
    !Transpose
    n0=[howmany,NLON]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, 1, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (NLON_LOCAL/=local_n1) & 
        call mpp_error('plan_fourier_to_grid', 'NLON_LOCAL/=local_n1', FATAL)

    f2gplans(n)%tcdat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%trin, [NLON, local_n0])
    call c_f_pointer(f2gplans(n)%tcdat, f2gplans(n)%rin, [howmany, local_n1])

    flags = plan_flags

    f2gplans(n)%plan = fftw_mpi_plan_many_transpose(howmany, NLON, 1, &
                            block0, block0, f2gplans(n)%trin, f2gplans(n)%rin, &
                            COMM_FFT, flags)
    if (.not.c_associated(f2gplans(n)%plan)) &
        call mpp_error('plan_fourier_to_grid: transpose1:',trim(null_plan_msg),FATAL)

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
    if (.not.c_associated(f2gplans(n)%splan)) &
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

    f2gplans(n)%cdat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%scouttr, [FTRUNC,local_n1])
    call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%srout, [TWO,FTRUNC,local_n1])
    call c_f_pointer(f2gplans(n)%cdat, f2gplans(n)%tsrout, [TWO,howmany,local_n0])

    f2gplans(n)%tplan = fftw_mpi_plan_many_transpose(FTRUNC, howmany, 2, &
                            block0, block0, f2gplans(n)%tsrout, f2gplans(n)%srout, &
                            COMM_FFT, flags) 
    if (.not.c_associated(f2gplans(n)%tplan)) &
        call mpp_error('plan_fourier_to_grid: transpose2:',trim(null_plan_msg),FATAL)

    f2gplans(n)%r2c = c_loc(f2gplans(n)%tsrout)
    call c_f_pointer(f2gplans(n)%r2c, f2gplans(n)%cout, [howmany, local_n0])

    call save_wisdom()

    call mpp_clock_end(clck_plan_f2g)

    return

end function plan_fourier_to_grid

!--------------------------------------------------------------------------------   
subroutine fourier_to_grid(coutp, rinp, id_in)
!--------------------------------------------------------------------------------   
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
        id = plan_fourier_to_grid(howmany)
    endif
                 
    call mpp_clock_begin(clck_fourier_to_grid)

    if (present(id_in)) id_in = id
    
    howmany = f2gplans(id)%howmany
    
    f2gplans(id)%cout = coutp

    !Transpose Back
    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gplans(id)%tplan, f2gplans(id)%tsrout, f2gplans(id)%srout)
    call mpp_clock_end(clck_f2g_tran)

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

    call mpp_clock_begin(clck_f2g_dft)
    call fftw_execute_dft_c2r(f2gplans(id)%splan, f2gplans(id)%scout, f2gplans(id)%srin) 
    call mpp_clock_end(clck_f2g_dft)

    !Transpose
    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gplans(id)%plan, f2gplans(id)%trin, f2gplans(id)%rin)
    call mpp_clock_end(clck_f2g_tran)
    
    rinp = f2gplans(id)%rin(1:howmany,:)
            
    call mpp_clock_end(clck_fourier_to_grid)
end subroutine fourier_to_grid

!--------------------------------------------------------------------------------   
subroutine end_grid_fourier()
!--------------------------------------------------------------------------------   
    implicit none

    call fftw_mpi_cleanup()

end subroutine end_grid_fourier


!--------------------------------------------------------------------------------   
subroutine fft_1dr2c_serial(rinp, coutp)
!--------------------------------------------------------------------------------   
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

    if (.not.c_associated(plan1d)) &
        call mpp_error('fft_1dr2c_serial:','NULL PLAN',FATAL)

    rinpl = rinp
    call fftw_execute_dft_r2c(plan1d, rinpl, coutpl)
    coutp = coutpl
    
    call fftw_destroy_plan(plan1d)
    call fftw_free(data1d)

    return
end subroutine fft_1dr2c_serial


!--------------------------------------------------------------------------------   
subroutine fft_1dc2c_serial(coutp)
!--------------------------------------------------------------------------------   

    complex :: coutp(:)
    complex(C_DOUBLE_COMPLEX), pointer :: coutpl(:)
    type(C_PTR) :: plan1d, data1d
    integer :: L
    L = size(coutp)

    data1d = fftw_alloc_complex(int(L, C_SIZE_T))
    call c_f_pointer(data1d, coutpl, [L])

    plan1d = fftw_plan_dft_1d(L, coutpl, coutpl, FFTW_FORWARD, FFTW_ESTIMATE) 
    if (.not.c_associated(plan1d)) &
        call mpp_error('fft_1dc2c_serial:','NULL PLAN',FATAL)

    coutpl = coutp
    call fftw_execute_dft(plan1d, coutpl, coutpl)
    coutp = coutpl
    
    call fftw_destroy_plan(plan1d)
    call fftw_free(data1d)

end subroutine fft_1dc2c_serial

end module grid_fourier_mod 

