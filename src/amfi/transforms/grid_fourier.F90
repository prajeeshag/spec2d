module grid_fourier_mod

!--------------------------------------------------------------------------------   
! Module for grid to fourier and fourier to grid transformation. 
!
! sF --c2r--> rsF --trnsps--> rsFT --c2r--> sFT --shffl--> FT --fft--> GT --trnsps--> G
!
use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_pe, mpp_npes, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
        mpp_sync, mpp_root_pe, CLOCK_MODULE, CLOCK_ROUTINE, mpp_gather, mpp_broadcast, &
        mpp_set_current_pelist, mpp_get_current_pelist, mpp_declare_pelist

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_layout, &
        mpp_get_domain_extents

use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, WARNING, NOTE, &
                    open_file, close_file, file_exist

use ocpack_mod, only : ocpack_typeP, ocpack_typeF, oc_nx, oc_ny, oc_nlat, oc_maxlon, &
                       npack=>oc_npack, oc_isreduced, oc_nfour, get_ocpackP

use spec_comm_mod, only : split_pelist, spec_comm_allgather, spec_comm_pe

use strman_mod, only : int2str

implicit none
include 'fftw3-mpi.f03'

private

integer :: nplang2f=0, nplanf2g=0
integer, parameter :: max_plans=10
integer, parameter :: rank=2
integer(C_INTPTR_T), parameter :: TWO=2, ONE=1

type ocplan_type
    type(C_PTR) :: plan
    type(C_PTR) :: iptr, optr
    real(C_DOUBLE), pointer :: oG(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: oF(:,:)
    integer :: sf0, ef0, sh0, eh0
end type ocplan_type

type plan_type
    integer(C_INTPTR_T) :: howmany
    type(C_PTR) :: tplan1, tplan2, fplan
    type(C_PTR) :: t1dat, t2dat, r2c
    real(C_DOUBLE), pointer :: G(:,:), GT(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: FT(:,:), sFT(:,:)
    real(C_DOUBLE), pointer :: rsFT(:,:,:), rsF(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: sF(:,:)
    type(ocplan_type), allocatable :: fp(:,:)
endtype plan_type

integer(C_INTPTR_T) :: MXLON, MXLON2, NFOUR, NFOUR2, MXFOUR, MXFOUR2
integer(C_INTPTR_T) :: NX_LOCAL, FLOCAL, FLOCAL2, NX, FBLOCK, FBLOCK2
integer(C_INTPTR_T) :: block0=FFTW_MPI_DEFAULT_BLOCK
integer :: COMM_FFT, COMM_PARENT
real, allocatable :: RSCALE(:,:)

integer, allocatable :: Tshuffle(:), Tshuffle2(:)
logical :: shuffle=.false.

logical :: debug=.false.
integer :: clck_grid_to_fourier, clck_fourier_to_grid, clck_g2f_tran, clck_g2f_dft, &
           clck_f2g_tran, clck_f2g_dft, clck_plan_g2f, clck_plan_f2g

integer :: plan_level = 3, plan_flags, id_grid2four, id_four2grid
logical :: initialized=.false.

character (len=256) :: wsdmfnm='fftw.wisdom'
character (len=256) :: null_plan_msg

type(plan_type) :: g2fp(max_plans), f2gp(max_plans)

integer :: jsp, jep, jlenp 
integer :: isp, iep, ilenp
integer :: layout(2)
integer, allocatable :: pes(:)

type(ocpack_typeP), allocatable :: ocP(:,:)

public :: init_grid_fourier, end_grid_fourier, grid_to_fourier, fourier_to_grid, &
          plan_grid_to_fourier, plan_fourier_to_grid

namelist/grid_fourier_nml/plan_level, debug

contains

!--------------------------------------------------------------------------------   
subroutine init_grid_fourier (domain, trunc, isf, flen, fextent, Tshuff)
!--------------------------------------------------------------------------------   
    implicit none
    type(domain2D), intent(in) :: domain
    integer, intent(in) :: trunc ! fourier truncation
    integer, intent(out) :: flen ! No: fouriers in this proc
    integer, intent(out) :: isf ! Starting of the fourier in this proc
    integer, intent(out) :: fextent(:)
    integer, intent(out), optional :: Tshuff(trunc+1) !if present shuffle the fourier for 
                                                         !load balance for triangular truncation
                                                         !and give the order of shuffled fouriers
                                                         !in Tshuff
    integer :: unit, i, k, stat, n, flags, t, j
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2)
    integer(C_INTPTR_T) :: howmany
    real :: rnlon
    integer :: jee, pe, ierr, npes
    integer, allocatable :: allpes(:), extent(:), xextent(:)

    unit = open_namelist_file()

    read(unit, nml=grid_fourier_nml,iostat=stat)

    call close_file(unit)

    call mpp_get_layout(domain,layout)

    allocate(extent(layout(1)))
    allocate(xextent(layout(2)))

    if (size(fextent)/=layout(2)) call mpp_error(FATAL,'init_grid_fourier: fextent '//&
                                    'array size must be equal to layout(2)')

    call mpp_get_compute_domain(domain,jsp,jep,isp,iep)

    call mpp_get_domain_extents(domain,extent,xextent)

    jlenp = jep - jsp + 1
    ilenp = iep - isp + 1

    allocate(allpes(mpp_npes()))
    allocate(pes(layout(2)))

    call mpp_get_current_pelist(commid=COMM_PARENT)

    COMM_FFT = 0    
    jee = 0
    do j = 1, size(extent)
        jee = jee + extent(j)
        call split_pelist(jee==jep, pes, npes, COMM_FFT)
    end do

    allocate(ocP(npack(),oc_ny()))
    call get_ocpackP(ocP)

    allocate(RSCALE(jsp:jep,isp:iep))

    do i = isp, iep
        do j = jsp, jep
            rnlon = 1./ocP(1,j)%ilen
            if (i>ocP(1,j)%ilen) rnlon = 1./ocP(2,j)%ilen
            RSCALE(j,i) = rnlon
        end do
    end do

    MXLON    = oc_maxlon()
    NX       = oc_nx()
    NFOUR    = trunc + 1
    MXFOUR   = MXLON/2+1
    NX_LOCAL = ilenp
    FBLOCK = ceiling(real(NFOUR)/layout(2))

    NFOUR2  = NFOUR  * npack()
    MXLON2  = MXLON  * npack()
    MXFOUR2 = MXFOUR * npack() 
    FBLOCK2 = FBLOCK * npack()    

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

    if(mod(NX,layout(2))/=0) &
        call mpp_error('grid_fourier_mod','npes in x-dir should be a factor of NX='//trim(int2str(NX)), FATAL)

    call fftw_mpi_init()

    clck_grid_to_fourier = mpp_clock_id('grid_to_fourier')
    clck_fourier_to_grid = mpp_clock_id('fourier_to_grid')
    clck_plan_g2f = mpp_clock_id('plan_grid_to_fourier')
    clck_plan_f2g = mpp_clock_id('plan_fourier_to_grid')
    clck_g2f_tran = mpp_clock_id('g2f_tran')
    clck_g2f_dft = mpp_clock_id('g2f_dft')
    clck_f2g_tran = mpp_clock_id('f2g_tran')
    clck_f2g_dft = mpp_clock_id('f2g_dft')

    if (import_wisdom(int(NX))==1) then
        plan_flags = ior(plan_flags,FFTW_WISDOM_ONLY)
    endif

    null_plan_msg = 'NULL PLAN: try rerunning after removing the wisdom file '//trim(wsdmfnm)

    allocate(Tshuffle(NFOUR))
    allocate(Tshuffle2(NFOUR2))
    Tshuffle = [(i, i=1,NFOUR)]

    if (present(Tshuff)) then
        if (mod(NFOUR,2)==0) then
            k = NFOUR + 1
            do i = 1, NFOUR/2
                k = k - 1 
                Tshuffle(2*(i-1)+1) = i
                Tshuffle(2*i) = k
            enddo
        else
            Tshuffle(NFOUR) = 1
            k = NFOUR + 1
            do i = 1, (NFOUR-1)/2
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

    Tshuffle2(1:NFOUR2:npack()) = Tshuffle
    if(npack()==2) Tshuffle2(2:NFOUR2:npack()) = Tshuffle+MXFOUR 

    howmany = 1
    n0=[NX,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, ONE, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    local_n1_prev = local_n1

    if (ilenp/=local_n0) & 
        call mpp_error('init_grid_fourier', 'ilenp/=local_n0 : '//trim(int2str(ilenp))//&
                        ', '//trim(int2str(local_n0)), FATAL)

    n0=[howmany,NFOUR]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, TWO, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n1_prev/=local_n0) &
        call mpp_error('init_grid_fourier', 'local_n1_prev/=local_n0', FATAL)

    isf     = local_1_start
    flen    = local_n1
    FLOCAL  = local_n1 
    FLOCAL2 = FLOCAL * npack()

    call mpp_set_current_pelist(pes,no_sync=.true.)
    call spec_comm_allgather([flen], fextent)

    call mpp_set_current_pelist()

    initialized = .true.
    call mpp_error('init_grid_fourier', '----Initialized----', NOTE)

    return
end subroutine init_grid_fourier


!--------------------------------------------------------------------------------   
function plan_grid_to_fourier(howmany,grid,four)
!--------------------------------------------------------------------------------   
    implicit none
    integer(C_INTPTR_T), intent(in) :: howmany
    real, pointer, optional :: grid(:,:)
    complex, pointer, optional :: four(:,:)
    integer :: plan_grid_to_fourier, n, flags, t, clck_transpose
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1, local_n1_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2), howmany2, alloc_local1 
    integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)
    

    if (.not.initialized) call mpp_error(FATAL,'grid_fourier_mod: not initialized')

    howmany2 = howmany * npack()

    call mpp_clock_begin(clck_plan_g2f)
    if (howmany<1) call mpp_error('plan_grid_to_fourier', 'howmany cannot be Zero', FATAL)

    nplang2f = nplang2f + 1
    
    if(nplang2f>max_plans) call mpp_error('grid_fourier_mod','No: g2f plans > Max_plans', FATAL)

    n = nplang2f
    plan_grid_to_fourier = n

    g2fp(n)%howmany = howmany
 
    n0=[howmany,NFOUR2]
    alloc_local1 = fftw_mpi_local_size_many_transposed(rank, n0, two, &
                   block0, FBLOCK2, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)
    !Transpose
    n0=[NX,howmany]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, one, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    alloc_local = max(alloc_local, alloc_local1)

    if (NX_LOCAL/=local_n0) & 
        call mpp_error('plan_grid_to_fourier', 'NX_LOCAL/=local_n0', FATAL)

    flags = plan_flags

    g2fp(n)%t1dat = fftw_alloc_complex(alloc_local)
    g2fp(n)%t2dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%G, [howmany, local_n0])
    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%GT, [NX, local_n1])

    if (present(grid)) call c_f_pointer(g2fp(n)%t1dat, grid, [howmany, local_n0])

    g2fp(n)%tplan1 = fftw_mpi_plan_many_transpose(NX, howmany, one, &
                            block0, block0, g2fp(n)%G, g2fp(n)%GT, &
                            COMM_FFT, flags)
    
    if (.not.c_associated(g2fp(n)%tplan1)) &
        call mpp_error('plan_grid_to_fourier: transpose1:',trim(null_plan_msg),FATAL)

    !multi-threaded shared memory fft
    call c_f_pointer(g2fp(n)%t2dat, g2fp(n)%FT, [MXFOUR2,local_n1])
 
    call set_ocplan_g2f(local_n1, local_1_start, howmany, n)

    local_n1_prev = local_n1
    !Transpose back

    n0=[howmany,NFOUR2]

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, two, &
                   block0, FBLOCK2, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n1_prev/=local_n0) &
        call mpp_error('plan_grid_to_fourier', 'local_n1_prev/=local_n0', FATAL)

    if (FLOCAL2/=local_n1) call mpp_error('plan_grid_to_fourier', 'FLOCAL2/=local_n1', FATAL)

    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%sFT, [NFOUR2,local_n0])
    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%rsFT, [TWO,NFOUR2,local_n0])
    call c_f_pointer(g2fp(n)%t1dat, g2fp(n)%rsF, [TWO,howmany,local_n1])

    g2fp(n)%tplan2 = fftw_mpi_plan_many_transpose(howmany, NFOUR2, two, &
                            block0, FBLOCK2, g2fp(n)%rsFT, g2fp(n)%rsF, &
                            COMM_FFT, flags) 

    if (.not.c_associated(g2fp(n)%tplan2)) &
        call mpp_error('plan_grid_to_fourier: transpose2:',trim(null_plan_msg),FATAL)

    g2fp(n)%r2c = c_loc(g2fp(n)%rsF)
    call c_f_pointer(g2fp(n)%r2c, g2fp(n)%sF, [howmany2, FLOCAL])

    if (present(four)) call c_f_pointer(g2fp(n)%r2c, four, [howmany2, FLOCAL])

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

    write(tmpc,*) layout(1)
    wsdmfnm = trim(wsdmfnm)//'_'//trim(adjustl(tmpc)) 

    write(tmpc,*) layout(2)
    wsdmfnm = trim(wsdmfnm)//'x'//trim(adjustl(tmpc)) 

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

    import_wisdom = isuccess
    return
end function import_wisdom 

!--------------------------------------------------------------------------------   
subroutine save_wisdom()
!--------------------------------------------------------------------------------   
    integer :: unitw

    call fftw_mpi_gather_wisdom(COMM_PARENT)

    if (spec_comm_pe(COMM_PARENT)==0) then
        unitw = open_file(file=trim(wsdmfnm),action='write',form='ascii') 
        call export_wisdom_to_file(unitw)
        call close_file(unitw)
    endif

    return
end subroutine save_wisdom 

function num_plans_for(hmlen, shm, howmany)
    integer(kind=C_INTPTR_T), intent(in) :: hmlen, shm, howmany
    integer :: num_plans_for
    integer :: hms, hme, nk, h, rem, remh

    hms = shm + 1 ! shmn is start of homany from fftw, which starts from 0.
    hme = hms + hmlen - 1
    nk = howmany/jlenp
    
    num_plans_for = 0
    h = hms
    remh = hmlen
    do while(h <= hme)
        rem = min((1+(h-1)/nk)*nk - h + 1,remh)
        remh = remh - rem
        h = h + rem
        num_plans_for = num_plans_for+1
    end do

    return
end function num_plans_for

!--------------------------------------------------------------------------------   
subroutine set_ocplan_g2f(hlen, sh, howmany, np)
!--------------------------------------------------------------------------------   
    integer(kind=C_INTPTR_T), intent(in) :: hlen, sh, howmany
    integer, intent(in) :: np
    integer(kind=C_INTPTR_T) :: nk
    integer :: hs, he, j, h, n, rem, i, flen, hh
    integer :: nn(1), idist, odist, istride, ostride, inembed(1), onembed(1)
    integer :: flags, is, numocpln, nlon, is_g, is_f, hmny, remh

    hs = sh + 1 
    he = hs + hlen - 1
    nk = howmany/jlenp 

    numocpln = num_plans_for(hlen, sh, howmany)

    allocate(g2fp(np)%fp(npack(),numocpln))

    n = 0
    h = hs
    remh = hlen
    do while(h<=he)
        j = jsp + (h-1)/nk
        rem = min((1+(h-1)/nk)*nk - h + 1,remh)
        remh = remh - rem
        hmny = rem
        n = n + 1
        hh = h - hs + 1
        do i = 1, npack()
            nlon = ocP(i,j)%ilen
            flen = ocP(i,j)%flen
            
            nn(1) = nlon
            idist = NX; odist = MXFOUR2
            istride = 1; ostride = 1
            inembed = nlon; onembed = flen
            flags = plan_flags

            is_g = ocP(i,j)%is
            g2fp(np)%fp(i,n)%iptr = c_loc(g2fp(np)%GT(is_g,hh))
            call c_f_pointer(g2fp(np)%fp(i,n)%iptr, &
                    g2fp(np)%fp(i,n)%oG, [int(NX),hmny])

            is_f = 1 + (i-1)*(MXFOUR)
            g2fp(np)%fp(i,n)%optr = c_loc(g2fp(np)%FT(is_f,hh))
            call c_f_pointer(g2fp(np)%fp(i,n)%optr, &
                    g2fp(np)%fp(i,n)%oF, [int(MXFOUR2),hmny])

            g2fp(np)%fp(i,n)%sf0 = is_f + flen
            g2fp(np)%fp(i,n)%ef0 = is_f + MXFOUR - 1
            g2fp(np)%fp(i,n)%sh0 = hh 
            g2fp(np)%fp(i,n)%eh0 = hh + hmny - 1

            !if(mpp_pe()==mpp_root_pe()) print '(10(I5,1x))', nk, j, h, hmny, hh, hh + hmny - 1

            g2fp(np)%fp(i,n)%plan = fftw_plan_many_dft_r2c(1, nn, hmny, &
                            g2fp(np)%fp(i,n)%oG, inembed, istride, idist, &
                            g2fp(np)%fp(i,n)%oF, onembed, ostride, odist, flags) 

            if (.not.c_associated(g2fp(np)%fp(i,n)%plan)) &
                call mpp_error('set_ocplan_g2f: '//trim(int2str(i))//' '//trim(int2str(n))//' :', &
                                trim(null_plan_msg), FATAL)

        end do
        h = h + rem
    end do
             
    return
end subroutine set_ocplan_g2f


!--------------------------------------------------------------------------------   
subroutine set_ocplan_f2g(hlen, sh, howmany, np)
!--------------------------------------------------------------------------------   
    integer(kind=C_INTPTR_T), intent(in) :: hlen, sh, howmany
    integer, intent(in) :: np
    integer(kind=C_INTPTR_T) :: nk
    integer :: hs, he, j, h, n, rem, i, flen, hh
    integer :: nn(1), idist, odist, istride, ostride, inembed(1), onembed(1)
    integer :: flags, is, numocpln, nlon, is_g, is_f, hmny, remh

    hs = sh + 1 
    he = hs + hlen - 1
    nk = howmany/jlenp 

    numocpln = num_plans_for(hlen, sh, howmany)

    allocate(f2gp(np)%fp(npack(),numocpln))

    n = 0
    h = hs
    remh = hlen
    do while(h<=he)
        j = jsp + (h-1)/nk
        rem = min((1+(h-1)/nk)*nk - h + 1,remh)
        remh = remh - rem
        hmny = rem
        n = n + 1
        hh = h - hs + 1
        do i = 1, npack()
            nlon = ocP(i,j)%ilen
            flen = ocP(i,j)%flen
            
            nn(1) = nlon
            idist = NX; odist = MXFOUR2
            istride = 1; ostride = 1
            inembed = nlon; onembed = flen
            flags = plan_flags

            is_g = ocP(i,j)%is
            f2gp(np)%fp(i,n)%iptr = c_loc(f2gp(np)%GT(is_g,hh))
            call c_f_pointer(f2gp(np)%fp(i,n)%iptr, &
                    f2gp(np)%fp(i,n)%oG, [int(NX),hmny])

            is_f = 1 + (i-1)*(MXFOUR)
            f2gp(np)%fp(i,n)%optr = c_loc(f2gp(np)%FT(is_f,hh))
            call c_f_pointer(f2gp(np)%fp(i,n)%optr, &
                    f2gp(np)%fp(i,n)%oF, [int(MXFOUR2),hmny])

            f2gp(np)%fp(i,n)%sf0 = is_f + flen
            f2gp(np)%fp(i,n)%ef0 = is_f + MXFOUR - 1
            f2gp(np)%fp(i,n)%sh0 = hh 
            f2gp(np)%fp(i,n)%eh0 = hh + hmny - 1

            f2gp(np)%fp(i,n)%plan = fftw_plan_many_dft_c2r(1, nn, hmny, &
                            f2gp(np)%fp(i,n)%oF, onembed, ostride, odist, & 
                            f2gp(np)%fp(i,n)%oG, inembed, istride, idist, flags)

            if (.not.c_associated(f2gp(np)%fp(i,n)%plan)) &
                call mpp_error('set_ocplan_f2g: '//trim(int2str(i))//' '//trim(int2str(n))//' :', &
                                trim(null_plan_msg), FATAL)

        end do
        h = h + rem
    end do
             
    return
end subroutine set_ocplan_f2g


!--------------------------------------------------------------------------------     
subroutine grid_to_fourier(Gp, sFp, id_in)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(in), optional :: Gp(:,:) ! nk*jlenp, ilenp 
    complex, intent(out), optional :: sFp(:,:) ! fourier, lat, lev
    integer, intent(inout), optional :: id_in
    integer :: id, i, j, ct, k
    integer(C_INTPTR_T) :: howmany, howmany2
    type(C_PTR) :: cptr1, cptr2
    real, pointer :: Gj(:,:,:)

    id = 0

    if (present(id_in)) id = id_in

    if (id<1) then
        if (.not.present(sFp).or..not.present(Gp)) then
            call mpp_error(FATAL,'grid_to_fourier: if sFp or Gp is not present then id_in '//&
                            'must be present and should be positive value')
        endif
        howmany = size(Gp,1)
        howmany2 = howmany*npack()

        do i = 1, nplang2f
            if (howmany==g2fp(i)%howmany) then
                id = i
                exit
            endif
        enddo
    else
        howmany = g2fp(id)%howmany
        howmany2 = howmany*npack()
    endif
    
    if (id<1) then
        id = plan_grid_to_fourier(howmany)
    end if
                 
    call mpp_clock_begin(clck_grid_to_fourier)

    if (present(id_in)) id_in = id
    
    howmany = g2fp(id)%howmany

    if (present(Gp)) g2fp(id)%G=Gp

    cptr1 = C_LOC(g2fp(id)%G)
    call c_f_pointer(cptr1,Gj,[int(howmany/jlenp),jlenp,ilenp])

    do k = 1, howmany/jlenp
        Gj(k,:,:) = Gj(k,:,:)*RSCALE(:,:)
    end do

    !Transpose
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fp(id)%tplan1, g2fp(id)%G, g2fp(id)%GT)
    
    call mpp_clock_end(clck_g2f_tran)

    !Serial FFT
    g2fp(id)%FT(:,:) = cmplx(0.,0.)
    call mpp_clock_begin(clck_g2f_dft)
    call execute_ocfft_r2c(id)
    !call fftw_execute_dft_r2c(g2fp(id)%fplan, g2fp(id)%GT, g2fp(id)%FT) 
    call mpp_clock_end(clck_g2f_dft)

    !Truncation
    g2fp(id)%sFT(:,:) = cmplx(0.,0.)
    do i = 1, NFOUR2
        j = Tshuffle2(i)
        g2fp(id)%sFT(i,:) = g2fp(id)%FT(j,:)!Truncation & Shuffle
    enddo

    !Transpose Back
    call mpp_clock_begin(clck_g2f_tran)
    call fftw_mpi_execute_r2r(g2fp(id)%tplan2, g2fp(id)%rsFT, g2fp(id)%rsF)
    call mpp_clock_end(clck_g2f_tran)

    if (present(sFp)) sFp = g2fp(id)%sF(1:howmany2,1:FLOCAL)

    call mpp_clock_end(clck_grid_to_fourier)
            
end subroutine grid_to_fourier


!--------------------------------------------------------------------------------   
function plan_fourier_to_grid(howmany,four,grid)
!--------------------------------------------------------------------------------   
    implicit none
    integer(C_INTPTR_T), intent(in) :: howmany ! howmany=jlenp*whatever
    real, pointer, optional :: grid(:,:)
    complex, pointer, optional :: four(:,:)
    integer :: plan_fourier_to_grid, n, flags, t, clck_transpose
    integer(C_INTPTR_T) :: local_n0, local_0_start, local_1_start, local_n1
    integer(C_INTPTR_T) :: local_n0_prev
    integer(C_INTPTR_T) :: alloc_local, n0(2), alloc_local1
    integer :: inembed(1), onembed(1), istride, ostride, idist, odist, nn(1)
    integer(C_INTPTR_T) :: howmany2

    if (.not.initialized) call mpp_error(FATAL,'grid_fourier_mod: not initialized')

    call mpp_clock_begin(clck_plan_f2g)

    howmany2 = howmany*npack()

    if (howmany<1) call mpp_error('plan_fourier_to_grid', 'howmany cannot be Zero', FATAL)

    nplanf2g = nplanf2g + 1
    
    if(nplanf2g>max_plans) call mpp_error('grid_fourier_mod','No: f2g plans > Max_plans', FATAL)

    n = nplanf2g
    plan_fourier_to_grid = n

    f2gp(n)%howmany = howmany
   
    !Transpose
    n0=[NFOUR2,howmany]
    alloc_local1 = fftw_mpi_local_size_many_transposed(rank, n0, two, &
                   FBLOCK2, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    n0=[howmany,NX]
    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, one, &
                   block0, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    alloc_local = max(alloc_local,alloc_local1)

    if (NX_LOCAL/=local_n1) & 
        call mpp_error('plan_fourier_to_grid', 'NX_LOCAL/=local_n1', FATAL)

    f2gp(n)%t1dat = fftw_alloc_complex(alloc_local)
    f2gp(n)%t2dat = fftw_alloc_complex(alloc_local)

    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%GT, [NX, local_n0])
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%G, [howmany, local_n1])
    if (present(grid)) call c_f_pointer(g2fp(n)%t1dat, grid, [howmany, local_n1])

    flags = plan_flags

    f2gp(n)%tplan1 = fftw_mpi_plan_many_transpose(howmany, NX, one, &
                            block0, block0, f2gp(n)%GT, f2gp(n)%G, &
                            COMM_FFT, flags)
    if (.not.c_associated(f2gp(n)%tplan1)) &
        call mpp_error('plan_fourier_to_grid: transpose1:',trim(null_plan_msg),FATAL)

    !multi-threaded shared memory fft
    call c_f_pointer(f2gp(n)%t2dat, f2gp(n)%FT, [MXFOUR2,local_n0])

    call set_ocplan_f2g(local_n0, local_0_start, howmany, n)
   
    local_n0_prev = local_n0 

    !Transpose back
    n0=[NFOUR2,howmany]
    alloc_local = fftw_mpi_local_size_many_transposed(rank, n0, two, &
                   FBLOCK2, block0, COMM_FFT, local_n0, &
                   local_0_start, local_n1, local_1_start)

    if (local_n0/=FLOCAL2) call mpp_error(FATAL,'plan_fourier_to_grid: local_n0/=FLOCAL2')

    if (local_n0_prev /= local_n1) &
        call mpp_error('plan_fourier_to_grid', 'local_n0_prev /= local_n1', FATAL)

    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%sFT, [NFOUR2,local_n1])
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%rsFT, [TWO,NFOUR2,local_n1])
    call c_f_pointer(f2gp(n)%t1dat, f2gp(n)%rsF, [TWO,howmany,FLOCAL2])
    
    
    f2gp(n)%tplan2 = fftw_mpi_plan_many_transpose(NFOUR2, howmany, two, &
                            FBLOCK2, block0, f2gp(n)%rsF, f2gp(n)%rsFT, &
                            COMM_FFT, flags) 

    if (.not.c_associated(f2gp(n)%tplan2)) &
        call mpp_error('plan_fourier_to_grid: transpose2:',trim(null_plan_msg),FATAL)

    f2gp(n)%r2c = c_loc(f2gp(n)%rsF)
    call c_f_pointer(f2gp(n)%r2c, f2gp(n)%sF, [howmany2, FLOCAL])
    if (present(four)) call c_f_pointer(g2fp(n)%r2c, four, [howmany2, FLOCAL])

    call mpp_clock_end(clck_plan_f2g)
    
    return

end function plan_fourier_to_grid



!--------------------------------------------------------------------------------   
subroutine fourier_to_grid(sFp, Gp, id_in)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(out), optional :: Gp(:,:) ! howmany, NX
    complex, intent(in), optional :: sFp(:,:) ! howmany2, FLOCAL
    integer, intent(inout), optional :: id_in
    integer :: id, i, j, ii, jj
    integer(C_INTPTR_T) :: howmany2, howmany 

    id = 0

    if (present(id_in)) id = id_in

    if (id<1) then
        if (.not.present(sFp).or..not.present(Gp)) then
            call mpp_error(FATAL,'fourier_to_grid: if sFp or Gp is not present then id_in '//&
                            'must be present and should be positive value')
        endif
        howmany = size(Gp,1)
        howmany2 = howmany*npack()

        do i = 1, nplanf2g
            if (howmany==f2gp(i)%howmany) then
                id = i
                exit
            endif
        enddo
    else
        howmany = f2gp(id)%howmany
        howmany2 = howmany*npack()
    endif
    
    if (id<1) then
        id = plan_fourier_to_grid(howmany)
    endif
                 
    call mpp_clock_begin(clck_fourier_to_grid)

    if (present(id_in)) id_in = id
   
    if (present(sFp)) f2gp(id)%sF = sFp

    !Transpose Back

    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gp(id)%tplan2, f2gp(id)%rsF, f2gp(id)%rsFT)
    call mpp_clock_end(clck_f2g_tran)

    !Serial FFT
    f2gp(id)%FT(:,:) = cmplx(0.,0.)
    do i = 1, NFOUR2
        j = Tshuffle2(i)
        f2gp(id)%FT(j,:) = f2gp(id)%sFT(i,:) !Truncation & Shuffle
    enddo

    call mpp_clock_begin(clck_f2g_dft)
    call execute_ocfft_c2r(id)
    !call fftw_execute_dft_c2r(f2gp(id)%fplan, f2gp(id)%FT, f2gp(id)%GT) 
    call mpp_clock_end(clck_f2g_dft)

    !Transpose
    call mpp_clock_begin(clck_f2g_tran)
    call fftw_mpi_execute_r2r(f2gp(id)%tplan1, f2gp(id)%GT, f2gp(id)%G)
    call mpp_clock_end(clck_f2g_tran)
    
    if (present(Gp)) Gp = f2gp(id)%G(1:howmany,1:NX_LOCAL)
            
    call mpp_clock_end(clck_fourier_to_grid)
end subroutine fourier_to_grid

!--------------------------------------------------------------------------------   
subroutine execute_ocfft_r2c(id)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: id
    integer :: i, n, sf, ef, sh, eh

    do i = 1, size(g2fp(id)%fp,1)
        do n = 1, size(g2fp(id)%fp,2)
            call fftw_execute_dft_r2c(g2fp(id)%fp(i,n)%plan, g2fp(id)%fp(i,n)%oG, g2fp(id)%fp(i,n)%oF)
            sf = g2fp(id)%fp(i,n)%sf0
            ef = g2fp(id)%fp(i,n)%ef0
            sh = g2fp(id)%fp(i,n)%sh0
            eh = g2fp(id)%fp(i,n)%eh0
            g2fp(id)%FT(sf:ef,sh:eh) = cmplx(0.,0.)
        end do
    end do

    return
end subroutine execute_ocfft_r2c

!--------------------------------------------------------------------------------   
subroutine execute_ocfft_c2r(id)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: id
    integer :: i, n

    do i = 1,  size(f2gp(id)%fp,1)
        do n = 1, size(f2gp(id)%fp,2)
            call fftw_execute_dft_c2r(f2gp(id)%fp(i,n)%plan, f2gp(id)%fp(i,n)%oF, f2gp(id)%fp(i,n)%oG)
        end do
    end do

    return
end subroutine execute_ocfft_c2r

!--------------------------------------------------------------------------------   
subroutine end_grid_fourier()
!--------------------------------------------------------------------------------   
    implicit none

    call save_wisdom()
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

