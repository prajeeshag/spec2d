module ocpack_mod
!--------------------------------------------------------------------------------   
! Module for grid specifications of octahedral packing
! ocpkP is the grid domain in (packed or unpacked)
! ocpkF is the fourier domain (always unpacked)
!

use mpp_mod, only : mpp_error, fatal, warning, note, mpp_init, mpp_npes, mpp_pe, &
        mpp_declare_pelist, mpp_get_current_pelist, mpp_set_current_pelist, mpp_root_pe, &
        mpp_gather, mpp_broadcast, mpp_sync
use strman_mod, only : int2str

implicit none
private

public :: ocpack_typeP, ocpack_typeF, init_ocpack, get_ocpackP, get_ocpackF, oc_npack, &
          oc_isreduced, oc_ny, oc_nx, oc_maxlon, oc_nlat, oc_nfour, get_hem, hem_type, x_block

type ocpack_typeP
    integer :: is, ie, ilen ! longitude start, end and length
    integer :: fs, fe, flen ! fourier, start, end and length
    integer :: f !-> connection to F-grid (ocpkF)
    integer :: g !-> index g is latitude index in a regular globe
end type 

type ocpack_typeF
    integer :: flen ! 
    integer :: i, p !-> connection to P-grid (ocpkP)
    integer :: g !-> index g is latitude index in a regular globe
    integer :: hs
end type 

type hem_type
    integer :: g
    integer :: s, n
end type

type npeblck
    integer :: npes
    integer :: blck
    integer :: rem
    integer, allocatable :: extent(:)
end type npeblck

type(npeblck), allocatable, public :: npesx(:), npesy(:)
integer :: nxpe=0, nype=0, xblock=-1

type(hem_type), allocatable :: jhem(:)

type(ocpack_typeF), dimension(:), allocatable :: ocpkF
type(ocpack_typeP), dimension(:,:), allocatable :: ocpkP

integer :: nplon = 20

integer :: ocny, ocnx, nlat, maxlon, nfour
integer :: num_pack = 2
logical :: reduced=.true.
character(len=1024) :: msg

logical :: initialized=.false., debug=.false.

contains

subroutine get_hem(hem) 
    type(hem_type), intent(out) :: hem(:)

    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')

    if(size(hem,1)==size(jhem,1)) then
        hem = jhem
    else
        call mpp_error(fatal, 'get_jhem: invalid size for argument "hem"!!!!') 
    end if

    return
end subroutine get_hem

subroutine get_ocpackF(ocpack) 
    type(ocpack_typeF), intent(out) :: ocpack(:)

    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')

    if(size(ocpack,1)==size(ocpkF,1)) then
        ocpack = ocpkF(:)
    else
        call mpp_error(fatal, 'get_ocpack: invalid size for argument "ocpack"!!!!') 
    end if

    return
end subroutine get_ocpackF

subroutine get_ocpackP(ocpack) 
    type(ocpack_typeP), intent(out) :: ocpack(:,:)

    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')

    if(size(ocpack,1)==size(ocpkP,1).and.size(ocpack,2)==size(ocpkP,2)) then
        ocpack = ocpkP
    else
        call mpp_error(fatal, 'get_ocpack: invalid size for argument "ocpack"!!!!') 
    end if

    return
end subroutine get_ocpackP

integer function x_block()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    x_block = xblock
    return
end function x_block

integer function oc_nx()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_nx = ocnx
    return
end function oc_nx

integer function oc_nfour()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_nfour = nfour
    return
end function oc_nfour

integer function oc_ny()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_ny = ocny
    return
end function oc_ny

integer function oc_nlat()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_nlat = nlat
    return
end function oc_nlat

integer function oc_maxlon()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_maxlon = maxlon
    return
end function oc_maxlon

integer function oc_npack()
    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')
    oc_npack = num_pack
    return
end function oc_npack

logical function oc_isreduced()
    oc_isreduced = reduced
    return
end function oc_isreduced

!--------------------------------------------------------------------------------   
subroutine init_ocpack(nlat_in, trunc, layout, yextent, xextent, isreduced, ispacked)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: nlat_in
    integer, intent(in) :: trunc
    integer, intent(in) :: layout(2)
    integer, intent(out), optional :: yextent(:), xextent(:)
    logical, intent(in), optional :: isreduced
    logical, intent(in), optional :: ispacked

    integer :: i, j, nplon_new
    integer, dimension(nlat_in) :: lonsperlat
    integer, dimension(:), allocatable :: ny_lcl, nlat_lcl
    integer :: maxpes, ny1, sumny, js, je, n, jj, jlen2, hs

    if(initialized) return

    call mpp_init()

    reduced = .true.
    num_pack = 2
    nlat = nlat_in
    nfour = trunc + 1

    if(present(ispacked).and..not.ispacked) then
        num_pack = 1
    endif

    if(present(isreduced)) reduced = isreduced

    if (mod(nlat,2).ne.0) then
        call mpp_error(fatal, 'ocpack_mod: nlat should be a even number')
    end if

    ocny = nlat/num_pack

    lonsperlat(1) = nplon
    lonsperlat(nlat) = nplon
    do i = 2, nlat/2
        lonsperlat(i) = lonsperlat(i-1) + 4
        lonsperlat(nlat-i+1) = lonsperlat(i)
    end do

    maxlon = maxval(lonsperlat)

    if (.not.reduced) then
        lonsperlat = maxlon
    end if
 
    ocnx = (num_pack-1)*minval(lonsperlat)+maxlon 
    
    call set_valid_layout()

    if (allocated(ny_lcl)) deallocate(ny_lcl)
    if (allocated(nlat_lcl)) deallocate(nlat_lcl)

    do i = 1, nype
        if (npesy(i)%npes==layout(1)) then
            allocate(ny_lcl(layout(1)))
            allocate(nlat_lcl(layout(1)))
            ny_lcl = npesy(i)%extent
            nlat_lcl = ny_lcl*2
            exit
        endif
    end do

    if (.not.allocated(ny_lcl)) then
        write(msg,*) layout(1)
        msg = trim(adjustl(msg))
        call mpp_error(NOTE,'ERROR: unsupported layout:: npes_y = '//trim(msg)//' is not supported')
        if (mpp_root_pe()==mpp_pe()) then
            write(*,*) 'Supported npes_y for this resolution are:'
            write(*,*) npesy(1:nype)%npes
        endif
        call mpp_error(FATAL,'ERROR: unsupported layout:: npes_y = '//trim(msg)//' is not supported')
    endif
        
    if (present(yextent)) then
        if (size(yextent)/=layout(1)) call mpp_error(FATAL,'init_ocpack:size(yextent)/=layout(1)')
        yextent = ny_lcl 
    endif

    i = 1
    do while(i <= nxpe)
        if (npesx(i)%npes==layout(2)) then
            exit
        endif
        i = i + 1
    end do

    if (i>nxpe) then
       write(msg,*) layout(2)
       msg = trim(adjustl(msg))
       call mpp_error(NOTE,'ERROR: unsupported layout:: npes_x = '//trim(msg)//' is not supported')
       if (mpp_root_pe()==mpp_pe()) then
           write(*,*) 'Supported npes_x for this resolution are:'
           write(*,*) npesx(1:nxpe)%npes
       endif
       call mpp_error(FATAL,'ERROR: unsupported layout:: npes_x = '//trim(msg)//' is not supported')
    endif
    
    if (present(xextent)) then
        if (size(xextent)/=layout(2)) call mpp_error(FATAL,'init_ocpack:size(xextent)/=layout(2)')
        xextent=npesx(i)%extent
    endif

    xblock = npesx(i)%blck

    allocate(ocpkF(nlat))
    allocate(ocpkP(num_pack,ocny))

    ocpkP(1,:)%is = 1
    ocpkP(:,:)%g = -1
    ocpkP(:,:)%f = -1

    ocpkP(1,:)%fs = 1
    ocpkP(1,:)%fe = ocpkP(1,:)%fs + maxlon/2 + 1 - 1

    if (num_pack==2) then
        ocpkP(2,:)%ie = ocnx
        ocpkP(2,:)%fs = ocpkP(1,:)%fe + 1
        ocpkP(2,:)%fe = ocpkP(2,:)%fs + maxlon/2 + 1 - 1
        do i = 1, nlat/4
            ocpkP(1,2*i-1)%g = i
            ocpkP(2,2*i-1)%g = nlat/2-i+1
            ocpkP(1,  2*i)%g = nlat-i+1
            ocpkP(2,  2*i)%g = nlat/2+i
        end do
        if (mod(nlat,4)/=0) then
            ocpkP(1,nlat/2)%g = nlat/4 + 1
            ocpkP(2,nlat/2)%g = nlat-nlat/4
        end if
    else
        do i = 1, nlat/4
            ocpkP(1,4*i-3)%g = i
            ocpkP(1,4*i-2)%g = nlat/2-i+1
            ocpkP(1,4*i-1)%g = nlat-i+1
            ocpkP(1,  4*i)%g = nlat/2+i
        end do
        if (mod(nlat,4)/=0) then
            ocpkP(1,nlat-1)%g = nlat/4
            ocpkP(1,nlat  )%g = nlat-nlat/4-1
        end if
    end if

    do i = 1, ocny
        do j = 1, num_pack
            if (j==2) ocpkP(2,i)%is = ocpkP(1,i)%ie + 1
            ocpkP(j,i)%ie = ocpkP(j,i)%is + lonsperlat(ocpkP(j,i)%g) - 1
            ocpkP(j,i)%ilen = lonsperlat(ocpkP(j,i)%g)
            ocpkP(j,i)%flen = ocpkP(j,i)%ilen/2 + 1
        end do
    end do

    ocpkF(:)%flen = nfour 
    ocpkF(:)%g = -1
    ocpkF(:)%i = 0
    ocpkF(:)%p = 0

    je = 0
    jj = 0
    do n = 1, layout(1)
        js = je + 1
        je = js + ny_lcl(n) - 1
        do i = 1, num_pack
            do j = js, je
                jj = jj + 1
                ocpkF(jj)%g = ocpkP(i,j)%g
                ocpkP(i,j)%f = jj
                ocpkF(jj)%i = i
                ocpkF(jj)%p = j
                hs = ocpkP(i,j)%g
                if (hs>nlat/2) hs = nlat - hs + 1
                ocpkF(jj)%hs = hs
            end do
        end do
    end do

    allocate(jhem(nlat/2))

    je = 0
    jj = 0
    do n = 1, layout(1)
        js = je + 1
        je = js + nlat_lcl(n) -1
        jlen2 = nlat_lcl(n)/2
        j = js
        do while (j<je)
            jj = jj + 1
            jhem(jj)%g = ocpkF(j)%g
            jhem(jj)%s = j
            jhem(jj)%n = j + 1
            if (mod(nlat,4)/=0.and.je==nlat) then
                if (j == js + jlen2 -1) then
                    j = j - 1
                    jhem(jj)%n = j + jlen2 + 1
                end if
            end if
            j = j + 2
        end do
    end do
   
    call mpp_error(NOTE,'ocpack initialized-----') 
    initialized = .true.

    return
end subroutine init_ocpack

!--------------------------------------------------------------------------------   
subroutine set_valid_layout()
!--------------------------------------------------------------------------------   
    integer :: npe, i, tmp, j, ocny1
    
    allocate(npesx(ocnx))
    nxpe = 0
    do i = 2, ocnx
        tmp = 0
        npe = 0
        do while (tmp<ocnx)
            npe = npe+1
            tmp = tmp+i
        end do
        nxpe = nxpe + 1
        npesx(nxpe)%npes = npe
        npesx(nxpe)%blck = i
        if ((tmp-ocnx)>0) then
            npesx(nxpe)%rem  = i-(tmp-ocnx)
        else
            npesx(nxpe)%rem = 0
        endif
    end do
    
    call unique(npesx,nxpe)
    
    do i = 1, nxpe
        allocate(npesx(i)%extent(npesx(i)%npes))
        npesx(i)%extent = npesx(i)%blck
        if (npesx(i)%rem>0) npesx(i)%extent(npesx(i)%npes) = npesx(i)%rem
        !if (mpp_pe()==mpp_root_pe()) then
        !    print *, 'npesx=', npesx(i)%npes, 'sum(extent)=', sum(npesx(i)%extent), 'extent=', npesx(i)%extent
        !endif
    end do

    ocny1=ocny
    if (mod(ocny,2)/=0) ocny1=ocny1-1 ! if odd
    ocny1=ocny1/2
    allocate(npesy(ocny1))
    
    do npe = 1, ocny1
        allocate(npesy(npe)%extent(npe))
        npesy(npe)%npes = npe
        npesy(npe)%extent = ocny1/npe
        tmp = mod(ocny1,npe)
        do i = 1, tmp
            npesy(npe)%extent(i)=npesy(npe)%extent(i)+1
        end do
        npesy(npe)%extent = npesy(npe)%extent*2
        if (mod(ocny,2)/=0) npesy(npe)%extent(npe)=npesy(npe)%extent(npe)+1
        !if (mpp_pe()==mpp_root_pe()) then
        !    print *, 'npesy=', npe, 'extent=', npesy(npe)%extent 
        !endif
    end do
    nype = ocny1

    return
end subroutine set_valid_layout

subroutine unique(val,n)
    type(npeblck) :: val(:)
    integer :: n, i, j, npes1, nn
    type(npeblck) :: uniq(size(val))

    npes1 = -1

    nn = n

    n = 0
    do i = 1, nn
        if (npes1 /= val(i)%npes) then
            n = n + 1
            uniq(n)%npes=val(i)%npes
            npes1=val(i)%npes
        endif
    enddo

    uniq(:)%blck = ocnx+1
    uniq(:)%rem = -1
    do j = 1, n
        do i = 1, nn
            if (uniq(j)%npes==val(i)%npes) then
                if (uniq(j)%blck>val(i)%blck) then
                    uniq(j)%blck=val(i)%blck
                    uniq(j)%rem=val(i)%rem
                endif
            endif
        enddo
    enddo

    val(1:n) = uniq(1:n)
    return
end subroutine unique

end  module ocpack_mod

#ifdef test_ocpack
program test
use ocpack_mod
implicit none
integer, parameter :: nlat=94, nfour=62
type(ocpack_typeP), allocatable :: ocpkP(:,:)
type(ocpack_typeF), allocatable :: ocpkF(:)
type(hem_type), allocatable :: jh(:)
integer :: i, j, npes

!call init_ocpack(nlat,nfour,maxlon=192,ispacked=.false.,isreduced=.true.)
!call init_ocpack(nlat,nfour,ispacked=.false.,isreduced=.true.)
print *, 'enter npes in y:'
read(*,*) npes
call init_ocpack(nlat,nfour,npes,ispacked=.true.,isreduced=.true.)

allocate(ocpkF(oc_nlat()))
allocate(ocpkP(oc_npack(),oc_ny()))
allocate(jh(oc_nlat()/2))

call get_ocpackF(ocpkF)
call get_ocpackP(ocpkP)
call get_hem(jh)

do i = 1, oc_ny()
    do j = 1, oc_npack()
        print '(8(i4,2x))', i, j, ocpkP(j,i)%is, ocpkP(j,i)%ie, &
                      ocpkP(j,i)%fs, ocpkP(j,i)%fe, ocpkP(j,i)%g, ocpkP(j,i)%f
    end do
end do

    print *, 'ocpkF+++++++++'
do i = 1, oc_nlat()
    print '(5(i4,2x))', ocpkF(i)%p, ocpkF(i)%i, ocpkF(i)%hs, ocpkF(i)%g
end do
    print *, 'Hemi+++++++++'
do i = 1, oc_nlat()/2
    print '(5(i4,2x))', jh(i)%g, jh(i)%s, jh(i)%n
end do

end program test
#endif

