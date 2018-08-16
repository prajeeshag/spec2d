module ocpack_mod
!--------------------------------------------------------------------------------   
! Module for grid specifications of octahedral packing
! ocpkP is the grid domain in (packed or unpacked)
! ocpkF is the fourier domain (always unpacked)
!

use mpp_mod, only : mpp_error, fatal, warning, note, mpp_init
use strman_mod, only : int2str

implicit none
private

public :: ocpack_typeP, ocpack_typeF, init_ocpack, get_ocpackP, get_ocpackF, oc_npack, &
          oc_isreduced, oc_ny, oc_nx, oc_maxlon, oc_nlat, oc_nfour

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
end type 

type(ocpack_typeF), dimension(:), allocatable :: ocpkF
type(ocpack_typeP), dimension(:,:), allocatable :: ocpkP

integer :: nplon = 20

integer :: ocny, ocnx, nlat, maxlon, nfour
integer :: num_pack = 2
logical :: reduced=.true.

logical :: initialized=.false.

contains

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

subroutine init_ocpack(nlat_in, trunc, max_lon, isreduced, ispacked)

    integer, intent(in) :: nlat_in
    integer, intent(in) :: trunc
    integer, intent(in), optional :: max_lon
    logical, intent(in), optional :: isreduced
    logical, intent(in), optional :: ispacked

    integer :: i, j, nplon_new
    integer, dimension(nlat_in) :: lonsperlat

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

    if(num_pack==2) reduced = .true.

    if (mod(nlat,2).ne.0) then
        call mpp_error(fatal, 'ocpack_mod: nlat should be a multiple of 2')
    end if

    if (reduced.and.present(max_lon)) then
        nplon_new = max_lon-(4*(nlat/2-1))
        if (nplon_new<8) call mpp_error(FATAL, 'nplon < 8, increase max_lon')
        if (nplon_new/=nplon) then
            nplon = nplon_new
            call mpp_error(warning, 'changing default number of (20) pole '// &
                      'longitudes (nplon) to: '//int2str(nplon))
        end if
    end if

    lonsperlat(1) = nplon
    lonsperlat(nlat) = nplon
    do i = 2, nlat/2
        lonsperlat(i) = lonsperlat(i-1) + 4
        lonsperlat(nlat-i+1) = lonsperlat(i)
    end do

    maxlon = maxval(lonsperlat)

    if (.not.reduced) then
        lonsperlat = maxlon
        if (present(max_lon)) then
            lonsperlat = max_lon
            maxlon = max_lon
        endif
    end if
 
    ocnx = (num_pack-1)*nplon+maxlon 
    ocny = nlat/num_pack
    
    allocate(ocpkF(nlat))
    allocate(ocpkP(num_pack,ocny))

    ocpkF(:)%flen = nfour 
    ocpkF(:)%g = -1
    ocpkF(:)%i = 0
    ocpkF(:)%p = 0

    do i = 1, nlat/4
        ocpkF(4*i-3)%g = i
        ocpkF(4*i-2)%g = nlat/2-i+1
        ocpkF(4*i-1)%g = nlat-i+1
        ocpkF(  4*i)%g = nlat/2+i

        ocpkF(4*i-3)%i = 1
        ocpkF(4*i-2)%i = 2
        ocpkF(4*i-1)%i = 1
        ocpkF(  4*i)%i = 2

        ocpkF(4*i-3)%p = 2*i-1 
        ocpkF(4*i-2)%p = 2*i-1
        ocpkF(4*i-1)%p = 2*i
        ocpkF(  4*i)%p = 2*i
    end do

    if (mod(nlat,4)/=0) then
        ocpkF(nlat-1)%g = nlat/4
        ocpkF(nlat  )%g = nlat-nlat/4-1
        ocpkF(nlat-1)%i = 1
        ocpkF(nlat  )%i = 2
        ocpkF(nlat-1)%p = nlat/2
        ocpkF(nlat  )%p = nlat/2
    end if

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

            ocpkP(1,2*i-1)%f = 4*i-3
            ocpkP(2,2*i-1)%f = 4*i-2
            ocpkP(1,  2*i)%f = 4*i-1
            ocpkP(2,  2*i)%f = 4*i
        end do
        if (mod(nlat,4)/=0) then
            ocpkP(1,nlat/2)%g = nlat/4
            ocpkP(2,nlat/2)%g = nlat-nlat/4-1
            ocpkP(1,nlat/2)%f = nlat-1
            ocpkP(2,nlat/2)%f = nlat
        end if
    else
        do i = 1, nlat/4
            ocpkP(1,4*i-3)%g = i
            ocpkP(1,4*i-2)%g = nlat/2-i+1
            ocpkP(1,4*i-1)%g = nlat-i+1
            ocpkP(1,  4*i)%g = nlat/2+i

            ocpkP(1,4*i-3)%f = 4*i-3  
            ocpkP(1,4*i-2)%f = 4*i-2
            ocpkP(1,4*i-1)%f = 4*i-1
            ocpkP(1,  4*i)%f = 4*i
        end do
        if (mod(nlat,4)/=0) then
            ocpkP(1,nlat-1)%g = nlat/4
            ocpkP(1,nlat  )%g = nlat-nlat/4-1
            ocpkP(1,nlat-1)%f = nlat-1
            ocpkP(1,nlat  )%f = nlat
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

    initialized = .true.

    return
end subroutine init_ocpack


function lpfac(ni) result(maxprime)
    implicit none
    integer, intent(in) :: ni
    integer :: maxprime
    integer :: n, i

    n = ni

    maxprime = -1

    do while (mod(n,2) == 0)
        maxprime = 2
        n = ishft(n,-1)
    enddo

    do i = 3, int(sqrt(real(n))), 2
        do while (mod(n,i) == 0)
            maxprime = i
            n = n / i
        enddo
    enddo

    if (n > 2) maxprime = n

    return
end function lpfac

end  module ocpack_mod

#ifdef test_ocpack
program test
use ocpack_mod
implicit none
integer, parameter :: nlat=94, nfour=62
type(ocpack_typeP), allocatable :: ocpkP(:,:)
type(ocpack_typeF), allocatable :: ocpkF(:)
integer :: i, j

!call init_ocpack(nlat,nfour,maxlon=192,ispacked=.false.,isreduced=.true.)
!call init_ocpack(nlat,nfour,ispacked=.false.,isreduced=.true.)
call init_ocpack(nlat,nfour,ispacked=.true.,isreduced=.true.)

allocate(ocpkF(oc_nlat()))
allocate(ocpkP(oc_npack(),oc_ny()))

call get_ocpackF(ocpkF)
call get_ocpackP(ocpkP)

do i = 1, oc_ny()
    do j = 1, oc_npack()
        print '(8(i4,2x))', i, j, ocpkP(j,i)%is, ocpkP(j,i)%ie, &
                      ocpkP(j,i)%fs, ocpkP(j,i)%fe, ocpkP(j,i)%g, ocpkP(j,i)%f
    end do
end do

    print *, 'ocpkF+++++++++'
do i = 1, oc_nlat()
    print '(5(i4,2x))', ocpkF(i)%p, ocpkF(i)%i, ocpkF(i)%flen, ocpkF(i)%g
end do

end program test
#endif
