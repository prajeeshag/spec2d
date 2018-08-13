module ocpack_mod

use mpp_mod, only : mpp_error, fatal, warning, note, mpp_init
use strman_mod, only : int2str

implicit none
private

public :: ocpack_type, init_ocpack, get_ocpack

type ocpack_type
    integer :: is
    integer :: ie
    integer :: ilen
    integer :: glat
end type 

type(ocpack_type), dimension(:,:), allocatable :: ocpk1
type(ocpack_type), dimension(:,:), allocatable :: ocpk2

integer :: nplon = 20

logical :: initialized=.false.

interface get_ocpack
    module procedure get_ocpack1, get_ocpack2
end interface get_ocpack

contains

subroutine get_ocpack1(ocpack) 
    type(ocpack_type), intent(out) :: ocpack(:)

    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')

    
    if(size(ocpack,1)==size(ocpk1,2)) then
        ocpack = ocpk1(1,:)
    else
        call mpp_error(fatal, 'get_ocpack: invalid size for argument "ocpack"!!!!') 
    end if

    return
end subroutine get_ocpack1

subroutine get_ocpack2(ocpack) 
    type(ocpack_type), intent(out) :: ocpack(:,:)

    if(.not.initialized) call mpp_error(fatal,'ocpack_mod: module not initiliazed!!!')

    
    if(size(ocpack,1)==1.and.size(ocpack,2)==size(ocpk1,2)) then
        ocpack = ocpk1
    elseif(size(ocpack,1)==2.and.size(ocpack,2)==size(ocpk2,2)) then
        ocpack = ocpk2
    else
        call mpp_error(fatal, 'get_ocpack: invalid size for argument "ocpack"!!!!') 
    end if

    return
end subroutine get_ocpack2


subroutine init_ocpack(nlat, maxlon, reduce)

    integer, intent(in) :: nlat
    integer, intent(in), optional :: maxlon
    logical, intent(in), optional :: reduce

    logical :: reduce1
    integer :: i, j, nplon_new
    integer, dimension(nlat) :: lonsperlat
    integer :: nlon, ocnx1, ocny1, ocnx2, ocny2


    call mpp_init()

    reduce1 = .true.

    if(present(reduce)) reduce1 = reduce

    if (mod(nlat,2).ne.0) then
        call mpp_error(fatal, 'ocpack_mod: nlat should be a multiple of 2')
    end if

    if (reduce1.and.present(maxlon)) then
        nplon_new = maxlon-(4*(nlat/2-1))
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

    nlon = maxval(lonsperlat)
  
    ocnx1 = nlon
    ocny1 = nlat

    ocnx2 = nplon+nlon
    ocny2 = nlat/2
   
    allocate(ocpk1(1,ocny1))
    allocate(ocpk2(2,ocny2))
 
    ocpk2(1,:)%is = 1
    ocpk2(2,:)%ie = ocnx2
    ocpk2(:,:)%glat = -1
  
    do i = 1, nlat/4
        ocpk2(1,2*i-1)%glat = i
        ocpk2(2,2*i-1)%glat = nlat/2-i+1
        ocpk2(1,2*i)%glat   = nlat-i+1
        ocpk2(2,2*i)%glat   = nlat/2+i
    end do

    if (mod(nlat,4)/=0) then
        ocpk2(1,nlat/2)%glat = nlat/4
        ocpk2(2,nlat/2)%glat = nlat-nlat/4-1
    end if

    do i = 1, ocny2
        ocpk2(1,i)%ie = ocpk2(1,i)%is + lonsperlat(ocpk2(1,i)%glat) - 1
        ocpk2(2,i)%is = ocpk2(1,i)%ie + 1
        ocpk2(1,i)%ilen = lonsperlat(ocpk2(1,i)%glat)
        ocpk2(2,i)%ilen = lonsperlat(ocpk2(2,i)%glat)
    end do

    
    ocpk1(1,:)%is = 1
    ocpk1(:,:)%glat = -1

    do i = 1, nlat/4
        ocpk1(1,4*i-3)%glat = i
        ocpk1(1,4*i-2)%glat = nlat/2-i+1
        ocpk1(1,4*i-1)%glat = nlat-i+1
        ocpk1(1,4*i)%glat   = nlat/2+i
    end do
    if (mod(nlat,4)/=0) then
        ocpk1(1,nlat-1)%glat = nlat/4
        ocpk1(1,nlat)%glat = nlat-nlat/4-1
    end if

    do i = 1, ocny1
        ocpk1(1,i)%ie = ocpk1(1,i)%is + lonsperlat(ocpk1(1,i)%glat) - 1
        ocpk1(1,i)%ilen = lonsperlat(ocpk1(1,i)%glat)
    end do

    initialized = .true.

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

type(ocpack_type), allocatable :: ocpk1(:,:), ocpk2(:,:)
integer :: ny=94

call init_ocpack(94)

allocate(ocpk1(1,ny))
allocate(ocpk2(2,ny/2))

call get_ocpack(ocpk1)
call get_ocpack(ocpk2)

do i = 1, 94/2
    print '(8(i4,2x))', ocpk2(1,i)%is, ocpk2(1,i)%ie, ocpk2(1,i)%ilen, &
                        ocpk2(1,i)%glat, ocpk2(2,i)%is, ocpk2(2,i)%ie,&
                        ocpk2(2,i)%ilen, ocpk2(2,i)%glat
end do

    print *, 'ocpk1+++++++++'
do i = 1, 94
    print '(4(i4,2x))', ocpk1(1,i)%is, ocpk1(1,i)%ie, ocpk1(1,i)%ilen, &
                        ocpk1(1,i)%glat
end do

end program test
#endif
