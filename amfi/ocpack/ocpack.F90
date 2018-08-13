module ocpack_mod

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init
use strman_mod, only : int2str

implicit none
private

public :: ocpack_type, init_ocpack, get_ocpack_info, get_ocpack

type ocpack_type
    integer :: is
    integer :: ie
    integer :: ilen
    integer :: glat
end type 

type(ocpack_type), dimension(:,:), allocatable :: OCPK

integer :: NLON, NLAT, OCNX, OCNY
integer :: NPLON = 20, NPACK=2
integer, dimension(:), allocatable :: LONSPERLAT

logical :: initialized=.false.

contains

subroutine get_ocpack_info(ny, np, nx)
    integer, intent(out) :: ny, np
    integer, intent(out), optional :: nx

    if(.not.initialized) call mpp_error(FATAL,'ocpack_mod: module not initiliazed!!!')

    ny = OCNY
    np = NPACK
    if(present(nx)) nx = OCNX

    return
end subroutine get_ocpack_info

subroutine get_ocpack(ocpack) 
    type(ocpack_type), intent(out) :: ocpack(:,:)

    if(.not.initialized) call mpp_error(FATAL,'ocpack_mod: module not initiliazed!!!')

    if(size(ocpack,1)/=size(OCPK,1).or.size(ocpack,2)/=size(OCPK,2)) then
        call mpp_error(FATAL, 'get_ocpack: argument size mismatch!!!!') 
    end if

    ocpack = OCPK

end subroutine get_ocpack

subroutine init_ocpack(num_lat, maxlon, num_pack, reduce)

    integer, intent(in) :: num_lat
    integer, intent(in), optional :: maxlon
    integer, intent(in), optional :: num_pack
    logical, intent(in), optional :: reduce

    logical :: reduce1
    integer :: i, j, nplon_new

    call mpp_init()

    reduce1 = .true.

    if(present(reduce)) reduce1 = reduce

    NLAT = num_lat

    if (mod(NLAT,2).ne.0) then
        call mpp_error(FATAL, 'ocpack_mod: NLAT should be a multiple of 2')
    end if

    if (reduce1.and.present(maxlon)) then
        nplon_new = maxlon-(4*(NLAT/2-1))
        if (nplon_new/=NPLON) then
            NPLON = nplon_new
            call mpp_error(WARNING, 'Changing default number of (20) pole '// &
                      'longitudes (NPLON) to: '//int2str(NPLON))
        end if
    end if

    if(present(num_pack)) then
        if (num_pack<1.or.num_pack>2) call mpp_error(FATAL, 'num_pack should be either 1 or 2')
        NPACK = num_pack
    end if

    allocate(LONSPERLAT(NLAT)) 

    LONSPERLAT(1) = NPLON
    LONSPERLAT(NLAT) = NPLON

    do i = 2, NLAT/2
        LONSPERLAT(i) = LONSPERLAT(i-1) + 4
        LONSPERLAT(NLAT-i+1) = LONSPERLAT(i)
    end do

    NLON = maxval(LONSPERLAT)
  
    OCNX = NPLON*(NPACK-1)+NLON
    OCNY = NLAT/NPACK
   
    allocate(OCPK(NPACK,OCNY))
 
    if (NPACK==2) then 

        OCPK(1,:)%IS = 1
        OCPK(2,:)%IE = OCNX
        OCPK(:,:)%GLAT = -1
  
        do i = 1, NLAT/4
            OCPK(1,2*i-1)%GLAT = i
            OCPK(2,2*i-1)%GLAT = NLAT/2-i+1
            OCPK(1,2*i)%GLAT   = NLAT-i+1
            OCPK(2,2*i)%GLAT   = NLAT/2+i
        end do

        if (mod(NLAT,4)/=0) then
            OCPK(1,NLAT/2)%GLAT = NLAT/4
            OCPK(2,NLAT/2)%GLAT = NLAT-NLAT/4-1
        end if

        do i = 1, OCNY
            OCPK(1,i)%IE = OCPK(1,i)%IS + LONSPERLAT(OCPK(1,i)%GLAT) - 1
            OCPK(2,i)%IS = OCPK(1,i)%IE + 1
            OCPK(1,i)%ILEN = LONSPERLAT(OCPK(1,i)%GLAT)
            OCPK(2,i)%ILEN = LONSPERLAT(OCPK(2,i)%GLAT)
        end do
    else

        OCPK(1,:)%IS = 1
        OCPK(:,:)%GLAT = -1

        do i = 1, NLAT/4
            OCPK(1,4*i-3)%GLAT = i
            OCPK(1,4*i-2)%GLAT = NLAT/2-i+1
            OCPK(1,4*i-1)%GLAT = NLAT-i+1
            OCPK(1,4*i)%GLAT   = NLAT/2+i
        end do

        do i = 1, OCNY
            OCPK(1,i)%IE = OCPK(1,i)%IS + LONSPERLAT(OCPK(1,i)%GLAT) - 1
            OCPK(1,i)%ILEN = LONSPERLAT(OCPK(1,i)%GLAT)
        end do

    end if

    initialized = .true.

end subroutine init_ocpack


function lpfac(ni) result(maxPrime)
    implicit none
    integer, intent(in) :: ni
    integer :: maxPrime
    integer :: n, i

    n = ni

    maxPrime = -1

    do while (mod(n,2) == 0)
        maxPrime = 2
        n = ishft(n,-1)
    enddo

    do i = 3, int(sqrt(real(n))), 2
        do while (mod(n,i) == 0)
            maxPrime = i
            n = n / i
        enddo
    enddo

    if (n > 2) maxPrime = n

    return
end function lpfac

end  module ocpack_mod

#ifdef test_ocpack
program test
use ocpack_mod

type(ocpack_type), allocatable :: ocpk(:,:)
integer :: ny, np

call init_ocpack(94,num_pack=1)
!call init_ocpack(94)

call get_ocpack_info(ny,np)

allocate(ocpk(np,ny))

call get_ocpack(ocpk)
do i = 1, ny
    print '(8(I4,2x))', ocpk(1,i)%is, ocpk(1,i)%ie, ocpk(1,i)%ilen, &
                    ocpk(1,i)%glat, ocpk(2,i)%is, ocpk(2,i)%ie,&
                    ocpk(2,i)%ilen, ocpk(2,i)%glat
end do

end program test
#endif
