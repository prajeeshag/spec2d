
module fourier_spherical_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe

use mpp_domains_mod, only : domain2D

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

use gauss_and_legendre_mod, only : compute_legendre, compute_gaussian 

use spherical_mod, only : spherical_init

use spherical_data_mod, only : nlat, legendrePol, specVar
use spherical_data_mod, only : ewlen4m, ews4m, ewe4m, owlen4m, ows4m, owe4m
use spherical_data_mod, only : js, je, js_hem, je_hem, ms, me, nwaves_oe
use spherical_data_mod, only : tshuffle, ne4m, ns4m, iseven, num_fourier
use spherical_data_mod, only : num_spherical, trunc, init_spherical_data

implicit none
private

public :: init_fourier_spherical, fourier_to_spherical, spherical_to_fourier

real, allocatable, dimension(:) :: sin_lat
real, allocatable, dimension(:) :: cos_lat
real, allocatable, dimension(:) :: cosm_lat
real, allocatable, dimension(:) :: cosm2_lat
real, allocatable, dimension(:) :: deg_lat
real, allocatable, dimension(:) :: wts_lat
real, allocatable, dimension(:) :: sin_hem

type(legendrePol(js=:,je=:,n=:)), allocatable :: legendre, legendre_wts

logical :: debug = .false.
logical :: initialized = .false.
   
contains

!------------------------------------------------------------------------------
subroutine init_fourier_spherical(num_fourier_in, num_spherical_in, nlat_in, &
                nwaves_oe_out, domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    integer, intent(in) :: num_fourier_in, num_spherical_in, nlat_in
    type(domain2d), optional :: domain_fourier_in
    integer, optional :: tshuffle_in(0:num_fourier_in)
    integer, intent(out) :: nwaves_oe_out

    call init_spherical_data(num_fourier_in, num_spherical_in, nlat_in, &
                            nwaves_oe_out, domain_fourier_in, tshuffle_in)

    call define_gaussian

    call define_legendre

    call spherical_init()

    initialized = .true.

end subroutine init_fourier_spherical



!--------------------------------------------------------------------------------   
subroutine spherical_to_fourier(waves,fourier)
!--------------------------------------------------------------------------------   
    complex, intent(out) :: fourier(:,js:,ms:) ! lat, lev, fourier
    type(specVar(nlev=*,n=*)), intent(in) :: waves

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m

    if (.not.initialized) call mpp_error('spherical_to_fourier', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)

    do m = ms, me
        if(ewlen4m(m)<1) cycle
        ews = ews4m(m); ewe = ewe4m(m)
        call do_matmul(waves%ev(ks:ke,ews:ewe), &
                       legendre%ev(js_hem:je_hem,ews:ewe), &
                       even(ks:ke,js_hem:je_hem,m),'T')
    enddo

    do m = ms, me
        if(owlen4m(m)<1) cycle
        ows = ows4m(m); owe = owe4m(m)
        call do_matmul(waves%od(ks:ke,ows:owe), &
                       legendre%od(js_hem:je_hem,ows:owe), &
                       odd(ks:ke,js_hem:je_hem,m),'T')
    enddo 
   
    fourier(ks:ke,js+1:je:2,ms:me) = even(ks:ke,js_hem:je_hem,ms:me) + odd(ks:ke,js_hem:je_hem,ms:me)
    fourier(ks:ke,js:je:2,ms:me)   = even(ks:ke,js_hem:je_hem,ms:me) - odd(ks:ke,js_hem:je_hem,ms:me)

    !odd(js_hem:je_hem,ks:ke,ms:me) = fourier(js+1:je:2,ks:ke,ms:me) - fourier(js:je:2,:,:) ! north_hem - south_hem
    !even(js_hem:je_hem,ks:ke,ms:me) = fourier(js+1:je:2,ks:ke,ms:me) + fourier(js:je:2,ks:ke,ms:me) ! north_hem + south_hem

    return 
end subroutine spherical_to_fourier



!--------------------------------------------------------------------------------   
subroutine fourier_to_spherical(fourier, waves)
!--------------------------------------------------------------------------------   
    complex, intent(in) :: fourier(:,js:,ms:) ! lat, lev, fourier
    type(specVar(nlev=*,n=*)), intent(out) :: waves

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m, k

    if (.not.initialized) call mpp_error('fourier_to_spherical', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)

    odd(ks:ke,js_hem:je_hem,ms:me)  = fourier(ks:ke,js+1:je:2,ms:me) - fourier(ks:ke,js:je:2,ms:me) ! north_hem - south_hem
    even(ks:ke,js_hem:je_hem,ms:me) = fourier(ks:ke,js+1:je:2,ms:me) + fourier(ks:ke,js:je:2,ms:me) ! north_hem + south_hem

    do m = ms, me
        if(ewlen4m(m)<1) cycle
        ews = ews4m(m); ewe = ewe4m(m)
        call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                       legendre_wts%ev(js_hem:je_hem,ews:ewe), &
                       waves%ev(ks:ke,ews:ewe),'N')
    enddo

    do m = ms, me
        if(owlen4m(m)<1) cycle
        ows = ows4m(m); owe = owe4m(m)
        call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                       legendre_wts%od(js_hem:je_hem,ows:owe), &
                       waves%od(ks:ke,ows:owe),'N')
    enddo 

    return 
end subroutine fourier_to_spherical


!--------------------------------------------------------------------------------   
subroutine do_matmul(A,B,C,TRANSB)
!--------------------------------------------------------------------------------
    complex, intent(in) :: A(:,:) !k,m transa=T or m,k transa=N
    real, intent(in) :: B(:,:) !k,n
    complex, intent(inout) :: C(:,:) !m,n
    character, intent(in) :: TRANSB

    type(C_PTR) :: APTR, CPTR
    real, pointer :: AP(:,:), CP(:,:)

    character, parameter :: TRANSA='N'
    real, parameter :: ALPHA=1., BETA=0.
    integer :: M, K, N, LDA, LDB, LDC, i

    select case (TRANSB)
    case('N','n')
        N = size(B,2)
    case('T','t')
        N = size(B,1)
    case default
        call mpp_error('do_matmul', 'TRANSB should be either T or N', FATAL)
    end select

    M=size(A,1)*2; K=size(A,2)

    APTR = C_LOC(A)
    call c_f_pointer(APTR, AP, [M,K])

    CPTR = C_LOC(C)
    call c_f_pointer(CPTR, CP, [M,N])

    LDA=size(AP,1); LDB=size(B,1); LDC=size(CP,1)

    call dgemm(TRANSA,TRANSB,M,N,K,ALPHA,AP(:,:),LDA, &
               B(:,:),LDB,BETA,CP(:,:),LDC)
    return
end subroutine do_matmul


!------------------------------------------------------------------------------
subroutine define_gaussian
!------------------------------------------------------------------------------

    integer :: j
    real, dimension(nlat/2) :: wts_hem

    allocate (sin_lat(js:je))
    allocate (cos_lat(js:je))
    allocate (cosm_lat(js:je))
    allocate (cosm2_lat(js:je))
    allocate (wts_lat(js:je))
    allocate (deg_lat(js:je))
    allocate (sin_hem(nlat/2))

    call compute_gaussian(sin_hem, wts_hem, nlat/2)

    sin_lat(js:je:2) = -sin_hem !Southern hemisphere
    sin_lat(js+1:je:2) = sin_hem !Northern hemisphere

    wts_lat(js:je:2) = wts_hem
    wts_lat(js+1:je:2) = wts_hem

    cos_lat = sqrt(1-sin_lat*sin_lat)
    cosm_lat = 1./cos_lat
    cosm2_lat = 1./(cos_lat*cos_lat)
    deg_lat = asin(sin_lat)*180.0/pi

    return
end subroutine define_gaussian

!--------------------------------------------------------------------------------   
subroutine define_legendre
!--------------------------------------------------------------------------------   
    integer :: j, m, w, wo, we, mshuff, n
    real, dimension(0:num_fourier,0:num_spherical,nlat/2) :: legendre_global
    character(len=8) :: suffix

    allocate(legendrePol(js=js_hem,je=je_hem,n=nwaves_oe) :: legendre, legendre_wts)

    call compute_legendre(legendre_global, num_fourier, 1, num_spherical, sin_hem, nlat/2)

    if (debug) then
        if (mpp_pe()==mpp_root_pe()) then
            write(suffix,'(I4.4)') mpp_pe()
            call write_data('debug_fourier_spherical_'//trim(suffix),'glegen', &
            reshape(legendre_global,shape=[nlat/2,num_spherical+1,num_fourier+1], order=[3,2,1]),no_domain=.true.)
        endif
    endif

    do j = js_hem, je_hem
        w = 0
        wo = 0
        we = 0
        do m = ms, me
            mshuff = tshuffle(m)
            do n = ns4m(m), ne4m(m)
                w = w + 1
                if (iseven(w)) then
                    we = we + 1
                    legendre%ev(j,we) = legendre_global(mshuff,n,j)
                    if (mshuff+n>num_fourier) legendre%ev(j,we) = 0.
                else
                    wo = wo + 1
                    legendre%od(j,wo) = legendre_global(mshuff,n,j)
                    if (mshuff+n>num_fourier) legendre%od(j,wo) = 0.
                endif
            enddo
        enddo

        legendre_wts%ev(j,:) = legendre%ev(j,:)*wts_lat(2*j)
        legendre_wts%od(j,:) = legendre%od(j,:)*wts_lat(2*j)
    enddo

    return
end subroutine define_legendre

end module fourier_spherical_mod

