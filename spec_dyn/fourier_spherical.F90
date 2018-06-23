
module fourier_spherical_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe

use mpp_domains_mod, only : domain2D

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

use spherical_mod, only : spherical_init, legendre_wts, legendre, legendredphi, triangle_mask

use spherical_data_mod, only : nlat, legendrePol, specVar
use spherical_data_mod, only : ewlen4m, ews4m, ewe4m, owlen4m, ows4m, owe4m
use spherical_data_mod, only : js, je, js_hem, je_hem, ms, me, nwaves_oe
use spherical_data_mod, only : tshuffle, ne4m, ns4m, iseven, num_fourier
use spherical_data_mod, only : num_spherical, trunc, init_spherical_data

implicit none
private

public :: init_fourier_spherical, fourier_to_spherical, spherical_to_fourier

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

    call spherical_init()

    initialized = .true.

end subroutine init_fourier_spherical


!--------------------------------------------------------------------------------   
subroutine spherical_to_fourier(waves,fourier,lat_deriv)
!--------------------------------------------------------------------------------   
    complex, intent(out) :: fourier(:,js:,ms:) ! lat, lev, fourier
    type(specVar(nlev=*,n=*)), intent(in) :: waves
    logical, intent(in), optional :: lat_deriv

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m, nj
    logical :: deriv

    deriv = .false.
    
    if (present(lat_deriv)) deriv = lat_deriv



    if (.not.initialized) call mpp_error('spherical_to_fourier', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)
    nj = legendre%nj

    if (deriv) then
        do m = ms, me
            if(ewlen4m(m)<1) cycle
            ews = ews4m(m); ewe = ewe4m(m)
            call do_matmul(waves%ev(ks:ke,ews:ewe), &
                           legendredphi%ev(:,ews:ewe), &
                           even(ks:ke,js_hem:je_hem,m),'T')
        enddo

        do m = ms, me
            if(owlen4m(m)<1) cycle
            ows = ows4m(m); owe = owe4m(m)
            call do_matmul(waves%od(ks:ke,ows:owe), &
                           legendredphi%od(:,ows:owe), &
                           odd(ks:ke,js_hem:je_hem,m),'T')
        enddo 

        fourier(ks:ke,js+1:je:2,ms:me) = even(ks:ke,js_hem:je_hem,ms:me) + odd(ks:ke,js_hem:je_hem,ms:me)
        fourier(ks:ke,js:je:2,ms:me)   = odd(ks:ke,js_hem:je_hem,ms:me)  - even(ks:ke,js_hem:je_hem,ms:me)

    else
        do m = ms, me
            if(ewlen4m(m)<1) cycle
            ews = ews4m(m); ewe = ewe4m(m)
            call do_matmul(waves%ev(ks:ke,ews:ewe), &
                           legendre%ev(1:nj,ews:ewe), &
                           even(ks:ke,js_hem:je_hem,m),'T')
        enddo

        do m = ms, me
            if(owlen4m(m)<1) cycle
            ows = ows4m(m); owe = owe4m(m)
            call do_matmul(waves%od(ks:ke,ows:owe), &
                           legendre%od(1:nj,ows:owe), &
                           odd(ks:ke,js_hem:je_hem,m),'T')
        enddo 

        fourier(ks:ke,js+1:je:2,ms:me) = even(ks:ke,js_hem:je_hem,ms:me) + odd(ks:ke,js_hem:je_hem,ms:me)
        fourier(ks:ke,js:je:2,ms:me)   = even(ks:ke,js_hem:je_hem,ms:me) - odd(ks:ke,js_hem:je_hem,ms:me)

    endif
   
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

    integer :: ks, ke, ews, ewe, ows, owe, m, k, nj

    if (.not.initialized) call mpp_error('fourier_to_spherical', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)
    nj = legendre_wts%nj

    odd(ks:ke,js_hem:je_hem,ms:me)  = fourier(ks:ke,js+1:je:2,ms:me) - fourier(ks:ke,js:je:2,ms:me) ! north_hem - south_hem
    even(ks:ke,js_hem:je_hem,ms:me) = fourier(ks:ke,js+1:je:2,ms:me) + fourier(ks:ke,js:je:2,ms:me) ! north_hem + south_hem

    do m = ms, me
        if(ewlen4m(m)<1) cycle
        ews = ews4m(m); ewe = ewe4m(m)
        call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                       legendre_wts%ev(1:nj,ews:ewe), &
                       waves%ev(ks:ke,ews:ewe),'N')
    enddo

    do m = ms, me
        if(owlen4m(m)<1) cycle
        ows = ows4m(m); owe = owe4m(m)
        call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                       legendre_wts%od(1:nj,ows:owe), &
                       waves%od(ks:ke,ows:owe),'N')
    enddo 

    do k = ks, ke
        waves%ev(k,:) = waves%ev(k,:)*triangle_mask%ev(:)
        waves%od(k,:) = waves%od(k,:)*triangle_mask%od(:)
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



end module fourier_spherical_mod

