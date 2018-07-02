
module fourier_spherical_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe

use mpp_domains_mod, only : domain2D

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

use spherical_mod, only : spherical_init, Pnm_wts, Pnm, Hnm, Hnm_wts, do_truncation

use spherical_data_mod, only : wlen4m, ws4m, we4m
use spherical_data_mod, only : js, je, js_hem, je_hem, ms, me, nwaves_oe
use spherical_data_mod, only : tshuffle, ne4m, ns4m, iseven
use spherical_data_mod, only : init_spherical_data
use spherical_data_mod, only : nlat, ev, od, jlen_hem

implicit none
private

public :: init_fourier_spherical, fourier_to_spherical, spherical_to_fourier

logical :: debug = .false.
logical :: initialized = .false.
   
contains

!------------------------------------------------------------------------------
subroutine init_fourier_spherical(trunc_in, nlat_in, &
                nwaves_oe_out, domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    integer, intent(in) :: trunc_in, nlat_in
    type(domain2d), optional :: domain_fourier_in
    integer, optional :: tshuffle_in(0:)
    integer, intent(out) :: nwaves_oe_out

    call init_spherical_data(trunc_in, nlat_in, &
              nwaves_oe_out, domain_fourier_in, tshuffle_in)

    call spherical_init()

    initialized = .true.

end subroutine init_fourier_spherical


!--------------------------------------------------------------------------------   
subroutine spherical_to_fourier(waves,fourier,useHnm)
!--------------------------------------------------------------------------------   
    complex, intent(out) :: fourier(:,js:,ms:) ! lat, lev, fourier
    complex,dimension(:,:,:), intent(in) :: waves
    logical, intent(in), optional :: useHnm

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m, nj
    logical :: deriv

    deriv = .false.
    
    if (present(useHnm)) deriv = useHnm

    if (.not.initialized) call mpp_error('spherical_to_fourier', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)

    if (deriv) then
        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(waves(ks:ke,ews:ewe,ev), &
                           Hnm(1:jlen_hem,ews:ewe,ev), &
                           even(ks:ke,js_hem:je_hem,m),'T')
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(waves(ks:ke,ows:owe,od), &
                           Hnm(1:jlen_hem,ows:owe,od), &
                           odd(ks:ke,js_hem:je_hem,m),'T')
        enddo 

        fourier(ks:ke,js+1:je:2,ms:me) = even(ks:ke,js_hem:je_hem,ms:me) + odd(ks:ke,js_hem:je_hem,ms:me)
        fourier(ks:ke,js:je:2,ms:me)   = odd(ks:ke,js_hem:je_hem,ms:me)  - even(ks:ke,js_hem:je_hem,ms:me)

    else
        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(waves(ks:ke,ews:ewe,ev), &
                           Pnm(1:jlen_hem,ews:ewe,ev), &
                           even(ks:ke,js_hem:je_hem,m),'T')
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(waves(ks:ke,ows:owe,od), &
                           Pnm(1:jlen_hem,ows:owe,od), &
                           odd(ks:ke,js_hem:je_hem,m),'T')
        enddo 

        fourier(ks:ke,js+1:je:2,ms:me) = even(ks:ke,js_hem:je_hem,ms:me) + odd(ks:ke,js_hem:je_hem,ms:me)
        fourier(ks:ke,js:je:2,ms:me)   = even(ks:ke,js_hem:je_hem,ms:me) - odd(ks:ke,js_hem:je_hem,ms:me)

    endif
   
    return 
end subroutine spherical_to_fourier


!--------------------------------------------------------------------------------   
subroutine fourier_to_spherical(fourier, waves, useHnm, do_trunc)
!--------------------------------------------------------------------------------   
    complex, intent(in) :: fourier(:,js:,ms:) ! lat, lev, fourier
    complex,dimension(:,:,:), intent(out) :: waves
    logical, optional :: useHnm, do_trunc

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    logical :: useHnm1, do_trunc1
    integer :: ks, ke, ews, ewe, ows, owe, m, k, nj

    useHnm1 = .false.
    if(present(useHnm)) useHnm1=useHnm

    do_trunc1 = .true.
    if(present(do_trunc)) do_trunc1=do_trunc

    if (.not.initialized) call mpp_error('fourier_to_spherical', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,1)

    if (useHnm1) then
        odd(ks:ke,js_hem:je_hem,ms:me)  = fourier(ks:ke,js+1:je:2,ms:me) + fourier(ks:ke,js:je:2,ms:me) ! 
        even(ks:ke,js_hem:je_hem,ms:me) = fourier(ks:ke,js+1:je:2,ms:me) - fourier(ks:ke,js:je:2,ms:me) !

        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                           Hnm_wts(1:jlen_hem,ews:ewe,ev), &
                           waves(ks:ke,ews:ewe,ev),'N')
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                           Hnm_wts(1:jlen_hem,ows:owe,od), &
                           waves(ks:ke,ows:owe,od),'N')
        enddo 

    else
        odd(ks:ke,js_hem:je_hem,ms:me)  = fourier(ks:ke,js+1:je:2,ms:me) - fourier(ks:ke,js:je:2,ms:me) ! 
        even(ks:ke,js_hem:je_hem,ms:me) = fourier(ks:ke,js+1:je:2,ms:me) + fourier(ks:ke,js:je:2,ms:me) ! 

        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                           Pnm_wts(1:jlen_hem,ews:ewe,ev), &
                           waves(ks:ke,ews:ewe,ev),'N')
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                           Pnm_wts(1:jlen_hem,ows:owe,od), &
                           waves(ks:ke,ows:owe,od),'N')
        enddo 
    endif

    call do_truncation(waves,do_trunc1)

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

