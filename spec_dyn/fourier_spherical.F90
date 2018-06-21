
module fourier_spherical_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_layout

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

use gauss_and_legendre_mod, only : compute_legendre, compute_gaussian 

implicit none
private

public :: init_fourier_spherical, fourier_to_spherical, spherical_to_fourier
public :: specVar

integer :: num_fourier
integer :: num_spherical
integer :: nlat
integer :: ms, me, mlen
integer :: js_hem, je_hem, jlen_hem, js, je, jlen

integer, allocatable :: ns4m(:) ! Starting of spherical for a particular fourier
integer, allocatable :: ne4m(:) ! Ending of spherical for a particulat fourier 
                                !(=num_spherical in case of Triangular truncation)
integer, allocatable :: nlen4m(:) ! number of spherical for a particular fourier
                                  !(Constant in case of rhomboidal truncation)

integer :: nwaves !Total number of spectral waves (local)
integer :: nwaves_oe !Total number of odd-even waves [=nwaves/2])

logical, allocatable :: iseven(:) ! oddeven flag (.true. = even); (.false. = odd)

integer, allocatable :: ws4m(:), we4m(:), wlen4m(:) !starting and ending index of waves for a particular m
integer, allocatable :: ews4m(:), ewe4m(:), ewlen4m(:) !starting and ending index of even waves for a particular m
integer, allocatable :: ows4m(:), owe4m(:), owlen4m(:) !starting and ending index of odd waves for a particular m

integer, allocatable :: tshuffle(:)

type(domain2d) :: domain_spherical

real, allocatable, dimension(:) :: sin_lat
real, allocatable, dimension(:) :: cos_lat
real, allocatable, dimension(:) :: cosm_lat
real, allocatable, dimension(:) :: cosm2_lat
real, allocatable, dimension(:) :: deg_lat
real, allocatable, dimension(:) :: wts_lat
real, allocatable, dimension(:) :: sin_hem

type legendrePol(js,je,n)
    integer, len :: js, je, n
    real, dimension(js:je,n) :: ev
    real, dimension(js:je,n) :: od
end type legendrePol

type specCoef(n)
    integer, len :: n
    real, dimension(n) :: ev
    real, dimension(n) :: od
end type specCoef

type specVar(n,nlev)
    integer, len :: n, nlev
    complex, dimension(nlev,n) :: ev
    complex, dimension(nlev,n) :: od
end type specVar

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

    integer :: m, w, n, unit, iostat, neadj
    integer :: noddwaves, nevenwaves
    character (len=8) :: suffix

    namelist/fourier_spherical_nml/debug

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=fourier_spherical_nml,iostat=iostat)
    call close_file(unit)

    if (mod(nlat_in,2)/=0) &
        call mpp_error('fourier_spherical', 'NLAT should be a even number', FATAL)
    
    nlat = nlat_in
    num_fourier = num_fourier_in
    num_spherical = num_spherical_in
    if (mod(num_spherical,2)==0) num_spherical = num_spherical + 1

    if (present(domain_fourier_in)) then
        call mpp_get_compute_domain(domain_fourier_in,js,je,ms,me)
    else
        js = 1
        je = nlat
        ms = 0
        me = num_fourier
    endif

    jlen = je - js + 1

    if (mod(js,2)==0) call mpp_error('fourier_spherical', 'js should be a odd number!!!', FATAL)
    if (mod(jlen,2)/=0) call mpp_error('fourier_spherical', 'jlen should be a even number!!!', FATAL)
    if (mod(je,2)/=0) call mpp_error('fourier_spherical', 'je should be a even number!!!', FATAL)

    js_hem = js/2 + 1
    je_hem = (je-1)/2 + 1
    jlen_hem = jlen/2

    mlen = me-ms+1

    allocate(tshuffle(ms:me))

    if(present(tshuffle_in)) then
        tshuffle(ms:me) = tshuffle_in(ms:me)
    else
        forall(m=ms:me) tshuffle(m) = m
    endif

    allocate(ns4m(ms:me))
    allocate(ne4m(ms:me))
    allocate(nlen4m(ms:me))
    allocate(ws4m(ms:me))
    allocate(we4m(ms:me))
    allocate(wlen4m(ms:me))
    allocate(ews4m(ms:me))
    allocate(ewe4m(ms:me))
    allocate(ewlen4m(ms:me))
    allocate(ows4m(ms:me))
    allocate(owe4m(ms:me))
    allocate(owlen4m(ms:me))

    ns4m(:) = 0
    do m = ms, me
       neadj = num_spherical - tshuffle(m)
       if (mod(neadj,2)==0) neadj = neadj + 1 
       ne4m(m) = neadj
    enddo 

    nlen4m(:) = ne4m - ns4m + 1

    nwaves = 0; noddwaves = 0; nevenwaves = 0

    do m = ms, me
        nwaves = nwaves + nlen4m(m)
        we4m(m) = nwaves
    enddo
    
    allocate(iseven(nwaves))
    iseven(:) = .false.

    w = 0
    do m = ms, me
        do n = ns4m(m), ne4m(m)
            w = w + 1
            if (mod(n,2)==0) then
                nevenwaves = nevenwaves + 1
                iseven(w) = .true.
            else
                noddwaves = noddwaves + 1
                iseven(w) = .false.
            endif 
        enddo
        ewe4m(m) = nevenwaves
        owe4m(m) = noddwaves 
    enddo

    if (noddwaves/=nevenwaves) call mpp_error('fourier_spherical', 'noddwaves/=nevenwaves', FATAL)

    nwaves_oe = noddwaves
    nwaves_oe_out = nwaves_oe

    ws4m(ms) = 1
    ews4m(ms) = 1
    ows4m(ms) = 1
    
    do m = ms+1, me
        ws4m(m) = we4m(m-1)+1
        ews4m(m) = ewe4m(m-1)+1
        ows4m(m) = owe4m(m-1)+1
    enddo

    wlen4m(:) = we4m(:) - ws4m(:) + 1
    ewlen4m(:) = ewe4m(:) - ews4m(:) + 1
    owlen4m(:) = owe4m(:) - ows4m(:) + 1

    if (debug) then
        print *, 'ews4m=', ews4m 
        print *, 'ows4m=', ows4m 
        print *, 'ewlen4m=', ewlen4m 
        print *, 'owlen4m=', owlen4m 
    endif

    call define_gaussian

    call define_legendre

    if (debug) then
        write(suffix,'(I4.4)') mpp_pe()
        print *, 'debug from fourier_spherical, pe= ', trim(suffix)
        call write_data('debug_fourier_spherical_'//trim(suffix),'olegen',legendre%od,no_domain=.true.)
        call write_data('debug_fourier_spherical_'//trim(suffix),'elegen',legendre%ev,no_domain=.true.)
        print *, 'pe, noddwaves, nevenwaves =', mpp_pe(), noddwaves, nevenwaves
        print *, 'pe, ns4m(:)=', ns4m(:)
        print *, 'pe, ne4m(:)=', ne4m(:)
    endif

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

