
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

public :: init_fourier_spherical, fourier_to_spherical

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

integer :: nwaves ! number of spectral waves (local)
integer :: noddwaves ! number of odd waves (local)
integer :: nevenwaves ! number of even waves (local)
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

real, allocatable, dimension(:,:) :: legendre ! nwaves, js_hem:he_hem 
real, allocatable, dimension(:,:) :: legendre_wts ! nwaves, js_hem:he_hem 

real, allocatable, dimension(:,:) :: elegendre ! 2, js_hem:je_hem nevenwaves
real, allocatable, dimension(:,:) :: elegendre_wts ! 2, js_hem:je_hem nevenwaves

real, allocatable, dimension(:,:) :: olegendre ! 2, js_hem:je_hem, noddwaves
real, allocatable, dimension(:,:) :: olegendre_wts ! 2, js_hem:je_hem, noddwaves

logical :: debug = .false.
logical :: initialized = .false.
   
contains

!------------------------------------------------------------------------------
subroutine init_fourier_spherical(num_fourier_in, num_spherical_in, nlat_in, &
                domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    integer, intent(in) :: num_fourier_in, num_spherical_in, nlat_in
    type(domain2d), optional :: domain_fourier_in
    integer, optional :: tshuffle_in(0:num_fourier_in)

    integer :: m, w, n, unit, iostat, neadj
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
        call write_data('debug_fourier_spherical_'//trim(suffix),'olegen',olegendre,no_domain=.true.)
        call write_data('debug_fourier_spherical_'//trim(suffix),'elegen',elegendre,no_domain=.true.)
        call write_data('debug_fourier_spherical_'//trim(suffix),'legen',legendre,no_domain=.true.)
        print *, 'pe, noddwaves, nevenwaves =', mpp_pe(), noddwaves, nevenwaves
        print *, 'pe, ns4m(:)=', ns4m(:)
        print *, 'pe, ne4m(:)=', ne4m(:)
    endif

    initialized = .true.

end subroutine init_fourier_spherical


!--------------------------------------------------------------------------------   
subroutine fourier_to_spherical(fourier, ewaves, owaves)
!--------------------------------------------------------------------------------   
    complex, intent(in) :: fourier(js:,:,ms:) ! lat, lev, fourier
    complex, intent(out) :: ewaves(:,:) ! nevenwaves, lev
    complex, intent(out) :: owaves(:,:) ! noddwaves, lev

    complex :: odd(js_hem:je_hem,size(fourier,2),ms:me)
    complex :: even(js_hem:je_hem,size(fourier,2),ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m

    if (.not.initialized) call mpp_error('fourier_to_spherical', 'call init_fourier_spherical first', FATAL)

    ks = 1; ke = size(fourier,2)

    odd(js_hem:je_hem,ks:ke,ms:me) = fourier(js+1:je:2,ks:ke,ms:me) - fourier(js:je:2,:,:) ! north_hem - south_hem
    even(js_hem:je_hem,ks:ke,ms:me) = fourier(js+1:je:2,ks:ke,ms:me) + fourier(js:je:2,ks:ke,ms:me) ! north_hem + south_hem

    do m = ms, me
        if(ewlen4m(m)<1) cycle
        ews = ews4m(m); ewe = ewe4m(m)
        call do_matmul(elegendre_wts(ews:ewe,js_hem:je_hem), &
                       even(js_hem:je_hem,ks:ke,m), &
                       ewaves(ews:ewe,ks:ke))
    enddo

    do m = ms, me
        if(owlen4m(m)<1) cycle
        ows = ows4m(m); owe = owe4m(m)
        call do_matmul(olegendre_wts(ows:owe,js_hem:je_hem), &
                       odd(js_hem:je_hem,ks:ke,m), &
                       owaves(ows:owe,ks:ke))
    enddo 
   
    return 
end subroutine fourier_to_spherical


!--------------------------------------------------------------------------------   
subroutine do_matmul(A,B,C)
!--------------------------------------------------------------------------------
    real, intent(in) :: A(:,:) !k,m
    complex, intent(in) :: B(:,:) !k,n
    complex, intent(inout) :: C(:,:) !k,n

    type(C_PTR) :: BPTR, CPTR
    real, pointer :: BP(:,:,:), CP(:,:,:)

    character :: TRANSA='T'
    character :: TRANSB='N'
    integer :: M, K, N, LDA, LDB, LDC, i
    real, parameter :: ALPHA=1., BETA=0.

    M=size(A,2); K=size(B,1); N=size(B,2)

    LDA=K; LDB=K; LDC=M

    BPTR = C_LOC(B)
    call c_f_pointer(BPTR, BP, [2,K,N])

    CPTR = C_LOC(C)
    call c_f_pointer(CPTR, CP, [2,M,N])

    do i = 1, 2
        call cgemm(TRANSA,TRANSB,M,N,K,ALPHA,A(1:K,1:M),LDA, &
                   BP(i,1:K,1:N),LDB,BETA,CP(i,1:M,1:N),LDC)
    enddo

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

    allocate(legendre(js_hem:je_hem,nwaves))
    allocate(legendre_wts(js_hem:je_hem,nwaves))

    allocate(olegendre(js_hem:je_hem,noddwaves))
    allocate(olegendre_wts(js_hem:je_hem,noddwaves))

    allocate(elegendre(js_hem:je_hem,nevenwaves))
    allocate(elegendre_wts(js_hem:je_hem,nevenwaves))

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
                legendre(j,w) = legendre_global(mshuff,n,j)
                if (iseven(w)) then
                    we = we + 1
                    elegendre(j,we) = legendre_global(mshuff,n,j)
                else
                    wo = wo + 1
                    olegendre(j,wo) = legendre_global(mshuff,n,j)
                endif
            enddo
        enddo

        legendre_wts(j,:) = legendre(j,:)*wts_lat(j)
        elegendre_wts(j,:) = elegendre(j,:)*wts_lat(j)
        olegendre_wts(j,:) = olegendre(j,:)*wts_lat(j)
    enddo

    return
end subroutine define_legendre

end module fourier_spherical_mod

