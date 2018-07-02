module spherical_data_mod

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe, mpp_sync

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_layout

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

implicit none
private

public :: init_spherical_data, get_wdecomp

integer, public :: num_fourier
integer, public :: num_spherical
integer, public :: trunc
integer, public :: truncadj
integer, public :: nlat
integer, public :: ms, me, mlen
integer, public :: js_hem, je_hem, jlen_hem, js, je, jlen

integer, public, allocatable :: ns4m(:) ! Starting of spherical for a particular fourier
integer, public, allocatable :: ne4m(:) ! Ending of spherical for a particulat fourier 
                                !(=num_spherical in case of Triangular
                                !truncation)
integer, public, allocatable :: nlen4m(:) ! number of spherical for a particular fourier
                                  !(Constant in case of rhomboidal truncation)

integer :: nwaves !Total number of waves
integer, public :: nwaves_oe !Total number of odd-even waves [=nwaves/2])
integer, public :: noddwaves, nevenwaves
integer :: noddwaves_g, nevenwaves_g, nwaves_oe_g

logical, public, allocatable :: iseven(:) ! oddeven flag (.true. = even); (.false. = odd)

integer, public, allocatable :: ws4m(:,:), we4m(:,:), wlen4m(:,:) !starting and ending index of 
                                                                  !even & odd waves for a particular m
integer, allocatable :: wdecomp(:,:)

integer, public, allocatable :: tshuffle(:)

logical :: debug

integer, parameter, public :: ev=1, od=2

contains

!--------------------------------------------------------------------------------   
subroutine get_wdecomp(wdom,neven_global,nodd_global)
!--------------------------------------------------------------------------------
    integer, intent(out) :: wdom(:,:)
    integer, intent(out), optional :: neven_global, nodd_global
    
    wdom(:,:) = wdecomp(:,:)

    if (present(neven_global)) neven_global = nevenwaves_g
    if (present(nodd_global)) nodd_global = noddwaves_g
end subroutine get_wdecomp

!------------------------------------------------------------------------------
subroutine init_spherical_data(trunc_in, nlat_in, &
                nwaves_oe_out, domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: trunc_in, nlat_in
    type(domain2d), optional :: domain_fourier_in
    integer, optional :: tshuffle_in(0:)
    integer, intent(out) :: nwaves_oe_out

    integer :: m, w, n, unit, iostat, neadj
    integer :: nsf4m(0:trunc_in), nef4m(0:trunc_in), nlenf4m(0:trunc_in)
    integer, allocatable :: wsf4m(:,:), wef4m(:,:), wlenf4m(:,:)
    integer :: we, wo
    character (len=8) :: suffix

    namelist/spherical_nml/debug

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=spherical_nml,iostat=iostat)
    call close_file(unit)

    if (mod(nlat_in,2)/=0) &
        call mpp_error('init_spherical_data', 'NLAT should be a even number', FATAL)
    
    nlat = nlat_in
    num_fourier = trunc_in
    num_spherical = trunc_in + 1
    trunc = trunc_in
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

    if (mod(js,2)==0) call mpp_error('init_spherical_data', 'js should be a odd number!!!', FATAL)
    if (mod(jlen,2)/=0) call mpp_error('init_spherical_data', 'jlen should be a even number!!!', FATAL)
    if (mod(je,2)/=0) call mpp_error('init_spherical_data', 'je should be a even number!!!', FATAL)

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

    allocate(wsf4m(0:num_fourier,2))
    allocate(wef4m(0:num_fourier,2))
    allocate(wlenf4m(0:num_fourier,2))


!--------------------------------------------------------------------------------   
    !global domain

    nsf4m(:) = 0
    do m = 0, num_fourier
       neadj = num_fourier + 1 - m
       nef4m(m) = neadj
    enddo 

    nlenf4m(:) = nef4m - nsf4m + 1

    nwaves = 0; noddwaves = 0; nevenwaves = 0

    do m = 0, num_fourier
        nwaves = nwaves + nlenf4m(m)
    enddo
    
    do m = 0, num_fourier
        do n = nsf4m(m), nef4m(m)
            if (mod(n,2)==0) then
                nevenwaves = nevenwaves + 1
            else
                noddwaves = noddwaves + 1
            endif 
        enddo
        wef4m(m,ev) = nevenwaves
        wef4m(m,od) = noddwaves 
    enddo

    if (noddwaves==nevenwaves) call mpp_error('init_spherical_data', 'global noddwaves==nevenwaves', FATAL)

    nevenwaves_g = nevenwaves
    noddwaves_g = noddwaves
    nwaves_oe_g = max(nevenwaves,noddwaves)
    
    wsf4m(0,:) = 1
    
    do m = 1, num_fourier
        wsf4m(m,:) = wef4m(m-1,:)+1
    enddo

    wlenf4m(:,:) = wef4m(:,:) - wsf4m(:,:) + 1

   nwaves = 0; noddwaves = 0; nevenwaves = 0
   !--------------------------------------------------------------------------------   
    !local

    allocate(ns4m(ms:me))
    allocate(ne4m(ms:me))
    allocate(nlen4m(ms:me))

    allocate(ws4m(ms:me,2))
    allocate(we4m(ms:me,2))
    allocate(wlen4m(ms:me,2))

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
        we4m(m,ev) = nevenwaves
        we4m(m,od) = noddwaves 
    enddo

    if (noddwaves/=nevenwaves) call mpp_error('init_spherical_data', 'noddwaves/=nevenwaves', FATAL)
    nwaves_oe = max(noddwaves,nevenwaves)
    
    nwaves_oe_out = nwaves_oe

    ws4m(ms,ev) = 1
    ws4m(ms,od) = 1
    
    do m = ms+1, me
        ws4m(m,ev) = we4m(m-1,ev)+1
        ws4m(m,od) = we4m(m-1,od)+1
    enddo

    wlen4m(:,ev) = we4m(:,ev) - ws4m(:,ev) + 1
    wlen4m(:,od) = we4m(:,od) - ws4m(:,od) + 1

    allocate(wdecomp(nwaves_oe,2))

    we = 0
    wo = 0
    wdecomp = 0
    do m = ms, me
        do n = 1, wlen4m(m,ev)
            we=we+1
            if(n<=wlenf4m(tshuffle(m),ev)) then
                wdecomp(we,ev) = wsf4m(tshuffle(m),ev) + n - 1
            endif
        enddo
        do n = 1, wlen4m(m,od)
            wo=wo+1
            if(n<=wlenf4m(tshuffle(m),od)) then
                wdecomp(wo,od) = wsf4m(tshuffle(m),od) + n - 1
            endif
        enddo
    enddo

    if (debug) then
        print *, 'ws4m=', ws4m 
        print *, 'wlen4m=', wlen4m 
        write(suffix,'(I4.4)') mpp_pe()
        print *, 'debug from fourier_spherical, pe= ', trim(suffix)
        print *, 'pe, noddwaves, nevenwaves =', mpp_pe(), noddwaves, nevenwaves
        print *, 'pe, ns4m(:)=', ns4m(:)
        print *, 'pe, ne4m(:)=', ne4m(:)

        !print *, 'pe, wdecomp(:,ev)=', wdecomp(:,ev)
        !print *, 'pe, wdecomp(:,od)=', wdecomp(:,od)
    endif

    deallocate(wlenf4m,wsf4m,wef4m)

end subroutine init_spherical_data

end module spherical_data_mod
