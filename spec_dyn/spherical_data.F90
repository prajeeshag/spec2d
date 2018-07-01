module spherical_data_mod

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_layout

use constants_mod, only : pi

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

implicit none
private

public :: init_spherical_data

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

integer, public :: nwaves !Total number of spectral waves (local)
integer, public :: nwaves_oe !Total number of odd-even waves [=nwaves/2])
integer, public :: noddwaves, nevenwaves

logical, public, allocatable :: iseven(:) ! oddeven flag (.true. = even); (.false. = odd)

integer, public, allocatable :: ws4m(:), we4m(:), wlen4m(:) !starting and ending index of 
                                                            !waves for a particular m
integer, public, allocatable :: ews4m(:), ewe4m(:), ewlen4m(:) !starting and ending index of 
                                                               !even waves for a particular m
integer, public, allocatable :: ows4m(:), owe4m(:), owlen4m(:) !starting and ending index of 
                                                               !odd waves for a particular m
integer, public, allocatable :: tshuffle(:)

logical :: debug

integer, parameter, public :: ev=1, od=2

contains

!------------------------------------------------------------------------------
subroutine init_spherical_data(trunc_in, nlat_in, &
                nwaves_oe_out, domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    integer, intent(in) :: trunc_in, nlat_in
    type(domain2d), optional :: domain_fourier_in
    integer, optional :: tshuffle_in(0:)
    integer, intent(out) :: nwaves_oe_out

    integer :: m, w, n, unit, iostat, neadj
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

    !if (noddwaves/=nevenwaves) call mpp_error('init_spherical_data', 'noddwaves/=nevenwaves', FATAL)
    nwaves_oe = max(noddwaves,nevenwaves)
    
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
        write(suffix,'(I4.4)') mpp_pe()
        print *, 'debug from fourier_spherical, pe= ', trim(suffix)
        print *, 'pe, noddwaves, nevenwaves =', mpp_pe(), noddwaves, nevenwaves
        print *, 'pe, ns4m(:)=', ns4m(:)
        print *, 'pe, ne4m(:)=', ne4m(:)
    endif

end subroutine init_spherical_data

end module spherical_data_mod
