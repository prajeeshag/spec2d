module spherical_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE, mpp_init, mpp_pe
use mpp_mod, only : mpp_root_pe, mpp_sync, mpp_clock_id, mpp_clock_begin
use mpp_mod, only : mpp_clock_end

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_layout, &
                            mpp_get_global_domain, mpp_get_domain_extents

use fms_mod, only : fms_init, open_namelist_file, close_file

use fms_io_mod, only : write_data

use constants_mod, only : RADIUS, PI

use gauss_and_legendre_mod, only : compute_legendre, compute_gaussian 

use ocpack_mod, only : ocpack_typeP, ocpack_typeF, get_ocpackP, get_ocpackF, oc_nx, oc_ny, &
                oc_nlat, npack=>oc_npack, hem_type, get_hem

use spec_comm_mod, only : split_pelist

!-------------------------------------------------------------------------
!   provides operations on spectral spherical harmonics fields that do not 
!      require transforms
!
!   spectral fields are complex with horizontal dimensions
!        (0:num_fourier,0:num_spherical)
!       where the zonal wavenumber of (m,n) is M = m*fourier_inc
!             the "meridional" wavenumber is n
!             the "total", 2D, spherical wavenumber L = M+n
!
!-----------------------------------------------------------------------


implicit none

include 'mpif.h'

private

public :: init_spherical, get_wdecomp, get_latsF, get_latsP, get_wdecompa
public :: compute_lon_deriv_cos, compute_lat_deriv_cos
public :: compute_ucos_vcos, compute_vor_div, compute_vor
public :: compute_div, do_truncation, get_spherical_wave
public :: fourier_to_spherical, spherical_to_fourier

real, dimension(:,:), allocatable :: eigen_laplacian
real, dimension(:,:), allocatable :: epsilon
real, dimension(:,:), allocatable :: r_epsilon
real, dimension(:,:), allocatable :: coef_uvm
real, dimension(:,:), allocatable :: coef_uvc
real, dimension(:,:), allocatable :: coef_uvp
real, dimension(:,:), allocatable :: coef_alpm
real, dimension(:,:), allocatable :: coef_alpp
real, dimension(:,:), allocatable :: coef_dym
real, dimension(:,:), allocatable :: coef_dx
real, dimension(:,:), allocatable :: coef_dyp

integer, dimension(:,:), allocatable :: triangle_mask
integer, dimension(:,:), allocatable :: triangle_mask1
integer, dimension(:,:), allocatable :: fourier_wave
integer, dimension(:,:), allocatable :: spherical_wave
integer, dimension(:,:), allocatable :: nnp1

real, allocatable, dimension(:) :: sin_latF
real, allocatable, dimension(:) :: cos_latF
real, allocatable, dimension(:) :: cosm_latF
real, allocatable, dimension(:) :: cosm2_latF
real, allocatable, dimension(:) :: deg_latF
real, allocatable, dimension(:) :: wts_latF

real, allocatable, dimension(:,:) :: sin_latP
real, allocatable, dimension(:,:) :: cos_latP
real, allocatable, dimension(:,:) :: cosm_latP
real, allocatable, dimension(:,:) :: cosm2_latP
real, allocatable, dimension(:,:) :: deg_latP
real, allocatable, dimension(:,:) :: wts_latP

real, allocatable, dimension(:) :: sin_hem
real, allocatable, dimension(:) :: wts_hem

real, dimension(:,:,:), allocatable :: Pnm, Pnm_wts
real, dimension(:,:,:), allocatable :: Hnm, Hnm_wts

integer :: num_fourier = 0
integer :: num_spherical = 0
integer :: trunc = 0
integer :: truncadj = 0
integer :: ms, me, mlen
integer :: js_hem, je_hem, jlen_hem
integer :: jsf, jef, jlenf, nlat

integer, allocatable :: ns4m(:)   ! Starting of spherical for a particular fourier
integer, allocatable :: ne4m(:)   ! Ending of spherical for a particulat fourier 
                                  !(=num_spherical in case of Triangular
                                  !truncation)
integer, allocatable :: nlen4m(:) ! number of spherical for a particular fourier
                                  !(Constant in case of rhomboidal truncation)

integer :: nwaves = 0 !Total number of waves
integer :: nwaves_oe = 0 !Total number of odd-even waves [=nwaves/2])
integer :: noddwaves = 0, nevenwaves = 0
integer :: noddwaves_g = 0, nevenwaves_g = 0, nwaves_oe_g = 0
integer :: nwavesglobala

logical, allocatable :: iseven(:) ! oddeven flag (.true. = even); (.false. = odd)

integer, allocatable :: ws4m(:,:), we4m(:,:), wlen4m(:,:) !starting and ending index of 
                                                                  !even & odd waves for a particular m
integer, allocatable :: wdecomp(:,:)

integer, allocatable :: wdecompa(:,:)

integer, allocatable :: tshuffle(:)

type(ocpack_typeP), allocatable :: ocP(:,:)
type(ocpack_typeF), allocatable :: ocF(:)

type(hem_type), allocatable :: jh(:)

logical :: debug

integer :: clck_f2s, clck_s2f

integer, parameter, public :: ev=1, od=2

integer :: f_x_comm, layout(2)

logical :: initialized=.false., notfpe=.false.

interface init_spherical
    module procedure init_spherical1 
    module procedure init_spherical2
end interface init_spherical

contains

!------------------------------------------------------------------------------
subroutine init_spherical1(trunc_in, nwaves_oe_out, &
                           domain_fourier_in, tshuffle_in)
!------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: trunc_in
    type(domain2d) :: domain_fourier_in
    integer, optional :: tshuffle_in(0:)
    integer, intent(out) :: nwaves_oe_out

    integer :: m, w, n, unit, iostat, neadj
    integer :: nsf4m(0:trunc_in), nef4m(0:trunc_in), nlenf4m(0:trunc_in)
    integer :: nsf4ma(0:trunc_in), nef4ma(0:trunc_in), nlenf4ma(0:trunc_in)
    integer, allocatable :: wsf4m(:,:), wef4m(:,:), wlenf4m(:,:)
    integer, allocatable :: wsf4ma(:,:), wef4ma(:,:), wlenf4ma(:,:)
    integer :: we, wo, j, mee, npes
    integer, allocatable :: xextent(:), yextent(:), pes(:)
    character (len=8) :: suffix

    namelist/spherical_nml/debug

    call mpp_init()
    call fms_init()

    unit = open_namelist_file()
    read(unit,nml=spherical_nml,iostat=iostat)
    call close_file(unit)

    nlat = oc_nlat()

    if (mod(nlat,2)/=0) &
        call mpp_error('init_spherical', 'NLAT should be a even number', FATAL)

    allocate(ocP(npack(),oc_ny()))
    allocate(ocF(nlat))

    call get_ocpackP(ocP)
    call get_ocpackF(ocF)
    allocate(jh(nlat/2))
    call get_hem(jh)
   
    num_fourier = trunc_in
    num_spherical = trunc_in + 1
    trunc = trunc_in

    if (mod(num_spherical,2)==0) num_spherical = num_spherical + 1

    call mpp_get_compute_domain(domain_fourier_in,jsf,jef,ms,me)
    jlenf = jef - jsf + 1
    mlen = me-ms+1

    if (mod(jsf,2)==0) call mpp_error('init_spherical', 'jsf should be a odd number!!!', FATAL)
    if (mod(jlenf,2)/=0) call mpp_error('init_spherical', 'jlenf should be a even number!!!', FATAL)
    if (mod(jef,2)/=0) call mpp_error('init_spherical', 'jef should be a even number!!!', FATAL)

    js_hem = jsf/2 + 1
    je_hem = (jef-1)/2 + 1
    jlen_hem = jlenf/2


    allocate(tshuffle(ms:me))

    if(present(tshuffle_in)) then
        tshuffle(ms:me) = tshuffle_in(ms:me)
    else
        forall(m=ms:me) tshuffle(m) = m
    endif

    call mpp_get_layout(domain_fourier_in,layout)

    allocate(yextent(layout(1)))
    allocate(xextent(layout(2)))
    allocate(pes(layout(1)))

    call mpp_get_domain_extents(domain_fourier_in, yextent, xextent)
 
    mee = -1
    do m = 1, layout(2)
        mee  = mee + xextent(m)
        call split_pelist(mee==me,pes,npes,f_x_comm)
    end do

    if (debug) then
        call mpp_error(NOTE,"after comm split f_x_comm")
        if (mpp_pe()==mpp_root_pe()) then
            print *, "--> spherical x npes: ", npes
            print *, "--> spherical x pelist: ", pes
            print *, "--> spherical x comm id: ", f_x_comm
        endif
    end if
 

    allocate(wsf4m(0:num_fourier,2))
    allocate(wef4m(0:num_fourier,2))
    allocate(wlenf4m(0:num_fourier,2))

    allocate(wsf4ma(0:num_fourier,2))
    allocate(wef4ma(0:num_fourier,2))
    allocate(wlenf4ma(0:num_fourier,2))

    !global domain

    nsf4m(:) = 0
    do m = 0, num_fourier
       neadj = num_fourier + 1 - m
       nef4m(m) = neadj
    enddo 

    nlenf4m(:) = nef4m - nsf4m + 1

    nsf4ma(:) = 0
    do m = 0, num_fourier
       neadj = num_fourier - m
       nef4ma(m) = neadj
    enddo 

    nlenf4ma(:) = nef4ma - nsf4ma + 1

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

    if (noddwaves==nevenwaves) call mpp_error('init_spherical', 'global noddwaves==nevenwaves', FATAL)

    nevenwaves_g = nevenwaves
    noddwaves_g = noddwaves
    nwaves_oe_g = max(nevenwaves,noddwaves)
    
    wsf4m(0,:) = 1
    
    do m = 1, num_fourier
        wsf4m(m,:) = wef4m(m-1,:)+1
    enddo

    wlenf4m(:,:) = wef4m(:,:) - wsf4m(:,:) + 1

    nwaves = 0; noddwaves = 0; nevenwaves = 0

    do m = 0, num_fourier
        do n = nsf4ma(m), nef4ma(m)
           nwaves = nwaves + 1
        enddo
        if (mod(nef4ma(m),2)==0) then
            wef4ma(m,ev) = nwaves
            wef4ma(m,od) = nwaves - 1
        else
            wef4ma(m,ev) = nwaves - 1
            wef4ma(m,od) = nwaves
        endif
    enddo

    nwavesglobala = nwaves 

    wsf4ma(0,1) = 1
    wsf4ma(0,2) = wsf4ma(0,1) + 1
    
    do m = 1, num_fourier
        if (mod(nef4ma(m-1),2) == 0) then
            wsf4ma(m,1) = wef4ma(m-1,1)+1
        else
            wsf4ma(m,1) = wef4ma(m-1,2)+1
        endif
        wsf4ma(m,2) = wsf4ma(m,1) + 1
    enddo

    wlenf4ma(:,:) = (wef4ma(:,:) - wsf4ma(:,:))/2 + 1

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

    if (noddwaves/=nevenwaves) call mpp_error('init_spherical', 'noddwaves/=nevenwaves', FATAL)
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
    allocate(wdecompa(nwaves_oe,2))

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

    we = 0
    wo = 0
    wdecompa = 0
    do m = ms, me
        do n = 1, wlen4m(m,ev)
            we=we+1
            if(n<=wlenf4ma(tshuffle(m),ev)) then
                wdecompa(we,ev) = wsf4ma(tshuffle(m),ev) + (n - 1)*2
            endif
        enddo
        do n = 1, wlen4m(m,od)
            wo=wo+1
            if(n<=wlenf4ma(tshuffle(m),od)) then
                wdecompa(wo,od) = wsf4ma(tshuffle(m),od) + (n - 1)*2
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

    call spherical_init()

    clck_f2s = mpp_clock_id('fourier2spherical')
    clck_s2f = mpp_clock_id('spherical2fourier')

    initialized = .true.

end subroutine init_spherical1


subroutine init_spherical2()
    notfpe = .true.

    nlat = oc_nlat()

    if (mod(nlat,2)/=0) &
        call mpp_error('init_spherical', 'NLAT should be a even number', FATAL)

    allocate(ocP(npack(),oc_ny()))
    allocate(ocF(nlat))

    call get_ocpackP(ocP)
    call get_ocpackF(ocF)

    ms = 0; me = -1
    jsf = 0; jef = -1
    jlenf = jef - jsf + 1

    js_hem = 0
    je_hem = -1
    jlen_hem = 0

    call define_gaussian()

    initialized = .true.
end subroutine init_spherical2

!--------------------------------------------------------------------------------   
subroutine get_latsP(sinlat,coslat,cosmlat,cosm2lat,deglat,wtslat)
!--------------------------------------------------------------------------------
    real, intent(out), optional, dimension(:,:) :: sinlat, coslat, cosmlat, &
                                                 cosm2lat, deglat, wtslat

    if(.not.initialized) &
        call mpp_error('spherical_mod', 'module not initialized', FATAL)

    if(present(sinlat)) then
        if (sum(shape(sinlat)-[oc_ny(),oc_nx()])==0) then
            sinlat = sin_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    if(present(coslat)) then
        if (sum(shape(coslat)-[oc_ny(),oc_nx()])==0) then
            coslat = cos_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    if(present(cosmlat)) then
        if (sum(shape(cosmlat)-[oc_ny(),oc_nx()])==0) then
            cosmlat = cosm_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    if(present(cosm2lat)) then
        if (sum(shape(cosm2lat)-[oc_ny(),oc_nx()])==0) then
            cosm2lat = cosm2_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    if(present(deglat)) then
        if (sum(shape(deglat)-[oc_ny(),oc_nx()])==0) then
            deglat = deg_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    if(present(wtslat)) then
        if (sum(shape(wtslat)-[oc_ny(),oc_nx()])==0) then
            wtslat = wts_latP
        else
            call mpp_error(FATAL,'get_latsP: array size mismatch!')
        endif
    endif

    return
end subroutine get_latsP


!--------------------------------------------------------------------------------   
subroutine get_latsF(sinlat,coslat,cosmlat,cosm2lat,deglat,wtslat)
!--------------------------------------------------------------------------------
    real, intent(out), optional, dimension(:) :: sinlat, coslat, cosmlat, &
                                                 cosm2lat, deglat, wtslat

    if(.not.initialized) &
        call mpp_error('spherical_mod', 'module not initialized', FATAL)

    if(present(sinlat)) then
        if (size(sinlat,1)==nlat) then
            sinlat = sin_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    if(present(coslat)) then
        if (size(coslat,1)==nlat) then
            coslat = cos_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    if(present(cosmlat)) then
        if (size(cosmlat,1)==nlat) then
            cosmlat = cosm_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    if(present(cosm2lat)) then
        if (size(cosm2lat,1)==nlat) then
            cosm2lat = cosm2_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    if(present(deglat)) then
        if (size(deglat,1)==nlat) then
            deglat = deg_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    if(present(wtslat)) then
        if (size(wtslat,1)==nlat) then
            wtslat = wts_latF
        else
            call mpp_error(FATAL,'get_latsF: array size mismatch!')
        endif
    endif

    return
end subroutine get_latsF


!--------------------------------------------------------------------------------   
subroutine get_wdecompa(wdom,nwavesglobal)
!--------------------------------------------------------------------------------
    integer, intent(out) :: wdom(:,:)
    integer, intent(out), optional :: nwavesglobal
    
    if (notfpe) return

    if(.not.initialized) &
        call mpp_error('get_wdecomp', 'module not initialized', FATAL)

    wdom(:,:) = wdecompa(:,:)

    if (present(nwavesglobal)) nwavesglobal = nwavesglobala

end subroutine get_wdecompa

!--------------------------------------------------------------------------------   
subroutine get_wdecomp(wdom,neven_global,nodd_global)
!--------------------------------------------------------------------------------
    integer, intent(out) :: wdom(:,:)
    integer, intent(out), optional :: neven_global, nodd_global
    
    if (notfpe) return

    if(.not.initialized) &
        call mpp_error('get_wdecomp', 'module not initialized', FATAL)

    wdom(:,:) = wdecomp(:,:)

    if (present(neven_global)) neven_global = nevenwaves_g
    if (present(nodd_global)) nodd_global = noddwaves_g

end subroutine get_wdecomp


!--------------------------------------------------------------------------------   
subroutine get_spherical_wave(spherical_wave_out,nnp1_out)
!--------------------------------------------------------------------------------   
    integer, intent(out), optional :: spherical_wave_out(:,:)
    integer, intent(out), optional :: nnp1_out(:,:)

    if (notfpe) return

    if(.not.initialized) &
        call mpp_error('get_spherical_wave', 'module not initialized', FATAL)

    if(present(spherical_wave_out)) spherical_wave_out = spherical_wave
    if(present(nnp1_out)) nnp1_out = nnp1

end subroutine get_spherical_wave


!--------------------------------------------------------------------------
subroutine spherical_init()
!--------------------------------------------------------------------------------   

    real, dimension(0:num_fourier,0:num_spherical) :: leigen_laplacian, &
                                                      lepsilon, &
                                                      r_lepsilon, &
                                                      lfourier_wave, &
                                                      lspherical_wave, &
                                                      lcoef_uvm, &
                                                      lcoef_uvc, &
                                                      lcoef_uvp, &
                                                      lcoef_alpm, &
                                                      lcoef_alpp, &
                                                      lcoef_dym, &
                                                      lcoef_dx, &
                                                      lcoef_dyp

    integer, dimension(0:num_fourier,0:num_spherical) :: ltriangle_mask, &
                                                         ltriangle_mask1, &
                                                         lnnp1

    integer :: m, n, we, wo, ma, np1, nm1, nums, numf

    if(initialized) return

    nums = num_spherical; numf = num_fourier

    ltriangle_mask = 1
    ltriangle_mask1 = 1
    do n=0,num_spherical
      do m=0,num_fourier
        lfourier_wave(m,n)   = m
        lspherical_wave(m,n) = lfourier_wave(m,n) + n
        lnnp1(m,n) = lspherical_wave(m,n)*(lspherical_wave(m,n)+1) 
        if(lspherical_wave(m,n)>trunc) ltriangle_mask(m,n) = 0
        if(lspherical_wave(m,n)>trunc+1) ltriangle_mask1(m,n) = 0
      end do
    end do
    
    lepsilon= (lspherical_wave**2 - lfourier_wave**2)/(4.0*lspherical_wave**2 - 1.0)
    lepsilon = sqrt(lepsilon)
    
    leigen_laplacian = lspherical_wave*(lspherical_wave + 1.0)/(radius*radius)
    
    where (lspherical_wave > 0) 
      lcoef_uvc = -radius*lfourier_wave/(lspherical_wave*(lspherical_wave + 1.0))
    else where 
      lcoef_uvc = 0.0
    end where
   
    lcoef_uvm = 0. 
    lcoef_uvm(:,0:nums-1) = -radius*lepsilon(:,1:nums) &
                            / lspherical_wave(:,1:nums)
    
    lcoef_uvp = 0. 
    lcoef_uvp(:,1:nums) = -radius*lepsilon(:,1:nums) &
                          / (lspherical_wave(:,0:nums-1) + 1.0)
    
    lcoef_alpm(:,:) = 0.
    lcoef_alpm(:,0:nums-1) = (lspherical_wave(:,1:nums) + 1.0) &
                            * lepsilon(:,1:nums) / radius

    lcoef_alpp(:,:) = 0.
    lcoef_alpp(:,1:nums) = lspherical_wave(:,0:nums-1)*lepsilon(:,1:nums)/radius
    
    lcoef_dym = 0.
    lcoef_dym(:,0:num_spherical-1) = (lspherical_wave(:,1:num_spherical) - 1.0) &
                                   * lepsilon(:,1:num_spherical) / radius

    lcoef_dx = lfourier_wave/radius

    lcoef_dyp(:,0:num_spherical-1) = 0.
    lcoef_dyp(:,1:num_spherical) =                                           &
      (lspherical_wave(:,0:num_spherical-1) + 2.0)*lepsilon(:,1:num_spherical)/radius

    allocate(eigen_laplacian(nwaves_oe,2), &
             epsilon(nwaves_oe,2), &
             r_epsilon(nwaves_oe,2), &
             fourier_wave(nwaves_oe,2), &
             spherical_wave(nwaves_oe,2), &
             coef_uvm(nwaves_oe,2), &
             coef_uvc(nwaves_oe,2), &
             coef_uvp(nwaves_oe,2), &
             coef_alpm(nwaves_oe,2), &
             coef_alpp(nwaves_oe,2), &
             coef_dym(nwaves_oe,2), &
             coef_dx(nwaves_oe,2), &
             coef_dyp(nwaves_oe,2), &
             triangle_mask(nwaves_oe,2), &
             triangle_mask1(nwaves_oe,2), &
             nnp1(nwaves_oe,2))


    wo = 0
    we = 0

    do m = ms, me
        ma = tshuffle(m)
        do n = ns4m(m), ne4m(m)
            np1 = n + 1
            nm1 = n - 1
            if (mod(n,2)==0) then
                we = we + 1
                eigen_laplacian(we,1) = leigen_laplacian(ma,n)
                epsilon(we,1)         = lepsilon(ma,n)
                r_epsilon(we,1)       = r_lepsilon(ma,n)
                fourier_wave(we,1)    = lfourier_wave(ma,n)
                spherical_wave(we,1)  = lspherical_wave(ma,n)
                coef_uvm(we,1)        = lcoef_uvm(ma,n)
                coef_uvc(we,1)        = lcoef_uvc(ma,n)
                coef_uvp(we,1)        = lcoef_uvp(ma,n)
                coef_alpm(we,1)       = lcoef_alpm(ma,n)
                coef_alpp(we,1)       = lcoef_alpp(ma,n)
                coef_dym(we,1)        = lcoef_dym(ma,n)
                coef_dx(we,1)         = lcoef_dx(ma,n)
                coef_dyp(we,1)        = lcoef_dyp(ma,n)
                triangle_mask(we,1)   = ltriangle_mask(ma,n)
                triangle_mask1(we,1)  = ltriangle_mask1(ma,n)
                nnp1(we,1)            = lnnp1(ma,n)
            else
                wo = wo + 1
                eigen_laplacian(wo,2) = leigen_laplacian(ma,n)
                epsilon(wo,2)         = lepsilon(ma,n)
                r_epsilon(wo,2)       = r_lepsilon(ma,n)
                fourier_wave(wo,2)    = lfourier_wave(ma,n)
                spherical_wave(wo,2)  = lspherical_wave(ma,n)
                coef_uvm(wo,2)        = lcoef_uvm(ma,n)
                coef_uvc(wo,2)        = lcoef_uvc(ma,n)
                coef_uvp(wo,2)        = lcoef_uvp(ma,n)
                coef_alpm(wo,2)       = lcoef_alpm(ma,n)
                coef_alpp(wo,2)       = lcoef_alpp(ma,n)
                coef_dym(wo,2)        = lcoef_dym(ma,n)
                coef_dx(wo,2)         = lcoef_dx(ma,n)
                coef_dyp(wo,2)        = lcoef_dyp(ma,n)
                triangle_mask(wo,2)   = ltriangle_mask(ma,n)
                triangle_mask1(wo,2)  = ltriangle_mask1(ma,n)
                nnp1(wo,2)            = lnnp1(ma,n)
            endif
        enddo
    enddo

    where (triangle_mask==0)
        eigen_laplacian = 0.
        coef_uvm = 0.
        coef_uvc = 0.
        coef_uvp = 0.
        coef_alpm = 0.
        coef_dym = 0.
        coef_dx = 0.
        coef_dyp = 0.
    end where

    call define_gaussian

    call define_legendre(lepsilon,lspherical_wave)

    initialized = .true.

return
end subroutine spherical_init


!------------------------------------------------------------------------------
subroutine define_gaussian
!------------------------------------------------------------------------------

    integer :: j, is, ie, i
    real :: sin_hem1(nlat), wts_hem1(nlat)

    allocate (sin_latF(nlat))
    allocate (cos_latF(nlat))
    allocate (cosm_latF(nlat))
    allocate (cosm2_latF(nlat))
    allocate (wts_latF(nlat))
    allocate (deg_latF(nlat))

    allocate (sin_latP(oc_ny(),oc_nx()))
    allocate (cos_latP(oc_ny(),oc_nx()))
    allocate (cosm_latP(oc_ny(),oc_nx()))
    allocate (cosm2_latP(oc_ny(),oc_nx()))
    allocate (wts_latP(oc_ny(),oc_nx()))
    allocate (deg_latP(oc_ny(),oc_nx()))

    allocate (sin_hem(nlat/2))
    allocate (wts_hem(nlat/2))

    call compute_gaussian(sin_hem, wts_hem, nlat/2)

    sin_hem1(1:nlat/2)         = -sin_hem !Southern hemisphere
    sin_hem1(nlat:nlat/2+1:-1) = sin_hem !Northern hemisphere

    wts_hem1(1:nlat/2)         = wts_hem
    wts_hem1(nlat:nlat/2+1:-1) = wts_hem
    
    do j = 1, nlat
        sin_latF(j) = sin_hem1(ocF(j)%g) 
        wts_latF(j) = wts_hem1(ocF(j)%g)
    end do

    cos_latF = sqrt(1-sin_latF*sin_latF)
    cosm_latF = 1./cos_latF
    cosm2_latF = 1./(cos_latF*cos_latF)
    deg_latF = asin(sin_latF)*180.0/pi

    do i = 1, npack()
        do j = 1, oc_ny()
            is = ocP(i,j)%is
            ie = ocP(i,j)%ie
            sin_latP(j,is:ie) = sin_hem1(ocP(i,j)%g) 
            wts_latP(j,is:ie) = wts_hem1(ocP(i,j)%g)
        end do
    end do

    cos_latP = sqrt(1-sin_latP*sin_latP)
    cosm_latP = 1./cos_latP
    cosm2_latP = 1./(cos_latP*cos_latP)
    deg_latP = asin(sin_latP)*180.0/pi

    return
end subroutine define_gaussian


!--------------------------------------------------------------------------------   
subroutine define_legendre(lepsilon,lspherical_wave)
!--------------------------------------------------------------------------------   
    real, dimension(0:num_fourier,0:num_spherical), intent(in) :: lepsilon, lspherical_wave
    integer :: j, m, w, wo, we, mshuff, n, jg
    real, dimension(0:num_fourier,0:num_spherical,nlat/2) :: Pnm_global
    real, dimension(0:num_fourier,0:num_spherical,nlat/2) :: Hnm_global
    real :: wgt
    character(len=8) :: suffix

    allocate(Pnm(jlen_hem,nwaves_oe,2))
    allocate(Pnm_wts(jlen_hem,nwaves_oe,2))
    allocate(Hnm(jlen_hem,nwaves_oe,2))
    allocate(Hnm_wts(jlen_hem,nwaves_oe,2))

    call compute_legendre(Pnm_global, num_fourier, 1, num_spherical, sin_hem, nlat/2)

    do m = 0, num_fourier
        do n = 0, num_spherical-1
            Hnm_global(m,n,:) = -lspherical_wave(m,n) &
                                          * lepsilon(m,n+1) &
                                          * Pnm_global(m,n+1,:)
        enddo 
        do n = 1, num_spherical
            Hnm_global(m,n,:) = Hnm_global(m,n,:) &
                                          + (lspherical_wave(m,n)+1) &
                                          * lepsilon(m,n) * Pnm_global(m,n-1,:)
        enddo
    enddo
  
    wgt = 1./RADIUS
    Hnm_global(:,:,:) = Hnm_global(:,:,:)*wgt

    w = 0
    wo = 0
    we = 0
    do m = ms, me
        mshuff = tshuffle(m)
        do n = ns4m(m), ne4m(m)
            w = w + 1
            if (iseven(w)) then
                we = we + 1
                do j = js_hem, je_hem
                    jg = jh(j)%g
                    Pnm(j-js_hem+1,we,1) = Pnm_global(mshuff,n,jg) 
                    Hnm(j-js_hem+1,we,1) = Hnm_global(mshuff,n,jg) 
                end do
            else
                wo = wo + 1
                do j = js_hem, je_hem
                    jg = jh(j)%g
                    Pnm(j-js_hem+1,wo,2) = Pnm_global(mshuff,n,jg) 
                    Hnm(j-js_hem+1,wo,2) = Hnm_global(mshuff,n,jg) 
                end do
            endif
        enddo
    enddo

    do j = js_hem, je_hem
        jg = jh(j)%g
        Pnm_wts(j-js_hem+1,:,:) = Pnm(j-js_hem+1,:,:)*wts_hem(jg)
        Hnm_wts(j-js_hem+1,:,:) = Hnm(j-js_hem+1,:,:)*wts_hem(jg)

        where(triangle_mask==0) Hnm(j-js_hem+1,:,:) = 0.
        
         where(triangle_mask1==0) 
             Hnm_wts(j-js_hem+1,:,:) = 0.
             Pnm_wts(j-js_hem+1,:,:) = 0.
             Pnm(j-js_hem+1,:,:) = 0.
         end where
    enddo

    return
end subroutine define_legendre


!---------------------------------------------------------------------------
subroutine compute_lon_deriv_cos(spherical, deriv_lon)
!---------------------------------------------------------------------------

    complex, dimension(:,:,:), intent(in) :: spherical
    complex, dimension(:,:,:), intent(out) :: deriv_lon
    
    integer :: k, n, i
    
    if (notfpe) return

    if(.not. initialized ) then
      call mpp_error('compute_lon_deriv','module spherical not initialized', FATAL)
    end if

    do i = 1, 2    
        do n = 1, size(spherical,2)
           deriv_lon(:,n,i) = coef_dx(n,i)*cmplx(-aimag(spherical(:,n,i)),real(spherical(:,n,i)))
        enddo
    enddo

    return
end subroutine compute_lon_deriv_cos


!--------------------------------------------------------------------------
subroutine compute_lat_deriv_cos(spherical,deriv_lat)
!--------------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in) :: spherical
    complex,dimension(:,:,:), intent(out) :: deriv_lat
    
    integer :: k, n, nw
    
    if (notfpe) return

    if(.not. initialized ) then
      call mpp_error('compute_lat_deriv','module spherical not initialized', FATAL)
    end if

    nw = size(spherical,2)
 
    deriv_lat = cmplx(0.,0.)

    do k = 1, size(spherical,1)

       deriv_lat(k,2:nw,ev) = &  
            - spherical(k,1:nw-1,od)*coef_dym(1:nw-1,od)

       deriv_lat(k,1:nw,ev) =  deriv_lat(k,1:nw,ev)   &
            + spherical(k,1:nw,od)*coef_dyp(1:nw,od)

       deriv_lat(k,1:nw,od) = &  
            - spherical(k,1:nw,ev)*coef_dym(1:nw,ev)

       deriv_lat(k,1:nw-1,od) =  deriv_lat(k,1:nw-1,od)   &
            + spherical(k,2:nw,ev)*coef_dyp(2:nw,ev)

    end do

    return
end subroutine compute_lat_deriv_cos


!--------------------------------------------------------------------------------   
subroutine do_truncation(spherical,full)
!--------------------------------------------------------------------------------   
    complex,dimension(:,:,:), intent(inout) :: spherical
    logical, optional :: full

    integer :: k
    logical :: full1

    if (notfpe) return

    full1 = .true.
    if(present(full)) full1=full

    if (full1) then
        do k = 1, size(spherical,1)
            spherical(k,:,:) = spherical(k,:,:) * triangle_mask(:,:)
        enddo
    else
        do k = 1, size(spherical,1)
            spherical(k,:,:) = spherical(k,:,:) * triangle_mask1(:,:)
        enddo
    endif

    return

end subroutine do_truncation


!----------------------------------------------------------------------
subroutine compute_ucos_vcos(vorticity , divergence, u_cos, v_cos, do_trunc)
!----------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in)  :: vorticity
    complex,dimension(:,:,:), intent(in)  :: divergence
    complex,dimension(:,:,:), intent(out) :: u_cos
    complex,dimension(:,:,:), intent(out) :: v_cos
    logical, intent(in), optional :: do_trunc

    logical :: do_trunc1
    integer :: k, nw

    if (notfpe) return

    if(.not. initialized ) then
      call mpp_error('compute_ucos_vcos','module spherical not initialized', FATAL)
    end if

    do_trunc1 = .true.
   
    if(present(do_trunc)) do_trunc1=do_trunc 

    nw = size(vorticity,2)

    do k=1,size(vorticity,1)

       u_cos(k,:,ev) = coef_uvc(:,ev)*                                     &
            cmplx(-aimag(divergence(k,:,ev)),real(divergence(k,:,ev)))

       u_cos(k,:,od) = coef_uvc(:,od)*                                     &
            cmplx(-aimag(divergence(k,:,od)),real(divergence(k,:,od)))

       u_cos(k,2:nw,ev) = u_cos(k,2:nw,ev) +         &
            coef_uvm(1:nw-1,od)*vorticity(k,1:nw-1,od)

       u_cos(k,1:nw,od) = u_cos(k,1:nw,od) +         &
            coef_uvm(1:nw,ev)*vorticity(k,1:nw,ev)

       u_cos(k,1:nw,ev) = u_cos(k,1:nw,ev) -     &
            coef_uvp(1:nw,od)*vorticity(k,1:nw,od)

       u_cos(k,1:nw-1,od) = u_cos(k,1:nw-1,od) -     &
            coef_uvp(2:nw,ev)*vorticity(k,2:nw,ev)

       v_cos(k,:,ev) = coef_uvc(:,ev)*                                     &
            cmplx(-aimag(vorticity(k,:,ev)),real(vorticity(k,:,ev)))

       v_cos(k,:,od) = coef_uvc(:,od)*                                     &
            cmplx(-aimag(vorticity(k,:,od)),real(vorticity(k,:,od)))

       v_cos(k,2:nw,ev) = v_cos(k,2:nw,ev) -         &
            coef_uvm(1:nw-1,od)*divergence(k,1:nw-1,od)

       v_cos(k,1:nw,od) = v_cos(k,1:nw,od) -         &
            coef_uvm(1:nw,ev)*divergence(k,1:nw,ev)

       v_cos(k,1:nw,ev) = v_cos(k,1:nw,ev) +     &
            coef_uvp(1:nw,od)*divergence(k,1:nw,od)

       v_cos(k,1:nw-1,od) = v_cos(k,1:nw-1,od) +     &
            coef_uvp(2:nw,ev)*divergence(k,2:nw,ev)

    end do

    call do_truncation(v_cos,do_trunc1)
    call do_truncation(u_cos,do_trunc1)
            
    return
end subroutine compute_ucos_vcos

!--------------------------------------------------------------------------------
subroutine compute_alpha_operator(spherical_a, spherical_b, rsign, alpha, do_trunc)
!--------------------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in)  :: spherical_a
    complex,dimension(:,:,:), intent(in)  :: spherical_b
    complex,dimension(:,:,:), intent(out) :: alpha
    integer, intent(in) :: rsign
    logical, intent(in), optional :: do_trunc

    logical :: do_trunc1
    integer :: k, nw

    if (notfpe) return

    if(.not. initialized ) then
      call mpp_error('compute_vor or div','module spherical not initialized', FATAL)
    end if

    do_trunc1 = .true.
    if(present(do_trunc)) do_trunc1=do_trunc
    
    alpha = cmplx(0.,0.)
   
    nw = size(spherical_a,2)

    do k = 1, size(spherical_a,1)
       alpha(k,:,ev) = coef_dx(:,ev)*    &
            cmplx(-aimag(spherical_a(k,:,ev)),real(spherical_a(k,:,ev)))

       alpha(k,:,od) = coef_dx(:,od)*    &
            cmplx(-aimag(spherical_a(k,:,od)),real(spherical_a(k,:,od)))

       alpha(k,2:nw,ev) = alpha(k,2:nw,ev) -  &
            rsign*coef_alpm(1:nw-1,od)  &
            *spherical_b(k,1:nw-1,od)

       alpha(k,1:nw,od) = alpha(k,1:nw,od) -  &
            rsign*coef_alpm(1:nw,ev)  &
            *spherical_b(k,1:nw,ev)

       alpha(k,1:nw,ev) = alpha(k,1:nw,ev) +  &
            rsign*coef_alpp(1:nw,od)*spherical_b(k,1:nw,od)

       alpha(k,1:nw-1,od) = alpha(k,1:nw-1,od) +  &
            rsign*coef_alpp(2:nw,ev)*spherical_b(k,2:nw,ev)

    end do

    call do_truncation(alpha,do_trunc1)

    return
end subroutine compute_alpha_operator

!-------------------------------------------------------------------------
subroutine compute_vor_div(u_cos, v_cos, vorticity, divergence, do_trunc)
!-------------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in)  :: u_cos
    complex,dimension(:,:,:), intent(in)  :: v_cos
    complex,dimension(:,:,:), intent(out) :: vorticity
    complex,dimension(:,:,:), intent(out) :: divergence
    logical, intent(in), optional :: do_trunc

    if (notfpe) return

    call compute_alpha_operator(v_cos, u_cos, -1, vorticity, do_trunc)
    call compute_alpha_operator(u_cos, v_cos, +1, divergence, do_trunc)

    return
end subroutine compute_vor_div

!-------------------------------------------------------------------------
subroutine compute_vor(u_cos, v_cos, vorticity, do_trunc)
!-------------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in)  :: u_cos
    complex,dimension(:,:,:), intent(in)  :: v_cos
    complex,dimension(:,:,:), intent(out) :: vorticity
    logical, intent(in), optional :: do_trunc

    if (notfpe) return

    call compute_alpha_operator(v_cos, u_cos, -1, vorticity, do_trunc)

    return
end subroutine compute_vor

!-------------------------------------------------------------------------
subroutine compute_div(u_cos, v_cos, divergence, do_trunc)
!-------------------------------------------------------------------------

    complex,dimension(:,:,:), intent(in)  :: u_cos
    complex,dimension(:,:,:), intent(in)  :: v_cos
    complex,dimension(:,:,:), intent(out) :: divergence
    logical, intent(in), optional :: do_trunc

    if (notfpe) return

    call compute_alpha_operator(u_cos, v_cos, +1, divergence, do_trunc)

    return
end subroutine compute_div

!--------------------------------------------------------------------------------   
subroutine spherical_to_fourier(waves,fourier,useHnm)
!--------------------------------------------------------------------------------   
    complex, intent(out) :: fourier(:,jsf:,ms:) ! lat, lev, fourier
    complex,dimension(:,:,:), intent(in) :: waves
    logical, intent(in), optional :: useHnm

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    integer :: ks, ke, ews, ewe, ows, owe, m, nj, j, js, jn
    logical :: deriv

    if (notfpe) return

    call mpp_clock_begin(clck_s2f)

    deriv = .false.
    
    if (present(useHnm)) deriv = useHnm

    if (.not.initialized) call mpp_error('spherical_to_fourier', 'module not initialized', FATAL)

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
 
        do j = js_hem, je_hem 
            js = jh(j)%s
            jn = jh(j)%n
            fourier(ks:ke,jn,ms:me) = even(ks:ke,j,ms:me) + odd(ks:ke,j,ms:me)
            fourier(ks:ke,js,ms:me) = odd(ks:ke,j,ms:me)  - even(ks:ke,j,ms:me)
        end do

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

        do j = js_hem, je_hem 
            js = jh(j)%s
            jn = jh(j)%n
            fourier(ks:ke,jn,ms:me) = even(ks:ke,j,ms:me) + odd(ks:ke,j,ms:me)
            fourier(ks:ke,js,ms:me) = even(ks:ke,j,ms:me) - odd(ks:ke,j,ms:me)
        end do

    endif
   
    call mpp_clock_end(clck_s2f)

    return 
end subroutine spherical_to_fourier


!--------------------------------------------------------------------------------   
subroutine fourier_to_spherical(fourier, waves, useHnm, do_trunc)
!--------------------------------------------------------------------------------   
    complex, intent(in) :: fourier(:,jsf:,ms:) ! lat, lev, fourier
    complex,dimension(:,:,:), intent(out) :: waves
    logical, optional :: useHnm, do_trunc

    complex :: odd(size(fourier,1),js_hem:je_hem,ms:me)
    complex :: even(size(fourier,1),js_hem:je_hem,ms:me)

    logical :: useHnm1, do_trunc1
    integer :: ks, ke, ews, ewe, ows, owe, m, k, nj, j, js, jn
    integer :: rqst(mlen*2), nrq, nbuff, klen, stts(MPI_STATUS_SIZE,mlen*2), ierr
    complex, dimension(size(waves,1),size(waves,2),size(waves,3)) :: buff

    if (notfpe) return

    call mpp_clock_begin(clck_f2s)

    useHnm1 = .false.
    if(present(useHnm)) useHnm1=useHnm

    do_trunc1 = .true.
    if(present(do_trunc)) do_trunc1=do_trunc

    if (.not.initialized) call mpp_error('fourier_to_spherical', 'module not initialized', FATAL)

    ks = 1; ke = size(fourier,1); klen = ke-ks+1
   
    nrq = 0

    if (useHnm1) then
        do j = js_hem, je_hem 
            js = jh(j)%s
            jn = jh(j)%n
            odd(ks:ke,j,ms:me)  = fourier(ks:ke,jn,ms:me) + fourier(ks:ke,js,ms:me) ! 
            even(ks:ke,j,ms:me) = fourier(ks:ke,jn,ms:me) - fourier(ks:ke,js,ms:me) !
        end do

        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                           Hnm_wts(1:jlen_hem,ews:ewe,ev), &
                           buff(ks:ke,ews:ewe,ev),'N')
#ifdef MPI3
            if (layout(1)>1) then
                nbuff = klen * wlen4m(m,ev)
                nrq = nrq + 1
                call MPI_Iallreduce(buff(ks,ews,ev), waves(ks,ews,ev), nbuff, &
                            MPI_DOUBLE_COMPLEX, MPI_SUM, f_x_comm, rqst(nrq), ierr)
            endif
#endif
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                           Hnm_wts(1:jlen_hem,ows:owe,od), &
                           buff(ks:ke,ows:owe,od),'N')
#ifdef MPI3
            if (layout(1)>1) then
                nbuff = klen * wlen4m(m,od)
                nrq = nrq + 1
                call MPI_Iallreduce(buff(ks,ows,od), waves(ks,ows,od), nbuff, MPI_DOUBLE_COMPLEX, &
                                    MPI_SUM, f_x_comm, rqst(nrq), ierr)
            end if
#endif
        enddo 

    else

        do j = js_hem, je_hem 
            js = jh(j)%s
            jn = jh(j)%n
            odd(ks:ke,j,ms:me)  = fourier(ks:ke,jn,ms:me) - fourier(ks:ke,js,ms:me) ! 
            even(ks:ke,j,ms:me) = fourier(ks:ke,jn,ms:me) + fourier(ks:ke,js,ms:me) ! 
        end do

        do m = ms, me
            if(wlen4m(m,ev)<1) cycle
            ews = ws4m(m,ev); ewe = we4m(m,ev)
            call do_matmul(even(ks:ke,js_hem:je_hem,m), &
                           Pnm_wts(1:jlen_hem,ews:ewe,ev), &
                           buff(ks:ke,ews:ewe,ev),'N')
#ifdef MPI3
            if (layout(1)>1) then
                nbuff = klen * wlen4m(m,ev)
                nrq = nrq + 1
                call MPI_Iallreduce(buff(ks,ews,ev), waves(ks,ews,ev), nbuff, MPI_DOUBLE_COMPLEX, &
                                    MPI_SUM, f_x_comm, rqst(nrq), ierr)
            endif
#endif
        enddo

        do m = ms, me
            if(wlen4m(m,od)<1) cycle
            ows = ws4m(m,od); owe = we4m(m,od)
            call do_matmul(odd(ks:ke,js_hem:je_hem,m), &
                           Pnm_wts(1:jlen_hem,ows:owe,od), &
                           buff(ks:ke,ows:owe,od),'N')
#ifdef MPI3
            if (layout(1)>1) then
                nbuff = klen * wlen4m(m,od)
                nrq = nrq + 1
                call MPI_Iallreduce(buff(ks,ows,od), waves(ks,ows,od), nbuff, MPI_DOUBLE_COMPLEX, &
                                    MPI_SUM, f_x_comm, rqst(nrq), ierr)
            endif
#endif
        enddo 
    endif

    if (layout(1)>1) then
#ifdef MPI3
       call MPI_WAITALL(nrq, rqst(1), stts(1,1), ierr)
#else
       nbuff = size(buff)
       call MPI_allreduce(buff, waves, nbuff, MPI_DOUBLE_COMPLEX, MPI_SUM, f_x_comm, ierr)
       if (ierr/=0) call mpp_error(FATAL,'fourier_to_spherical: error in MPI_allreduce')
#endif
    else 
        waves = buff
    endif

    call do_truncation(waves,do_trunc1)

    call mpp_clock_end(clck_f2s)
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

end module spherical_mod

