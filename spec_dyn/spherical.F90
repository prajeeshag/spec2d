module spherical_mod

use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_error, FATAL, WARNING, NOTE

use fms_io_mod, only: write_data

use spherical_data_mod, only: nwaves_oe, ms, me, nlat
use spherical_data_mod, only: ns4m, ne4m, num_spherical, num_fourier, trunc
use spherical_data_mod, only: tshuffle, js, je, js_hem, je_hem, jlen_hem
use spherical_data_mod, only: iseven, ev, od, get_wdecomp

use constants_mod, only : RADIUS, PI

use gauss_and_legendre_mod, only : compute_legendre, compute_gaussian 


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
private

public :: compute_lon_deriv_cos, compute_lat_deriv_cos, get_wdecomp
public :: nwaves_oe, num_spherical, num_fourier, trunc, ev, od
public :: spherical_init, compute_ucos_vcos, compute_vor_div, compute_vor
public :: compute_div, triangle_mask, do_truncation, nnp1, get_spherical_wave

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
real, dimension(:,:), allocatable :: nnp1
integer, dimension(:,:), allocatable :: triangle_mask
integer, dimension(:,:), allocatable :: triangle_mask1
integer, dimension(:,:), allocatable :: fourier_wave
integer, dimension(:,:), allocatable :: spherical_wave

real, allocatable, dimension(:), public :: sin_lat
real, allocatable, dimension(:), public :: cos_lat
real, allocatable, dimension(:), public :: cosm_lat
real, allocatable, dimension(:), public :: cosm2_lat
real, allocatable, dimension(:), public :: deg_lat
real, allocatable, dimension(:), public :: wts_lat
real, allocatable, dimension(:), public :: sin_hem

real, dimension(:,:,:), allocatable, public :: Pnm, Pnm_wts
real, dimension(:,:,:), allocatable, public :: Hnm, Hnm_wts

logical :: module_is_initialized = .false.

contains


!--------------------------------------------------------------------------------   
subroutine get_spherical_wave(spherical_wave_out)
!--------------------------------------------------------------------------------   
    integer, intent(out) :: spherical_wave_out(:,:)

    if(.not.module_is_initialized) &
        call mpp_error('get_spherical_wave', 'module not initialized', FATAL)

    spherical_wave_out = spherical_wave

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
                                                      lcoef_dyp, &
                                                      lnnp1

    integer, dimension(0:num_fourier,0:num_spherical) :: ltriangle_mask, &
                                                         ltriangle_mask1

    integer :: m, n, we, wo, ma, np1, nm1, nums, numf

    if(module_is_initialized) return

    nums = num_spherical; numf = num_fourier

    ltriangle_mask = 1
    ltriangle_mask1 = 1
    do n=0,num_spherical
      do m=0,num_fourier
        lfourier_wave(m,n)   = m
        lspherical_wave(m,n) = lfourier_wave(m,n) + n
        lnnp1(m,n) = lspherical_wave(m,n)*(lspherical_wave(m,n)+1.) 
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

    module_is_initialized = .true.

return
end subroutine spherical_init


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
subroutine define_legendre(lepsilon,lspherical_wave)
!--------------------------------------------------------------------------------   
    real, dimension(0:num_fourier,0:num_spherical), intent(in) :: lepsilon, lspherical_wave
    integer :: j, m, w, wo, we, mshuff, n
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
                Pnm(:,we,1) = Pnm_global(mshuff,n,js_hem:je_hem)
                Hnm(:,we,1) = Hnm_global(mshuff,n,js_hem:je_hem)
            else
                wo = wo + 1
                Pnm(:,wo,2) = Pnm_global(mshuff,n,js_hem:je_hem)
                Hnm(:,wo,2) = Hnm_global(mshuff,n,js_hem:je_hem)
            endif
        enddo
    enddo

    do j = js_hem, je_hem
        Pnm_wts(j-js_hem+1,:,:) = Pnm(j-js_hem+1,:,:)*wts_lat(2*j)
        Hnm_wts(j-js_hem+1,:,:) = Hnm(j-js_hem+1,:,:)*wts_lat(2*j)

        where(triangle_mask==0) Hnm(j,:,:) = 0.
        
         where(triangle_mask1==0) 
             Hnm_wts(j,:,:) = 0.
             Pnm_wts(j,:,:) = 0.
             Pnm(j,:,:) = 0.
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
    
    if(.not. module_is_initialized ) then
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
    
    if(.not. module_is_initialized ) then
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

    if(.not. module_is_initialized ) then
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

    if(.not. module_is_initialized ) then
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

    call compute_alpha_operator(u_cos, v_cos, +1, divergence, do_trunc)

    return
end subroutine compute_div

end module spherical_mod
