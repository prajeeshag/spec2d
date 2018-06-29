module spherical_mod

use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_error, FATAL, WARNING, NOTE

use fms_io_mod, only: write_data

use spherical_data_mod, only: nwaves_oe, specVar, specCoef, ms, me, nlat
use spherical_data_mod, only: ns4m, ne4m, num_spherical, num_fourier, trunc
use spherical_data_mod, only: tshuffle, legendrePol, js, je, js_hem, je_hem
use spherical_data_mod, only: iseven
use spherical_data_mod, only: operator(+), operator(-), operator(*)
use spherical_data_mod, only: operator(/), assignment(=)

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

public :: compute_lon_deriv_cos, compute_lat_deriv_cos
public :: nwaves_oe, specVar, specCoef, num_spherical, num_fourier, trunc
public :: spherical_init, compute_ucos_vcos, compute_vor_div, compute_vor
public :: compute_div, triangle_mask, do_truncation, nnp1, spherical_wave
public :: operator(+), operator(-), operator(*), operator(/), assignment(=)

type(specCoef(n=:)), allocatable :: eigen_laplacian
type(specCoef(n=:)), allocatable :: epsilon
type(specCoef(n=:)), allocatable :: r_epsilon
type(specCoef(n=:)), allocatable :: fourier_wave
type(specCoef(n=:)), allocatable :: spherical_wave
type(specCoef(n=:)), allocatable :: coef_uvm
type(specCoef(n=:)), allocatable :: coef_uvc
type(specCoef(n=:)), allocatable :: coef_uvp
type(specCoef(n=:)), allocatable :: coef_alpm
type(specCoef(n=:)), allocatable :: coef_alpp
type(specCoef(n=:)), allocatable :: coef_dym
type(specCoef(n=:)), allocatable :: coef_dx
type(specCoef(n=:)), allocatable :: coef_dyp
type(specCoef(n=:)), allocatable :: triangle_mask
type(specCoef(n=:)), allocatable :: triangle_mask1
type(specCoef(n=:)), allocatable :: nnp1

real, allocatable, dimension(:), public :: sin_lat
real, allocatable, dimension(:), public :: cos_lat
real, allocatable, dimension(:), public :: cosm_lat
real, allocatable, dimension(:), public :: cosm2_lat
real, allocatable, dimension(:), public :: deg_lat
real, allocatable, dimension(:), public :: wts_lat
real, allocatable, dimension(:), public :: sin_hem

type(legendrePol(nj=:,n=:)), allocatable, public :: Pnm, Pnm_wts
type(legendrePol(nj=:,n=:)), allocatable, public :: Hnm, Hnm_wts

logical :: module_is_initialized = .false.

contains

!--------------------------------------------------------------------------
subroutine spherical_init()

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
                                                      ltriangle_mask, &
                                                      ltriangle_mask1, &
                                                      lnnp1

    integer :: m, n, we, wo, ma, np1, nm1, nums, numf

    if(module_is_initialized) return

    nums = num_spherical; numf = num_fourier

    ltriangle_mask = 1.0
    ltriangle_mask1 = 1.0
    do n=0,num_spherical
      do m=0,num_fourier
        lfourier_wave(m,n)   = m
        lspherical_wave(m,n) = lfourier_wave(m,n) + n
        lnnp1(m,n) = lspherical_wave(m,n)*(lspherical_wave(m,n)+1.) 
        if(lspherical_wave(m,n)>trunc) ltriangle_mask(m,n) = 0.0
        if(lspherical_wave(m,n)>trunc+1) ltriangle_mask1(m,n) = 0.0
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

    allocate(specCoef(n=nwaves_oe) :: eigen_laplacian, &
                                      epsilon, &
                                      r_epsilon, &
                                      fourier_wave, &
                                      spherical_wave, &
                                      coef_uvm, &
                                      coef_uvc, &
                                      coef_uvp, &
                                      coef_alpm, &
                                      coef_alpp, &
                                      coef_dym, &
                                      coef_dx, &
                                      coef_dyp, &
                                      triangle_mask, &
                                      triangle_mask1, &
                                      nnp1)


    wo = 0
    we = 0

    do m = ms, me
        ma = tshuffle(m)
        do n = ns4m(m), ne4m(m)
            np1 = n + 1
            nm1 = n - 1
            if (mod(n,2)==0) then
                we = we + 1
                eigen_laplacian%ev(we) = leigen_laplacian(ma,n)
                epsilon%ev(we)         = lepsilon(ma,n)
                r_epsilon%ev(we)       = r_lepsilon(ma,n)
                fourier_wave%ev(we)    = lfourier_wave(ma,n)
                spherical_wave%ev(we)  = lspherical_wave(ma,n)
                coef_uvm%ev(we)        = lcoef_uvm(ma,n)
                coef_uvc%ev(we)        = lcoef_uvc(ma,n)
                coef_uvp%ev(we)        = lcoef_uvp(ma,n)
                coef_alpm%ev(we)       = lcoef_alpm(ma,n)
                coef_alpp%ev(we)       = lcoef_alpp(ma,n)
                coef_dym%ev(we)        = lcoef_dym(ma,n)
                coef_dx%ev(we)         = lcoef_dx(ma,n)
                coef_dyp%ev(we)        = lcoef_dyp(ma,n)
                triangle_mask%ev(we)   = ltriangle_mask(ma,n)
                triangle_mask1%ev(we)  = ltriangle_mask1(ma,n)
                nnp1%ev(we)            = lnnp1(ma,n)
            else
                wo = wo + 1
                eigen_laplacian%od(wo) = leigen_laplacian(ma,n)
                epsilon%od(wo)         = lepsilon(ma,n)
                r_epsilon%od(wo)       = r_lepsilon(ma,n)
                fourier_wave%od(wo)    = lfourier_wave(ma,n)
                spherical_wave%od(wo)  = lspherical_wave(ma,n)
                coef_uvm%od(wo)        = lcoef_uvm(ma,n)
                coef_uvc%od(wo)        = lcoef_uvc(ma,n)
                coef_uvp%od(wo)        = lcoef_uvp(ma,n)
                coef_alpm%od(wo)       = lcoef_alpm(ma,n)
                coef_alpp%od(wo)       = lcoef_alpp(ma,n)
                coef_dym%od(wo)        = lcoef_dym(ma,n)
                coef_dx%od(wo)         = lcoef_dx(ma,n)
                coef_dyp%od(wo)        = lcoef_dyp(ma,n)
                triangle_mask%od(wo)   = ltriangle_mask(ma,n)
                triangle_mask1%od(wo)  = ltriangle_mask1(ma,n)
                nnp1%od(wo)            = lnnp1(ma,n)
            endif
        enddo
    enddo

    where (triangle_mask%ev==0.)
        eigen_laplacian%ev = 0.
        coef_uvm%ev = 0.
        coef_uvc%ev = 0.
        coef_uvp%ev = 0.
        coef_alpm%ev = 0.
    !    coef_alpp%ev = 0.
        coef_dym%ev = 0.
        coef_dx%ev = 0.
        coef_dyp%ev = 0.
    end where

    where (triangle_mask%od==0.)
        eigen_laplacian%od = 0.
        coef_uvm%od = 0.
        coef_uvc%od = 0.
        coef_uvp%od = 0.
        coef_alpm%od = 0.
    !    coef_alpp%od = 0.
        coef_dym%od = 0.
        coef_dx%od = 0.
        coef_dyp%od = 0.
    end where

    call define_gaussian

    call define_legendre(lepsilon,lspherical_wave)

    print *, 'max(wavenumber)=', maxval(spherical_wave%od(:)), maxval(spherical_wave%ev(:))

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

    allocate(legendrePol(nj=je_hem-js_hem+1,n=nwaves_oe) :: Pnm, Pnm_wts, Hnm, Hnm_wts)

    call compute_legendre(Pnm_global, num_fourier, 1, num_spherical, sin_hem, nlat/2)

    !call write_data('legeglob','lege', Pnm_global)

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
                Pnm%ev(:,we) = Pnm_global(mshuff,n,js_hem:je_hem)
                Hnm%ev(:,we) = Hnm_global(mshuff,n,js_hem:je_hem)
            else
                wo = wo + 1
                Pnm%od(:,wo) = Pnm_global(mshuff,n,js_hem:je_hem)
                Hnm%od(:,wo) = Hnm_global(mshuff,n,js_hem:je_hem)
            endif
        enddo
    enddo

    do j = js_hem, je_hem
        Pnm_wts%ev(j-js_hem+1,:) = Pnm%ev(j-js_hem+1,:)*wts_lat(2*j)
        Pnm_wts%od(j-js_hem+1,:) = Pnm%od(j-js_hem+1,:)*wts_lat(2*j)
        Hnm_wts%ev(j-js_hem+1,:) = Hnm%ev(j-js_hem+1,:)*wts_lat(2*j)
        Hnm_wts%od(j-js_hem+1,:) = Hnm%od(j-js_hem+1,:)*wts_lat(2*j)

        where(triangle_mask%ev==0.) Hnm%ev(j-js_hem+1,:) = 0.
        where(triangle_mask%od==0.) Hnm%od(j-js_hem+1,:) = 0.
        
        where(triangle_mask1%ev==0.) 
            Hnm_wts%ev(j-js_hem+1,:) = 0.
            Pnm_wts%ev(j-js_hem+1,:) = 0.
            Pnm%ev(j-js_hem+1,:) = 0.
        end where

        where(triangle_mask1%od==0.)
            Hnm_wts%od(j-js_hem+1,:) = 0.
            Pnm_wts%od(j-js_hem+1,:) = 0.
            Pnm%od(j-js_hem+1,:) = 0.
        end where
    enddo

    return
end subroutine define_legendre


!---------------------------------------------------------------------------
subroutine compute_lon_deriv_cos(spherical, deriv_lon)
!---------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in) :: spherical
    type(specVar(nlev=*,n=*)), intent(out) :: deriv_lon
    
    integer :: k, n
    
    if(.not. module_is_initialized ) then
      call mpp_error('compute_lon_deriv','module spherical not initialized', FATAL)
    end if
    
    do n = 1, spherical%n
       deriv_lon%ev(:,n) = coef_dx%ev(n)*cmplx(-aimag(spherical%ev(:,n)),real(spherical%ev(:,n)))
       deriv_lon%od(:,n) = coef_dx%od(n)*cmplx(-aimag(spherical%od(:,n)),real(spherical%od(:,n)))
    end do

    return
end subroutine compute_lon_deriv_cos



!--------------------------------------------------------------------------
subroutine compute_lat_deriv_cos(spherical,deriv_lat)
!--------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in) :: spherical
    type(specVar(nlev=*,n=*)), intent(out) :: deriv_lat
    
    integer :: k, n, nw
    
    if(.not. module_is_initialized ) then
      call mpp_error('compute_lat_deriv','module spherical not initialized', FATAL)
    end if

    nw = spherical%n   
 
    deriv_lat%ev(:,:) = cmplx(0.,0.)
    deriv_lat%od(:,:) = cmplx(0.,0.)

    do k = 1, spherical%nlev

       deriv_lat%ev(k,2:nw) = &  
            - spherical%od(k,1:nw-1)*coef_dym%od(1:nw-1)

       deriv_lat%ev(k,1:nw) =  deriv_lat%ev(k,1:nw)   &
            + spherical%od(k,1:nw)*coef_dyp%od(1:nw)

       deriv_lat%od(k,1:nw) = &  
            - spherical%ev(k,1:nw)*coef_dym%ev(1:nw)

       deriv_lat%od(k,1:nw-1) =  deriv_lat%od(k,1:nw-1)   &
            + spherical%ev(k,2:nw)*coef_dyp%ev(2:nw)

       !deriv_lat(:,1:num_spherical,k) = &  
       !     - spherical(:,0:num_spherical-1,k)*coef_dym(:,1:num_spherical)

       !deriv_lat(:,0:num_spherical-1,k) =  deriv_lat(:,0:num_spherical-1,k)   &
       !     + spherical(:,1:num_spherical,k)*coef_dyp(:,0:num_spherical-1)

    end do

    return
end subroutine compute_lat_deriv_cos


!--------------------------------------------------------------------------------   
subroutine do_truncation(spherical)
!--------------------------------------------------------------------------------   
    type(specVar(nlev=*,n=*)), intent(inout) :: spherical

    integer :: k

    do k = 1, spherical%nlev
        where(triangle_mask%ev==0.) spherical%ev(k,:) = 0.
        where(triangle_mask%od==0.) spherical%od(k,:) = 0.
    enddo

    return

end subroutine do_truncation



!----------------------------------------------------------------------
subroutine compute_ucos_vcos(vorticity , divergence, u_cos, v_cos, do_trunc)
!----------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in)  :: vorticity
    type(specVar(nlev=*,n=*)), intent(in)  :: divergence
    type(specVar(nlev=*,n=*)), intent(out) :: u_cos
    type(specVar(nlev=*,n=*)), intent(out) :: v_cos
    logical, intent(in), optional :: do_trunc

    logical :: do_trunc1
    integer :: k, nw

    if(.not. module_is_initialized ) then
      call mpp_error('compute_ucos_vcos','module spherical not initialized', FATAL)
    end if

    do_trunc1 = .true.
   
    if(present(do_trunc)) do_trunc1=do_trunc 

    nw = vorticity%n

    do k=1,vorticity%nlev

       u_cos%ev(k,:) = coef_uvc%ev(:)*                                     &
            cmplx(-aimag(divergence%ev(k,:)),real(divergence%ev(k,:)))

       u_cos%od(k,:) = coef_uvc%od(:)*                                     &
            cmplx(-aimag(divergence%od(k,:)),real(divergence%od(k,:)))

       u_cos%ev(k,2:nw) = u_cos%ev(k,2:nw) +         &
            coef_uvm%od(1:nw-1)*vorticity%od(k,1:nw-1)

       u_cos%od(k,1:nw) = u_cos%od(k,1:nw) +         &
            coef_uvm%ev(1:nw)*vorticity%ev(k,1:nw)

       u_cos%ev(k,1:nw) = u_cos%ev(k,1:nw) -     &
            coef_uvp%od(1:nw)*vorticity%od(k,1:nw)

       u_cos%od(k,1:nw-1) = u_cos%od(k,1:nw-1) -     &
            coef_uvp%ev(2:nw)*vorticity%ev(k,2:nw)

       v_cos%ev(k,:) = coef_uvc%ev(:)*                                     &
            cmplx(-aimag(vorticity%ev(k,:)),real(vorticity%ev(k,:)))

       v_cos%od(k,:) = coef_uvc%od(:)*                                     &
            cmplx(-aimag(vorticity%od(k,:)),real(vorticity%od(k,:)))

       v_cos%ev(k,2:nw) = v_cos%ev(k,2:nw) -         &
            coef_uvm%od(1:nw-1)*divergence%od(k,1:nw-1)

       v_cos%od(k,1:nw) = v_cos%od(k,1:nw) -         &
            coef_uvm%ev(1:nw)*divergence%ev(k,1:nw)

       v_cos%ev(k,1:nw) = v_cos%ev(k,1:nw) +     &
            coef_uvp%od(1:nw)*divergence%od(k,1:nw)

       v_cos%od(k,1:nw-1) = v_cos%od(k,1:nw-1) +     &
            coef_uvp%ev(2:nw)*divergence%ev(k,2:nw)

    end do

    if (do_trunc1) then
       call do_truncation(v_cos)
       call do_truncation(u_cos)
    endif
            
    return
end subroutine compute_ucos_vcos

!--------------------------------------------------------------------------------
subroutine compute_alpha_operator(spherical_a, spherical_b, rsign, alpha, do_trunc)
!--------------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in)  :: spherical_a
    type(specVar(nlev=*,n=*)), intent(in)  :: spherical_b
    type(specVar(nlev=*,n=*)), intent(out) :: alpha
    integer, intent(in) :: rsign
    logical, intent(in), optional :: do_trunc

    logical :: do_trunc1
    integer :: k, nw

    if(.not. module_is_initialized ) then
      call mpp_error('compute_vor or div','module spherical not initialized', FATAL)
    end if

    do_trunc1 = .true.
    if(present(do_trunc)) do_trunc1=do_trunc
    
    alpha%ev = cmplx(0.,0.)
    alpha%od = cmplx(0.,0.)
   
    nw = spherical_a%n 

    do k = 1, spherical_a%nlev
       alpha%ev(k,:) = coef_dx%ev(:)*    &
            cmplx(-aimag(spherical_a%ev(k,:)),real(spherical_a%ev(k,:)))

       alpha%od(k,:) = coef_dx%od(:)*    &
            cmplx(-aimag(spherical_a%od(k,:)),real(spherical_a%od(k,:)))

       alpha%ev(k,2:nw) = alpha%ev(k,2:nw) -  &
            rsign*coef_alpm%od(1:nw-1)  &
            *spherical_b%od(k,1:nw-1)

       alpha%od(k,1:nw) = alpha%od(k,1:nw) -  &
            rsign*coef_alpm%ev(1:nw)  &
            *spherical_b%ev(k,1:nw)

       alpha%ev(k,1:nw) = alpha%ev(k,1:nw) +  &
            rsign*coef_alpp%od(1:nw)*spherical_b%od(k,1:nw)

       alpha%od(k,1:nw-1) = alpha%od(k,1:nw-1) +  &
            rsign*coef_alpp%ev(2:nw)*spherical_b%ev(k,2:nw)

    end do

    if (do_trunc1) call do_truncation(alpha)

    return
end subroutine compute_alpha_operator

!-------------------------------------------------------------------------
subroutine compute_vor_div(u_cos, v_cos, vorticity, divergence, do_trunc)
!-------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in)  :: u_cos
    type(specVar(nlev=*,n=*)), intent(in)  :: v_cos
    type(specVar(nlev=*,n=*)), intent(out) :: vorticity
    type(specVar(nlev=*,n=*)), intent(out) :: divergence
    logical, intent(in), optional :: do_trunc

    call compute_alpha_operator(v_cos, u_cos, -1, vorticity, do_trunc)
    call compute_alpha_operator(u_cos, v_cos, +1, divergence, do_trunc)

    return
end subroutine compute_vor_div

!-------------------------------------------------------------------------
subroutine compute_vor(u_cos, v_cos, vorticity, do_trunc)
!-------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in)  :: u_cos
    type(specVar(nlev=*,n=*)), intent(in)  :: v_cos
    type(specVar(nlev=*,n=*)), intent(out) :: vorticity
    logical, intent(in), optional :: do_trunc

    call compute_alpha_operator(v_cos, u_cos, -1, vorticity, do_trunc)

    return
end subroutine compute_vor

!-------------------------------------------------------------------------
subroutine compute_div(u_cos, v_cos, divergence, do_trunc)
!-------------------------------------------------------------------------

    type(specVar(nlev=*,n=*)), intent(in)  :: u_cos
    type(specVar(nlev=*,n=*)), intent(in)  :: v_cos
    type(specVar(nlev=*,n=*)), intent(out) :: divergence
    logical, intent(in), optional :: do_trunc

    call compute_alpha_operator(u_cos, v_cos, +1, divergence, do_trunc)

    return
end subroutine compute_div

!!-------------------------------------------------------------------
!subroutine compute_gradient_cos_3d(spherical, deriv_lon, deriv_lat) 
!!-------------------------------------------------------------------
!
!complex, intent(in), dimension (:,:,:) :: spherical
!complex, intent(out), dimension (:,:,:) :: deriv_lat
!complex, intent(out), dimension (:,:,:) :: deriv_lon
!
!deriv_lon = compute_lon_deriv_cos(spherical)
!deriv_lat = compute_lat_deriv_cos(spherical)
!
!return
!end subroutine compute_gradient_cos_3d
!
!!----------------------------------------------------------------------
!function compute_laplacian_3d(spherical, power) result(laplacian)
!!----------------------------------------------------------------------
!
!  complex, intent(in), dimension (:,:,:) :: spherical
!  integer, optional :: power
!
!  complex, dimension (size(spherical,1), size(spherical,2), size(spherical,3)) :: laplacian
!
!  integer :: k
!  real, dimension(size(spherical,1), size(spherical,2)) :: factor
!
!  if(.not. module_is_initialized ) then
!      call mpp_error('compute_laplacian','module spherical not initialized', FATAL)
!  end if
!
!  if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
!      if(present(power)) then 
!          if(power >= 0) then
!              factor = (-eigen_laplacian)**power
!          else 
!              where (eigen_laplacian .ne. 0.0) 
!                  factor = (-eigen_laplacian)**power
!              else where
!                  factor = 0.0
!              end where
!          end if
!      else
!          factor = -eigen_laplacian
!      end if
!  else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
!      if(present(power)) then 
!          if(power >= 0) then
!              factor = (-eigen_laplacian(ms:me,ns:ne))**power
!          else 
!              where (eigen_laplacian(ms:me,ns:ne) .ne. 0.0) 
!                  factor = (-eigen_laplacian(ms:me,ns:ne))**power
!              else where
!                  factor = 0.0
!              end where
!          end if
!      else
!          factor = -eigen_laplacian(ms:me,ns:ne)
!      end if
!  else
!      call mpp_error( 'compute_laplacian', 'invalid argument size', FATAL )
!  endif
!
!  do k= 1,size(spherical,3)
!     laplacian(:,:,k) = spherical(:,:,k)*factor
!  end do
!
!  return
!end function compute_laplacian_3d

!!-----------------------------------------------------------------------
!subroutine triangular_truncation_3d(spherical, trunc)
!!-----------------------------------------------------------------------
!
!complex, intent(inout), dimension (:,:,:) :: spherical
!integer, intent(in), optional :: trunc
!integer :: k
!
!if( size(spherical,1).EQ.num_fourier+1 .AND. size(spherical,2).EQ.num_spherical+1 )then
!    if(present(trunc)) then
!        do k=1, size(spherical,3)
!           where (spherical_wave > trunc)
!               spherical(:,:,k) = cmplx(0.,0.)
!           end where
!        end do
!    else
!        do k=1, size(spherical,3)
!           spherical(:,:,k) = spherical(:,:,k)*triangle_mask
!        end do
!    end if
!else if( size(spherical,1).EQ.me-ms+1 .AND. size(spherical,2).EQ.ne-ns+1 )then
!    if(present(trunc)) then
!        do k=1, size(spherical,3)
!           where (spherical_wave(ms:me,ns:ne) > trunc)
!               spherical(:,:,k) = cmplx(0.,0.)
!           end where
!        end do
!    else
!        do k=1, size(spherical,3)
!           spherical(:,:,k) = spherical(:,:,k)*triangle_mask(ms:me,ns:ne)
!        end do
!    end if
!else
!    call mpp_error( 'triang_trunc', 'invalid argument size', FATAL )
!endif
!
!return
!end subroutine triangular_truncation_3d


!rest are 2d versions of 3d routines

!!-------------------------------------------------------------------
!subroutine compute_gradient_cos_2d(spherical, deriv_lon, deriv_lat) 
!!-------------------------------------------------------------------
!
!complex, intent(in),  dimension(:,:) :: spherical
!complex, intent(out), dimension(size(spherical,1), size(spherical,2)) :: deriv_lat
!complex, intent(out), dimension(size(spherical,1), size(spherical,2)) :: deriv_lon
!
!complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
!complex, dimension(size(spherical,1), size(spherical,2), 1) :: deriv_lat_3d
!complex, dimension(size(spherical,1), size(spherical,2), 1) :: deriv_lon_3d
!
!spherical_3d(:,:,1) = spherical(:,:)
!call compute_gradient_cos_3d(spherical_3d, deriv_lon_3d, deriv_lat_3d)
!deriv_lon(:,:) = deriv_lon_3d(:,:,1)
!deriv_lat(:,:) = deriv_lat_3d(:,:,1)
!
!return
!end subroutine compute_gradient_cos_2d
!
!!----------------------------------------------------------------------
!function compute_laplacian_2d(spherical, power) result(laplacian)
!!----------------------------------------------------------------------
!
!complex, intent(in), dimension (:,:) :: spherical
!integer, optional :: power
!
!complex, dimension (size(spherical,1), size(spherical,2)) :: laplacian
!
!complex, dimension(size(spherical,1), size(spherical,2), 1) :: spherical_3d
!complex, dimension(size(spherical,1), size(spherical,2), 1) :: laplacian_3d
!
!spherical_3d(:,:,1) = spherical(:,:)
!laplacian_3d = compute_laplacian_3d(spherical_3d, power)
!laplacian(:,:) = laplacian_3d(:,:,1)
!
!return
!end function compute_laplacian_2d
!
!!----------------------------------------------------------------------
!subroutine compute_ucos_vcos_2d(vorticity , divergence, u_cos, v_cos)
!!----------------------------------------------------------------------
!
!complex, intent(in),  dimension (:,:) :: vorticity
!complex, intent(in),  dimension (size(vorticity,1), size(vorticity,2)) :: divergence
!complex, intent(out), dimension (size(vorticity,1), size(vorticity,2)) :: u_cos
!complex, intent(out), dimension (size(vorticity,1), size(vorticity,2)) :: v_cos
!
!complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: vorticity_3d
!complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: divergence_3d
!complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: u_cos_3d
!complex, dimension (size(vorticity,1), size(vorticity,2), 1) :: v_cos_3d
!
!vorticity_3d(:,:,1)  = vorticity(:,:)
!divergence_3d(:,:,1) = divergence(:,:)
!call compute_ucos_vcos_3d(vorticity_3d, divergence_3d, u_cos_3d, v_cos_3d)
!u_cos(:,:) = u_cos_3d(:,:,1)
!v_cos(:,:) = v_cos_3d(:,:,1)
!
!return
!end subroutine compute_ucos_vcos_2d
!
!!-------------------------------------------------------------------------
!subroutine compute_vor_div_2d(u_div_cos, v_div_cos, vorticity, divergence)
!!-------------------------------------------------------------------------
!
!complex, intent(in),  dimension (:,:) :: u_div_cos
!complex, intent(in),  dimension (size(u_div_cos,1), size(u_div_cos,2)) :: v_div_cos
!complex, intent(out), dimension (size(u_div_cos,1), size(u_div_cos,2)) :: vorticity
!complex, intent(out), dimension (size(u_div_cos,1), size(u_div_cos,2)) :: divergence
!
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: vorticity_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: divergence_3d
!
!u_div_cos_3d(:,:,1) = u_div_cos(:,:)
!v_div_cos_3d(:,:,1) = v_div_cos(:,:)
!call compute_vor_div_3d(u_div_cos_3d, v_div_cos_3d, vorticity_3d, divergence_3d)
!vorticity(:,:)  = vorticity_3d(:,:,1)
!divergence(:,:) = divergence_3d(:,:,1)
!
!return
!end subroutine compute_vor_div_2d
!
!!-------------------------------------------------------------------------
!function compute_vor_2d(u_div_cos, v_div_cos) result(vorticity)
!!-------------------------------------------------------------------------
!
!complex, intent(in),  dimension (:,:) :: u_div_cos
!complex, intent(in),  dimension (:,:) :: v_div_cos
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2)) :: vorticity
!
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: vorticity_3d
!
!u_div_cos_3d(:,:,1) = u_div_cos(:,:)
!v_div_cos_3d(:,:,1) = v_div_cos(:,:)
!vorticity_3d = compute_vor_3d(u_div_cos_3d, v_div_cos_3d)
!vorticity(:,:) = vorticity_3d(:,:,1)
!
!return
!end function compute_vor_2d
!
!!-------------------------------------------------------------------------
!function compute_div_2d(u_div_cos, v_div_cos) result(divergence)
!!-------------------------------------------------------------------------
!
!complex, intent(in),  dimension (:,:) :: u_div_cos
!complex, intent(in),  dimension (:,:) :: v_div_cos
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2)) :: divergence
!
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: u_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: v_div_cos_3d
!complex, dimension (size(u_div_cos,1), size(u_div_cos,2), 1) :: divergence_3d
!
!u_div_cos_3d(:,:,1) = u_div_cos(:,:)
!v_div_cos_3d(:,:,1) = v_div_cos(:,:)
!divergence_3d = compute_div_3d(u_div_cos_3d, v_div_cos_3d)
!divergence(:,:) = divergence_3d(:,:,1)
!
!return
!end function compute_div_2d
!
!!--------------------------------------------------------------------------------
!function compute_alpha_operator_2d(spherical_a, spherical_b, rsign) result(alpha)
!!--------------------------------------------------------------------------------
!
!complex, intent(in),  dimension (:,:) :: spherical_a
!complex, intent(in),  dimension (:,:) :: spherical_b
!integer, intent(in) :: rsign
!
!complex, dimension (size(spherical_a,1), size(spherical_a,2)) :: alpha
!
!complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: spherical_a_3d
!complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: spherical_b_3d
!complex, dimension (size(spherical_a,1), size(spherical_a,2), 1) :: alpha_3d
!
!spherical_a_3d(:,:,1) = spherical_a(:,:)
!spherical_b_3d(:,:,1) = spherical_b(:,:)
!alpha_3d = compute_alpha_operator_3d(spherical_a_3d, spherical_b_3d, rsign)
!alpha(:,:) = alpha_3d(:,:,1)
!
!return
!end function compute_alpha_operator_2d
!
!!-----------------------------------------------------------------------
!subroutine triangular_truncation_2d(spherical, trunc)
!!-----------------------------------------------------------------------
!
!complex, intent(inout), dimension (:,:) :: spherical
!integer, intent(in), optional :: trunc
!complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d
!
!spherical_3d(:,:,1) = spherical(:,:)
!call triangular_truncation_3d(spherical_3d, trunc)
!spherical(:,:) = spherical_3d(:,:,1)
!
!return
!end subroutine triangular_truncation_2d
!
!!--------------------------------------------------------------------------------
!subroutine rhomboidal_truncation_2d(spherical, trunc_fourier, trunc_spherical)
!!--------------------------------------------------------------------------------
!
!complex, intent(inout), dimension (:,:) :: spherical
!integer, intent(in), optional :: trunc_fourier, trunc_spherical
!
!complex, dimension (size(spherical,1), size(spherical,2), 1) :: spherical_3d
!
!spherical_3d(:,:,1) = spherical(:,:)
!call rhomboidal_truncation_3d(spherical_3d, trunc_fourier, trunc_spherical)
!spherical(:,:) = spherical_3d(:,:,1)
!
!return
!end subroutine rhomboidal_truncation_2d
!
!!-----------------------------------------------------------------------
!subroutine spherical_end
!
!if(.not.module_is_initialized) return
!
!deallocate(eigen_laplacian, epsilon, r_epsilon, fourier_wave, spherical_wave)
!deallocate(coef_uvm, coef_uvc, coef_uvp, coef_alpm, coef_alpp, coef_dym, coef_dx, coef_dyp, triangle_mask)
!module_is_initialized = .false.
!
!return
!end subroutine spherical_end
!!-----------------------------------------------------------------------

end module spherical_mod
