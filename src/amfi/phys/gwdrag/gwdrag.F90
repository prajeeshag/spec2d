module gwdrag_mod

use constants_mod, only : GRAV, CP_AIR, RDGAS, RVGAS

use mpp_domains_mod, only : domain2D, mpp_get_compute_domain, mpp_get_global_domain

use fms_mod, only : read_data, mpp_error, FATAL, mpp_pe, mpp_root_pe, file_exist, WARNING, NOTE

use gwdps_mod, only : gwdps

implicit none
private

public :: init_gwdrag, gwdrag

integer, parameter :: nmtvr = 14
integer :: lonf, latf
integer :: js, je, is, ie

real, allocatable :: clx(:,:,:), oa4(:,:,:), theta(:,:), sigma(:,:), gamma(:,:), &
                     elvmax(:,:), oc(:,:), hprime(:,:)

logical :: initialized=.false.

contains

subroutine init_gwdrag(domain)
    type(domain2D) :: domain
    real, allocatable :: tmp(:,:,:)
    character(len=16) :: fld
    integer :: i, k

    call mpp_get_compute_domain(domain,js,je,is,ie)
    call mpp_get_global_domain(domain,ysize=lonf,xsize=latf)

    allocate(clx(4,js:je,is:ie), oa4(4,js:je,is:ie), &
             theta(js:je,is:ie), sigma(js:je,is:ie), &
             gamma(js:je,is:ie), elvmax(js:je,is:ie), &
             oc(js:je,is:ie), hprime(js:je,is:ie))

    clx = 0.; oa4 = 0.; theta = 0.; sigma = 0.
    gamma = 0.; elvmax = 0.; elvmax = 0.; oc = 0.
    hprime = 0.

#ifndef AQUAPLANET
    if (file_exist('INPUT/mtn.nc')) then
        allocate(tmp(latf,lonf,nmtvr))   
        do i = 1, nmtvr
            write(fld,*) i-1
            fld = 'mtn'//trim(adjustl(fld))
            call read_data('INPUT/mtn.nc',trim(fld),tmp(:,:,i))
        enddo

        hprime(js:je,is:ie) = tmp(js:je,is:ie,1) 
        oc(js:je,is:ie) = tmp(js:je,is:ie,2) 

        do k = 1, 4
            oa4(k,js:je,is:ie) = tmp(js:je,is:ie,k+2) 
            clx(k,js:je,is:ie) = tmp(js:je,is:ie,k+6)
        enddo 

        theta(js:je,is:ie)  = tmp(js:je,is:ie,11)
        gamma(js:je,is:ie)  = tmp(js:je,is:ie,12)
        sigma(js:je,is:ie)  = tmp(js:je,is:ie,13)
        elvmax(js:je,is:ie) = tmp(js:je,is:ie,14)

        deallocate(tmp)
    else
        call mpp_error(FATAL,'gwdrag_mod: file INPUT/mtn.nc does not exist')
    endif
#else
    call mpp_error(NOTE,'gwdrag_mod: AQUAPLANET MODE') 
#endif
 
    initialized = .true.
 
end subroutine init_gwdrag

!--------------------------------------------------------------------------------   
subroutine gwdrag(dvdt, dudt, ugrs, vgrs, tgrs, &
                  qgrs, kpbl, prsi, del, prsl, prslk, phii, &
                  phil, dtp, dusfcg, dvsfcg)
!--------------------------------------------------------------------------------   
    real, intent(in) :: dtp
    real, dimension(:,:,:), intent(out) :: dvdt, dudt
    real, dimension(:,:,:), intent(in) :: ugrs, vgrs, tgrs, qgrs, prsi, del, prsl, &
                                        prslk, phii, phil 
    integer, dimension(:,:), intent(in) :: kpbl
    real, dimension(:,:), intent(out) :: dusfcg, dvsfcg
    
    integer :: imax, levs

    if (.not.initialized) call mpp_error(FATAL,'gwdrag_mod: not initialized!')

    levs = size(dvdt,1)
    imax = size(dvdt,2)*size(dvdt,3)

    call  gwdps_drv(imax, levs, dvdt, dudt, ugrs, vgrs, tgrs, &
            qgrs, kpbl, prsi, del, prsl, prslk, phii, &
            phil, dtp, hprime, oc, oa4, clx, theta, &
            sigma, gamma, elvmax, dusfcg, dvsfcg)

    return
end subroutine gwdrag

!--------------------------------------------------------------------------------   
subroutine gwdps_drv(imax, levs, dvdt, dudt, ugrs, vgrs, tgrs, &
            qgrs, kpbl, prsi, del, prsl, prslk, phii, &
            phil, dtp, hprime1, oc1, oa41, clx1, theta1, &
            sigma1, gamma1, elvmax1, dusfcg, dvsfcg)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: imax, levs
    real, intent(in) :: dtp
    real, dimension(levs,imax), intent(out) :: dvdt, dudt
    real, dimension(levs,imax), intent(in) :: ugrs, vgrs, tgrs, qgrs, del, prsl
    real, dimension(levs,imax), intent(in) :: prslk, phil
    real, dimension(levs+1,imax), intent(in) :: prsi, phii
    real, dimension(imax), intent(in) :: hprime1, oc1, theta1, sigma1, gamma1, elvmax1
    real, dimension(4,imax), intent(in) :: oa41, clx1
    real, dimension(imax), intent(out) :: dusfcg, dvsfcg
    integer, dimension(imax), intent(in) :: kpbl

    integer :: me, ipr, i
    logical :: lprnt

    if (.not.initialized) call mpp_error(FATAL,'gwdrag_mod: not initialized!')

    me = 1
    if(mpp_pe()==mpp_root_pe()) me = 0
    lprnt = .false.
    ipr = 0

    do i = 1, imax
        call gwdps(levs, dvdt(:,i), dudt(:,i), ugrs(:,i), vgrs(:,i), tgrs(:,i), &
            qgrs(:,i), kpbl(i), prsi(:,i), del(:,i), prsl(:,i), prslk(:,i), phii(:,i), &
            phil(:,i), dtp, hprime1(i), oc1(i), oa41(1:4,i), clx1(1:4,i), theta1(i), &
            sigma1(i), gamma1(i), elvmax1(i), dusfcg(i), dvsfcg(i), GRAV, CP_AIR, RDGAS, &
            RVGAS, lonf, nmtvr, me, lprnt, ipr)
    enddo

    return

end subroutine gwdps_drv

end module gwdrag_mod
