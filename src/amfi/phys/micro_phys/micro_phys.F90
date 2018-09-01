module micro_phys_mod

use mpp_mod, only : mpp_error, FATAL, NOTE

use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_global_domain

use fms_io_mod, only : restart_file_type, reg_rf=>register_restart_field, file_exist, &
    restore_state, save_restart

use gscond_mod, only : gscond

use precpd_mod, only : precpd 

use ocpack_mod, only : ocpack_typeP, oc_ny, oc_nlat, oc_npack, get_ocpackP

implicit none
private

public :: init_micro_phys, micro_phys, end_micro_phys, save_restart_micro_phys 

real, allocatable, dimension(:,:,:) :: tp, qp, tp1, qp1
real, allocatable, dimension(:,:) :: psp, psp1
real, allocatable, dimension(:,:) :: work1, work2

type(restart_file_type) :: rstrt
character (len=32) :: resfnm = 'amfi_mp_res'

integer :: is, ie, js, je, nlev

real :: rhbbot = 0.85, rhbtop = 0.85, rhc_max = 0.9999

logical :: initialized=.false.

contains

subroutine init_micro_phys(domain,nlev_in,lat)
    type(domain2d) :: domain
    integer, intent(in) :: nlev_in
    real, intent(in) :: lat(:,:)

    real, parameter :: dxmax=-16.118095651, dxmin=-9.800790154, &
                       dxinv=1.0/(dxmax-dxmin)
    integer :: indx, j, nlons, i
    real :: tmp1
    type(ocpack_typeP), allocatable :: ocpkP(:,:)

    nlev = nlev_in

    call mpp_get_compute_domain(domain,js,je,is,ie)

    allocate(tp(nlev,js:je,is:ie), tp1(nlev,js:je,is:ie), &
             qp(nlev,js:je,is:ie), qp1(nlev,js:je,is:ie), &
             psp(js:je,is:ie), psp1(js:je,is:ie))

    tp = 0.;  qp = 0.;  tp1 = 0.;  qp1 = 0. 
    psp = 0; psp1 = 0.;
        
    indx = reg_rf(rstrt, resfnm, 'tp', tp, domain, mandatory=.false.)
    indx = reg_rf(rstrt, resfnm, 'tp1', tp1, domain, mandatory=.false.)
    indx = reg_rf(rstrt, resfnm, 'qp', qp, domain, mandatory=.false.)
    indx = reg_rf(rstrt, resfnm, 'qp1', qp1, domain, mandatory=.false.)
    indx = reg_rf(rstrt, resfnm, 'psp', psp, domain, mandatory=.false.)
    indx = reg_rf(rstrt, resfnm, 'psp1', psp1, domain, mandatory=.false.)

    allocate(work1(js:je,is:ie), work2(js:je,is:ie))

    allocate(ocpkP(oc_npack(),oc_ny()))
    call get_ocpackP(ocpkP)

    do i = is, ie
        do j = js, je
            nlons = ocpkP(1,j)%ilen
            if (i>nlons) nlons = ocpkP(2,j)%ilen
            tmp1 = (log(cos(lat(j,i)) / (nlons*oc_nlat())) - dxmin) * dxinv
            tmp1 = max( 0.0, min( 1.0, tmp1 ) )
            work1(j,i) = tmp1
        end do
    end do
    work2 = 1.0 - work1

    if (file_exist('INPUT/'//trim(resfnm))) call restore_state(rstrt)

    initialized = .true. 
end subroutine init_micro_phys

subroutine micro_phys(dt, prsl, ps, prslk, del, q, cwm, t, rn, sn)
    real, intent(in) :: dt
    real, dimension(:,js:,is:), intent(in) :: prsl, prslk, del
    real, dimension(js:,is:), intent(in) :: ps
    real, dimension(:,js:,is:), intent(inout) :: q, cwm, t
    real, dimension(js:,is:), intent(out) :: rn, sn

    real, dimension(size(q,1),size(q,2),size(q,3)) :: rhc
   
    integer :: km, imax, k

    if(.not.initialized) call mpp_error(FATAL,'micro_phys_mod: module not initialized!')

    km = size(q,1)
    imax = size(q,2)*size(q,3)

    do k = 1, km
        rhc(k,:,:) = rhbbot - (rhbbot - rhbtop) * (1.0 - prslk(k,:,:))
        rhc(k,:,:) = rhc_max*work1(:,:) + rhc(k,:,:)*work2(:,:)
    enddo

    where(rhc>1.0) rhc = 1.0
    where(rhc<0.0) rhc = 0.0

    rn = 0.; sn = 0.

    call zhao_carr(imax, km, dt, prsl, ps, del, q, cwm, t, tp, qp, psp, &
                     tp1, qp1, psp1, rhc, rn, sn)

    return
        
end subroutine micro_phys


subroutine zhao_carr(imax, km, dt, prsl, ps, del, q, cwm, t, tpc, qpc, pspc, &
                     tpc1, qpc1, pspc1, rhc, rn, sn)

      integer, intent(in) :: imax, km
      real, intent(in) :: dt
      real, dimension(km,imax), intent(in) :: prsl, del, rhc
      real, dimension(imax), intent(in) :: ps
      real, dimension(km,imax), intent(inout) :: q, cwm, t, tpc, qpc, tpc1, qpc1
      real, dimension(imax), intent(inout) :: pspc, pspc1, rn, sn
      
      integer :: i
          
      do i = 1, imax
            call gscond (km, dt, prsl(:,i), ps(i), q(:,i), cwm(:,i), t(:,i), &
                         tpc(:,i), qpc(:,i), pspc(i), tpc1(:,i), qpc1(:,i), &
                         pspc1(i), rhc(:,i))

            call precpd (km, dt, del(:,i), prsl(:,i), ps(i), q(:,i), cwm(:,i), &
                         t(:,i), rn(i), sn(i), rhc(:,i))

      end do

      return

end subroutine zhao_carr

subroutine save_restart_micro_phys(tstamp)
    character(len=*), optional :: tstamp

    call save_restart(rstrt,tstamp)
   
end subroutine save_restart_micro_phys 

subroutine end_micro_phys()
    call save_restart(rstrt)
    return
end subroutine end_micro_phys

end module micro_phys_mod
