module implicit_mod

use, intrinsic :: iso_c_binding

use spherical_mod, only : spherical_wave, nnp1

use constants_mod, only : rd => rdgas, cp => cp_air
use constants_mod, only : rearth => radius

use mpp_mod, only : mpp_error, FATAL

use fms_io_mod, only : write_data
implicit none
private

public :: init_implicit, do_implicit 

real, allocatable :: bmhyb(:,:), amhyb(:,:)
real, allocatable :: svhyb(:), d_hyb_m(:,:,:)
real, allocatable :: tor_hyb(:), dm205_hyb(:,:,:)
real, allocatable :: ak5(:), bk5(:)

integer :: levs, levp1

logical :: initialized=.false.

contains

!--------------------------------------------------------------------------------   
subroutine init_implicit(ak, bk, ref_temp, dt, trunc)
!--------------------------------------------------------------------------------   
    real,    intent(in) :: ak(:), bk(:), ref_temp, dt
    integer, intent(in) :: trunc
 
    integer :: jcap1
   
    levs = size(ak,1)-1
    levp1 = levs + 1
    jcap1 = trunc + 3
 
    allocate(tor_hyb(levs))
    allocate(svhyb(levs))
    allocate(amhyb(levs,levs))
    allocate(bmhyb(levs,levs))
 
    allocate(d_hyb_m(levs,levs,jcap1))
    allocate(dm205_hyb(jcap1,levs,levs))
    allocate(ak5(levp1),bk5(levp1))

    ak5 = ak; bk5 = bk
 
    call am_bm_hyb
    
    call get_cd_hyb(dt,trunc+3)

    call write_data('rimplicit','bmhyb',bmhyb)
    call write_data('rimplicit','amhyb',amhyb)
    call write_data('rimplicit','svhyb',svhyb)
    call write_data('rimplicit','tor_hyb',tor_hyb)
    call write_data('rimplicit','d_hyb_m',d_hyb_m)

    initialized = .true.

end subroutine init_implicit


!--------------------------------------------------------------------------------   
subroutine do_implicit(div, tem, ps, div_n, tem_n, ps_n, &
                               div_dt, tem_dt, ps_dt, dt)
!--------------------------------------------------------------------------------   
    implicit none
    complex, dimension(:,:,:), intent(in) :: div, tem, div_n, tem_n, ps, ps_n 
    complex, dimension(:,:,:), intent(inout) :: div_dt, tem_dt, ps_dt
    real,                      intent(in) :: dt

    integer :: i

    if (.not.initialized) call mpp_error('do_implicit', 'module not initialized', FATAL)

    do i = 1, 2
        call implicit_corr_drv(div(:,:,i), tem(:,:,i), ps(:,:,i), div_n(:,:,i), tem_n(:,:,i), ps_n(:,:,i), &
                 div_dt(:,:,i), tem_dt(:,:,i), ps_dt(:,:,i), nnp1(:,i), spherical_wave(:,i), dt)
    enddo

end subroutine do_implicit


!--------------------------------------------------------------------------------   
subroutine implicit_corr_drv(div, tem, ps, div_n, tem_n, ps_n, &
                             div_dt, tem_dt, ps_dt, snnp1ev, ndexev, dt)
!--------------------------------------------------------------------------------   
    implicit none
    complex, dimension(:,:), intent(in) :: div, tem, div_n, tem_n
    complex, dimension(:,:), intent(in) :: ps, ps_n
    integer, dimension(:),   intent(in) :: ndexev
    real,    dimension(:),   intent(in) :: snnp1ev 
    real,                    intent(in) :: dt
    complex, dimension(:,:), intent(inout) :: div_dt, tem_dt
    complex, dimension(:,:), intent(inout) :: ps_dt

    integer :: nwaves
    type(C_PTR) :: pcdiv, pctem, pcps
    type(C_PTR) :: pcdiv_n, pctem_n, pcps_n
    type(C_PTR) :: pcdiv_dt, pctem_dt, pcps_dt

    real, pointer, dimension(:,:,:) :: pdiv, ptem, pdiv_n, ptem_n, pdiv_dt, ptem_dt
    real, pointer, dimension(:,:) :: pps, pps_n, pps_dt

    real, dimension(size(div,2),2,levs) :: de, te, de_n, te_n
    real, dimension(size(div,2),2)      :: qe, qe_n
    real, dimension(size(div,2),2,levs) :: xe, ye
    real, dimension(size(div,2),2)      :: ze

    integer :: i, j, k

    nwaves=size(div,2)

    pcdiv = C_LOC(div)
    pctem = C_LOC(tem)
    pcps = C_LOC(ps)
    pcdiv_n = C_LOC(div_n)
    pctem_n = C_LOC(tem_n)
    pcps_n = C_LOC(ps_n)
    pcdiv_dt = C_LOC(div_dt)
    pctem_dt = C_LOC(tem_dt)
    pcps_dt = C_LOC(ps_dt)

    call c_f_pointer(pcdiv,pdiv,[2,levs,nwaves])
    call c_f_pointer(pctem,ptem,[2,levs,nwaves])
    call c_f_pointer(pcps,pps,[2,nwaves])
    call c_f_pointer(pcdiv_n,pdiv_n,[2,levs,nwaves])
    call c_f_pointer(pctem_n,ptem_n,[2,levs,nwaves])
    call c_f_pointer(pcps_n,pps_n,[2,nwaves])
    call c_f_pointer(pcdiv_dt,pdiv_dt,[2,levs,nwaves])
    call c_f_pointer(pctem_dt,ptem_dt,[2,levs,nwaves])
    call c_f_pointer(pcps_dt,pps_dt,[2,nwaves])

    do i = 1, nwaves 
        de(i,1:2,1:levs)   = pdiv(1:2,1:levs,i)
        te(i,1:2,1:levs)   = ptem(1:2,1:levs,i)
        qe(i,1:2)          = pps(1:2,i)
        de_n(i,1:2,1:levs) = pdiv_n(1:2,1:levs,i)
        te_n(i,1:2,1:levs) = ptem_n(1:2,1:levs,i)
        qe_n(i,1:2)        = pps_n(1:2,i)
        xe(i,1:2,1:levs)   = pdiv_dt(1:2,1:levs,i)
        ye(i,1:2,1:levs)   = ptem_dt(1:2,1:levs,i)
        ze(i,1:2)          = pps_dt(1:2,i)
    enddo

    call implicit_corr(de,te,qe,de_n,te_n,qe_n,xe,ye,ze,snnp1ev,ndexev,nwaves,dt)

    do i = 1, nwaves 
        pdiv_dt(1:2,1:levs,i) = xe(i,1:2,1:levs)
        ptem_dt(1:2,1:levs,i) = ye(i,1:2,1:levs)
        pps_dt(1:2,i)         = ze(i,1:2)
    enddo
    
    return

end subroutine implicit_corr_drv


!--------------------------------------------------------------------------------   
subroutine implicit_corr(de,te,qe,de_n,te_n,qe_n,xe,ye,ze,snnp1ev,ndexev,nwaves,dt)
!--------------------------------------------------------------------------------   
    implicit none

    real, dimension(nwaves,2,levs), intent(in) :: de,te,de_n,te_n
    real, dimension(nwaves,2),      intent(in) :: qe,qe_n
    integer, dimension(nwaves),     intent(in) :: ndexev
    real, dimension(nwaves),        intent(in) :: snnp1ev 
    real,                           intent(in) :: dt
    integer,                        intent(in) :: nwaves
    real, dimension(nwaves,2,levs), intent(inout) :: xe,ye
    real, dimension(nwaves,2),      intent(inout) :: ze
    
    real :: ue(nwaves,2,levs), ve(nwaves,2,levs)
    real :: qdtze(nwaves,2), elne(nwaves,2,levs)
    real :: svhybdt, u1, u2
    real, parameter :: cons0 = 0., cons1 = 1., cons2 = 2.
    integer :: i,indev,indev1,indev2,j,k,l,locl,n

    indev1 = 1; indev2 = nwaves

    do j=1,levs
        do k=1,levs,2
            do indev = indev1 , indev2
                ye(indev,1,j) = ye(indev,1,j) + de_n(indev,1,k  )*bmhyb(j,k  ) + de_n(indev,1,k+1)*bmhyb(j,k+1)
                ye(indev,2,j) = ye(indev,2,j) + de_n(indev,2,k  )*bmhyb(j,k  ) + de_n(indev,2,k+1)*bmhyb(j,k+1)
            enddo
         enddo
    enddo

    do k=1,levs
        do indev = indev1 , indev2
            ze(indev,1) = ze(indev,1) + de_n(indev,1,k)*svhyb(k)
            ze(indev,2) = ze(indev,2) + de_n(indev,2,k)*svhyb(k)
        enddo
    enddo

    do indev = indev1 , indev2
        qdtze(indev,1)   =  qe(indev,1)-qe_n(indev,1) + dt*ze(indev,1)
        qdtze(indev,2)   =  qe(indev,2)-qe_n(indev,2) + dt*ze(indev,2)
    enddo
 
    do k=1,levs
        do indev = indev1 , indev2
            elne(indev,1,k) = te(indev,1,k)-te_n(indev,1,k) + dt*  ye(indev,1,k)
            elne(indev,2,k) = te(indev,2,k)-te_n(indev,2,k) + dt*  ye(indev,2,k)
        enddo
    enddo

    call dgemm ('n', 't', indev2-indev1+1, levs, levs, cons1, elne(indev1,1,1), nwaves*2, &
                    amhyb(1,1), levs, cons0, ve(indev1,1,1), nwaves*2)
    call dgemm ('n', 't', indev2-indev1+1, levs, levs, cons1, elne(indev1,2,1), nwaves*2, &
                    amhyb(1,1), levs, cons0, ve(indev1,2,1), nwaves*2)

    do j=1,levs
        do indev = indev1 , indev2
            ve(indev,1,j) = ve(indev,1,j) + tor_hyb(j)*qdtze(indev,1)
            ve(indev,1,j) = ve(indev,1,j) * snnp1ev(indev)
            ve(indev,1,j) = ve(indev,1,j) + xe(indev,1,j)
 
            ue(indev,1,j) = de(indev,1,j) + ve(indev,1,j)*dt
 
            ve(indev,2,j) = ve(indev,2,j) + tor_hyb(j)*qdtze(indev,2)
 
            ve(indev,2,j) = ve(indev,2,j) * snnp1ev(indev)
 
            ve(indev,2,j) = ve(indev,2,j) + xe(indev,2,j)
 
            ue(indev,2,j) = de(indev,2,j) + ve(indev,2,j)*dt
        enddo
    enddo

    do indev = indev1 , indev2
        call dgemm ('n', 't', 1, levs, levs, cons1, ue(indev,1,1), nwaves*2, &
                       d_hyb_m(1,1,ndexev(indev)+1), levs, cons0, ve(indev,1,1), nwaves*2)
        call dgemm ('n', 't', 1, levs, levs, cons1, ue(indev,2,1), nwaves*2, &
                       d_hyb_m(1,1,ndexev(indev)+1), levs, cons0, ve(indev,2,1), nwaves*2)
    enddo

    call dgemm ('n', 't', indev2-indev1+1, 1, levs, dt, ve(indev1,1,1), nwaves*2, &
                       svhyb(1), 1, cons0, ue(indev1,1,1), nwaves*2)
    call dgemm ('n', 't', indev2-indev1+1, 1, levs, dt, ve(indev1,2,1), nwaves*2, &
                       svhyb(1), 1, cons0, ue(indev1,2,1), nwaves*2)

    do indev = indev1 , indev2
        qdtze(indev,1) = qdtze(indev,1) + qe_n(indev,1) - ue(indev,1,1)
 
        ze(indev,1) = cons2*qdtze(indev,1) - qe(indev,1)
 
        qdtze(indev,2) = qdtze(indev,2) + qe_n(indev,2) - ue(indev,2,1)
 
        ze(indev,2) = cons2*qdtze(indev,2) - qe(indev,2)
    enddo

    call dgemm ('n', 't', indev2-indev1+1, levs, levs, cons1, ve(indev1,1,1), nwaves*2, &
                    bmhyb(1,1), levs, cons0, ue(indev1,1,1), nwaves*2)
    call dgemm ('n', 't', indev2-indev1+1, levs, levs, cons1, ve(indev1,2,1), nwaves*2, &
                    bmhyb(1,1), levs, cons0, ue(indev1,2,1), nwaves*2)

    do j=1,levs
        do indev = indev1 , indev2
            u1 = elne(indev,1,j) - dt * ue(indev,1,j) + te_n(indev,1,j)
 
            ye(indev,1,j) = cons2*u1-te(indev,1,j) 
 
            xe(indev,1,j) = cons2*ve(indev,1,j) - de(indev,1,j)
 
            u2 = elne(indev,2,j) - dt * ue(indev,2,j) + te_n(indev,2,j)
 
            ye(indev,2,j) = cons2*u2-te(indev,2,j) 
 
            xe(indev,2,j) = cons2*ve(indev,2,j) - de(indev,2,j)
        enddo
    enddo

    return
end subroutine implicit_corr



!--------------------------------------------------------------------------------   
SUBROUTINE AM_BM_hyb
!--------------------------------------------------------------------------------   
 
    implicit none 
    
    real :: pk5ref(levp1),beta,dpkref(levs), tref(levs),psref,kappa,factor, &
            alfaref(levs), vecm(levs),   yecm(levs,levs),tecm(levs,levs)
    integer :: k,kk,j,irow,icol,icolbeg,icolend
    real :: ref_temp = 300.
    
    do k=1,levs
     tref(k)=ref_temp
    enddo
    psref=80.
    beta=1.
    kappa=rd/cp
    
    
    
    
    do k=1,levp1
     pk5ref(k)=ak5(k)+bk5(k)*psref
    enddo
    
    
    do k=1,levs
     dpkref(k)=pk5ref(k+1)-pk5ref(k)
     tor_hyb(k)=beta*rd*tref(k)/(rearth*rearth)
    enddo
    
    
    alfaref(1)=log(2.) ! could also be=1.  but watch for layer values
    
    do k=2,levs
     alfaref(k)=1.-(pk5ref(k)/dpkref(k))*log(pk5ref(k+1)/pk5ref(k))
    enddo
    
    
    do k=1,levs
    enddo
    
     yecm=0.
     do irow=1,levs
        yecm(irow,irow)=alfaref(irow)*rd
        icolbeg=irow+1
        if(icolbeg.le.levs)then
         do icol=icolbeg,levs
          yecm(irow,icol)=rd*log( pk5ref(icol+1)/pk5ref(icol) )
         enddo
        endif
     enddo
    
    tecm=0.
    
    do irow=1,levs
       tecm(irow,irow)=kappa*tref(irow)*alfaref(irow)
       icolend=irow-1
    
    
    do icol=1,icolend
    factor=(kappa*tref(irow)/ dpkref(irow))*log(pk5ref(irow+1)/pk5ref(irow))
    tecm(irow,icol)=factor*dpkref(icol)
    enddo
    enddo
     
    do icol=1,levs
     vecm(icol)=dpkref(icol)/psref
    enddo
     
    do j=1,levs
     svhyb(j)=vecm(levs+1-j)
    do k=1,levs
      amhyb(k,j)=yecm(levs+1-k,levs+1-j)
      bmhyb(k,j)=tecm(levs+1-k,levs+1-j)
    enddo
    enddo
    
    do j=1,levs
    do k=1,levs
      amhyb(k,j)=amhyb(k,j)*beta/(rearth*rearth)
    enddo
    enddo
     
     
    return
end subroutine

!--------------------------------------------------------------------------------   
subroutine get_cd_hyb(dt,jcap1)
!--------------------------------------------------------------------------------   
    
    implicit none
    integer :: jcap1
    integer :: i,j,k,n,nn
    real :: dt,rnn1
    real :: ym(levs,levs)
    real :: rim(levs,levs)
    
    real :: ddd(jcap1),ppp(jcap1),rrr(jcap1)
    integer              lu(levs),mu(levs)
    real :: cons0,cons1  
    cons0 = 0.d0 
    cons1 = 1.d0
    
    do k=1,levs
    do j=1,levs
    rim(j,k)=cons0 
    enddo
    enddo
    
    do k=1,levs
    rim(k,k) = cons1     !constant
    enddo
    
    do i=1,levs
    
    do  j=1,levs
     ym(i,j) = tor_hyb(i)*svhyb(j)
    enddo
    
    do k=1,levs
    do j=1,levs
    ym(i,j) = ym(i,j) + amhyb(i,k)*bmhyb(k,j)
    enddo
    enddo
    
    enddo
    do nn=1,jcap1
    
     n = nn-1
     rnn1 =       n*(n+1)
    
     do i=1,levs
     do j=1,levs
      dm205_hyb(nn,i,j) = rim(i,j) + rnn1*dt*dt*ym(i,j)
     enddo
     enddo
    
    enddo
    !.............................................................
    call matinv(dm205_hyb,jcap1,levs,ddd,ppp,rrr)
    do nn=1,jcap1
    do i=1,levs
    do j=1,levs
    d_hyb_m(i,j,nn)=dm205_hyb(nn,i,j)
    enddo
    enddo
    enddo
    return
end subroutine

!--------------------------------------------------------------------------------   
subroutine matinv(a,m,n,d,p,r)
!--------------------------------------------------------------------------------   
      implicit none
      integer :: i,j,k,l,m,n
      real :: a(m,n,n),d(m),p(m),r(m)
      real :: cons0,cons1     !constant
      cons0 = 0.d0     !constant
      cons1 = 1.d0     !constant
      do l=1,m
      d(l)=cons1     !constant
      enddo
      do k=1,n
      do l=1,m
      p(l)=a(l,k,k)
      enddo
      do l=1,m
      r(l)=-cons1/p(l)     !constant
      enddo
      do l=1,m
      a(l,k,k)=cons0       !constant
      enddo
      do i=1,n
      do l=1,m
      a(l,i,k)=a(l,i,k)*r(l)
      enddo
      enddo
      do i=1,n
      if(i.eq.k) cycle
      do j=1,n
      do l=1,m
      a(l,i,j)=a(l,i,k)*a(l,k,j)+a(l,i,j)
      enddo
      enddo
      enddo 
      do l=1,m
      r(l)=-r(l)
      enddo
      do j=1,n
      do l=1,m
      a(l,k,j)=a(l,k,j)*r(l)
      enddo
      enddo
      do l=1,m
      d(l)=d(l)*p(l)
      enddo
      do l=1,m
      a(l,k,k)=r(l)
      enddo
      enddo
      return
end subroutine
end module implicit_mod
