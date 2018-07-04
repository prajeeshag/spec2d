module horiz_diffusion_mod

use constants_mod, only : rerth => RADIUS,  rd => RDGAS, cp => CP_AIR 

use mpp_mod, only : mpp_broadcast, mpp_pe, mpp_root_pe, mpp_gather, mpp_npes
use mpp_mod, only : mpp_error, FATAL
    
implicit none
private

public :: init_horiz_diffusion, horiz_diffusion


real, parameter :: rkappa = cp / rd

integer :: levs
real, allocatable :: rfact(:,:,:), rfactrd(:,:,:)
real, allocatable :: bkly(:), sf(:), rrfact(:,:,:)

logical :: coef_pe = .false.

integer :: coef_idx = -1, coef_from_pe = -1

contains

!--------------------------------------------------------------------------------   
subroutine init_horiz_diffusion(jcap,deltim,sl,sph_wave,bk5)
!--------------------------------------------------------------------------------   
    implicit none
    integer, intent(in) :: jcap, sph_wave(:,:)
    real, intent(in) :: sl(:), deltim, bk5(:)

    real :: rsl(size(sl,1)), rthk(size(sl,1)), rtrd(size(sl,1))
    real :: dn(size(sph_wave,1),2), realval, fshk, rtnp, slrd0
    integer, allocatable :: coefpelist(:)
    real :: rtrd1, dn1, fact
    integer :: n0, jdel, np, jdelh, npd, nd
    integer :: i, k, j, kd, ku

    levs = size(sl,1)

    do i = 1, size(sph_wave,1)
        if (sph_wave(i,1) == 0) then
            coef_pe = .true.
            coef_idx = i
            coef_from_pe = mpp_pe()
            exit
        endif
    enddo

    allocate(coefpelist(mpp_npes()))

    call mpp_gather([coef_from_pe],coefpelist)

    if (mpp_pe()==mpp_root_pe()) then
        if (count(coefpelist>=0)/=1) &
            call mpp_error('init_horiz_diffusion','Error in finding coef_from_pe', FATAL)

        coef_from_pe=maxval(coefpelist)
    endif

    call mpp_broadcast(coef_from_pe,mpp_root_pe())

    allocate(bkly(levs),sf(levs)) 
    allocate(rfact(levs,size(sph_wave,1),2))
    allocate(rrfact(levs,size(sph_wave,1),2))
    allocate(rfactrd(levs,size(sph_wave,1),2))

    rsl(:) = 0.
    where(sl/=0.) rsl = 1./sl

    bkly(:) = 1.0
    do  k=1,levs
        bkly(k)=0.5*(bk5(levs-k+1)+bk5(levs-k+2))*rsl(k)
    enddo

    n0=0             ! maximum wavenumber for zero diffusion
    jdel=8           ! order of diffusion (even power to raise del)
    np=jcap
    fshk=1.0         ! extra height-dependent diffusion factor per scale height
    if(jcap.gt.170) then
      rtnp=(jcap/170.)**4*1.1/3600 !reciprocal of time scale of diffusion at reference wavenumber np
      fshk=2.2 ! extra height-dependent diffusion factor per scale height
    elseif(jcap.eq.170) then
      rtnp=4*3.e15/(rerth**4)*float(80*81)**2
    elseif(jcap.eq.126) then 
      fshk=2.2 ! extra height-dependent diffusion factor per scale height
      rtnp=4*3.e15/(rerth**4)*float(80*81)**2 
    else
      rtnp=1*3.e15/(rerth**4)*float(80*81)**2
    endif

    slrd0=0.002        ! sigma level at which to begin rayleigh damping
    rtrd1=1./(5*86400) ! reciprocal of time scale per scale height
                       !  above beginning sigma level for rayleigh damping
    do k=1,levs
      rtrd(k) = 0.
      if(sl(k).lt.slrd0) then
        if (rsl(k)/=0.) rtrd(k)=rtrd1*log(slrd0*rsl(k))
      endif
      rthk(k)=(sl(k))**log(1/fshk)
    enddo

    jdelh=jdel/2
    npd=max(np-n0,0)
    realval=npd*(npd+1)
    dn1=2.*rtnp/realval**jdelh

    do i = 1, size(sph_wave,1)
        nd=max(sph_wave(i,1)-n0,0)
        realval=nd*(nd+1)
        dn(i,1)=dn1*realval**jdelh
    enddo

    do i = 1, size(sph_wave,1)
        nd=max(sph_wave(i,2)-n0,0)
        realval=nd*(nd+1)
        dn(i,2)=dn1*realval**jdelh
    enddo

    do j = 1, 2
        do i = 1, size(sph_wave,1)
            do k = 1, levs
              fact = deltim*dn(i,j)*rthk(k)
              rfact(k,i,j) = 1./(1.+fact)
              rfactrd(k,i,j) = 1./(1.+fact+deltim*rtrd(k))
              rrfact(k,i,j) = fact*rfact(k,i,j)
            enddo
        enddo
    enddo

    do k=1,levs
        kd=max(k-1,1)
        ku=min(k+1,levs)
        sf(k)=sl(k)/(sl(ku)-sl(kd))/sqrt(2.)
    enddo

    return

end subroutine init_horiz_diffusion


!--------------------------------------------------------------------------------   
subroutine horiz_diffusion(rte,we,xe,ye,qme)
!--------------------------------------------------------------------------------   
    implicit none
    
    complex, intent(inout) :: rte(:,:,:,:)
    complex, intent(inout) :: we(:,:,:)
    complex, intent(inout) :: xe(:,:,:)
    complex, intent(inout) :: ye(:,:,:)
    complex, intent(in) :: qme(:,:)

    real :: coef00(size(ye,1))
    integer :: i, j, k, n

    coef00(:) = 0.
    if (coef_pe) then
        coef00(:) = real(ye(:,coef_idx,1))
    endif 

    call mpp_broadcast(coef00,size(coef00),coef_from_pe)
    call updown(coef00)

    we(:,:,:) = we(:,:,:)*rfactrd(:,:,:)
    xe(:,:,:) = xe(:,:,:)*rfactrd(:,:,:)

    do j = 1, size(ye,3) 
        do i = 1, size(ye,2) 
            ye(:,i,j) = ye(:,i,j)*rfact(:,i,j) + rrfact(:,i,j) * coef00(:) * bkly(:) * qme(i,j)
        enddo
    enddo

    do n = 1, size(rte,4)
        rte(:,:,:,n) = rte(:,:,:,n)*rfact(:,:,:)
    enddo
    
    return
end subroutine horiz_diffusion

!--------------------------------------------------------------------------------   
subroutine updown(zoncoef)
!--------------------------------------------------------------------------------   
   implicit none
   real, intent(inout) :: zoncoef(:)
   real :: xd(size(zoncoef,1)), xu(size(zoncoef,1))
   integer :: k

   xu(levs)=zoncoef(levs)
   do k=1,levs-1
     xu(k)=zoncoef(k+1)
   enddo

   xd(1)=zoncoef(1)
   do k=2,levs
     xd(k)=zoncoef(k-1)
   enddo

    do k=1,levs
      zoncoef(k)=(xu(k)-xd(k))*sf(k)
    enddo

    return
end subroutine updown

end module horiz_diffusion_mod
