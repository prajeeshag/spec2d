module gfidi_mod

use, intrinsic :: iso_c_binding

use vertical_levels_mod, only: ak5=>ak, bk5=>bk, ck, dbk 
use constants_mod, only : rerth => RADIUS, OMEGA, rd=>RDGAS
use constants_mod, only : cp => CP_AIR, cvap => CP_VAPOR

implicit none
private

public :: gfidi_drv

integer, parameter :: lons_lat=1

contains

!--------------------------------------------------------------------------------   
subroutine gfidi_drv(levs, ntrac, ni, nj, deltim, sinlat_in, rcl_in, & 
           dg,tg,ug,vg,rqg,dphi,dlam,qg, &
           dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf, &
           dqdt,dtdt,drdt,dudt,dvdt,spdmax)
!--------------------------------------------------------------------------------   
    implicit none

    integer, intent(in) :: levs, ntrac, ni, nj
    real, intent(in) :: sinlat_in(nj), deltim, rcl_in(nj)
    real, intent(in)  :: dg(levs,nj*ni), tg(levs,nj*ni), &
                         ug(levs,nj*ni), vg(levs,nj*ni), &
                         rqg(levs,nj*ni,ntrac), &
                         dphi(nj*ni), dlam(nj*ni), qg(nj*ni)
    real, intent(in)  :: dtdf(levs,nj*ni), dtdl(levs,nj*ni), &
                         dudf(levs,nj*ni), dudl(levs,nj*ni), &
                         dvdf(levs,nj*ni), dvdl(levs,nj*ni), &
                         drdf(levs,nj*ni,ntrac), drdl(levs,nj*ni,ntrac)
    real, intent(out) :: dudt(levs,nj*ni), dvdt(levs,nj*ni), &
                         dqdt(nj*ni), dtdt(levs,nj*ni), &
                         drdt(levs,nj*ni,ntrac)
    real, intent(inout) :: spdmax(levs)

    real :: rcl(nj*ni), sinlat(nj*ni)
    type(c_ptr) :: ptr1, ptr2
    real, pointer :: ptmp1(:,:), ptmp2(:,:)

    real, dimension(levs,ntrac,nj*ni) :: rqg1, drdf1, drdl1, drdt1

    integer :: i, np

    np = nj*ni

    ptr1 = c_loc(rcl)
    call c_f_pointer(ptr1,ptmp1,[nj,ni])

    ptr2 = c_loc(sinlat)
    call c_f_pointer(ptr2,ptmp2,[nj,ni])

    do i = 1, ni
        ptmp1(1:nj,i) = rcl_in(1:nj) 
        ptmp2(1:nj,i) = sinlat_in(1:nj)
    enddo

    do i = 1, ntrac
        rqg1(:,i,:) = rqg(:,:,i) 
        drdf1(:,i,:) = drdf(:,:,i) 
        drdl1(:,i,:) = drdl(:,:,i) 
    enddo

    do i = 1, np
        call gfidi_hyb(levs, ntrac, deltim, sinlat(i), rcl(i), & 
           dg(1:levs,i),tg(1:levs,i),ug(1:levs,i),vg(1:levs,i), &
           rqg1(1:levs,1:ntrac,i),dphi(i),dlam(i),qg(i), &
           dtdf(1:levs,i),dtdl(1:levs,i),drdf1(1:levs,1:ntrac,i), &
           drdl1(1:levs,1:ntrac,i),dudl(1:levs,i),dvdl(1:levs,i), &
           dudf(1:levs,i),dvdf(1:levs,i), &
           !output
           dqdt(i),dtdt(1:levs,i),drdt1(1:levs,1:ntrac,i), &
           dudt(1:levs,i),dvdt(1:levs,i),spdmax(1:levs))
    enddo

    do i = 1, ntrac
        drdt(:,:,i) = drdt1(:,i,:)
    enddo 
    
    return
end subroutine gfidi_drv

!--------------------------------------------------------------------------------   
subroutine gfidi_hyb(levs, ntrac, deltim, sinlat, rcl, &
           dg,tg,ug,vg,rqg,dphi,dlam,qg, &
           dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf, &
           dqdt,dtdt,drdt,dudt,dvdt,spdmax)
!--------------------------------------------------------------------------------   
 
! sela  add calculation of spdlat as in gfidiu
! Fanglin Yang, June 2007, use flux-limited scheme for vertical advection of tracers
    implicit none

    integer, intent(in) :: levs, ntrac
    real, intent(in) :: sinlat, deltim, rcl
    real, intent(in)  :: dg(lons_lat,levs), tg(lons_lat,levs), &
                         ug(lons_lat,levs), vg(lons_lat,levs), &
                         rqg(lons_lat,levs,ntrac), &
                         dphi(lons_lat), dlam(lons_lat), qg(lons_lat)
    real, intent(in)  :: dtdf(lons_lat,levs), dtdl(lons_lat,levs), &
                         dudf(lons_lat,levs), dudl(lons_lat,levs), &
                         dvdf(lons_lat,levs), dvdl(lons_lat,levs), &
                         drdf(lons_lat,levs,ntrac), drdl(lons_lat,levs,ntrac)
    real, intent(out) :: dudt(lons_lat,levs), dvdt(lons_lat,levs), &
                         dqdt(lons_lat), dtdt(lons_lat,levs), &
                         drdt(lons_lat,levs,ntrac)
    real, intent(inout) :: spdmax(levs)

    !local 
    real :: coriol,rk,sinra
    real :: pk5(lons_lat,levs+1), dpk(lons_lat,levs)
    real :: zadv(lons_lat,levs,3+ntrac)
    real :: rqg_half(lons_lat,0:levs,ntrac), rqg_d(lons_lat,0:levs,ntrac)
 
    real :: dot(lons_lat,levs+1), dotinv(lons_lat,levs+1), &
            ek(lons_lat,levs),      cg(lons_lat,levs), &
            cb(lons_lat,levs),      db(lons_lat,levs), &
            zlam(lons_lat,levs),    zphi(lons_lat,levs), &
            worka(lons_lat,levs),   workb(lons_lat,levs), &
            workc(lons_lat,levs), &
            phiu(lons_lat,levs),    phiv(lons_lat,levs), &
            uprs(lons_lat,levs),    vprs(lons_lat,levs), &
            cofa(lons_lat,levs),    cofb(lons_lat,levs), &
            alfa(lons_lat,levs),    rlnp(lons_lat,levs), &
            px1u(lons_lat,levs),    px1v(lons_lat,levs), &
            px2u(lons_lat,levs),    px2v(lons_lat,levs), &
            px2(lons_lat,levs), &
            px3u(lons_lat,levs),    px3v(lons_lat,levs), &
            px4u(lons_lat,levs),    px4v(lons_lat,levs), &
            px5u(lons_lat,levs),    px5v(lons_lat,levs), &
            uphi(lons_lat,levs),    vphi(lons_lat,levs), &
            expq(lons_lat), &
            rdel(lons_lat,levs),   rdel2(lons_lat,levs), &
            sumdel(lons_lat),          del(lons_lat,levs), &
            si(lons_lat,levs+1),     sl(lons_lat,levs), &
            rk1,rkr,                 dif(lons_lat)
 
    real :: cons0,cons0p5,cons1,cons2,clog2   !constant
    real :: rmin,rmax,delta,delta1
    real :: rrkp,rrk1m,phkp,phk1m,bb,cc,tmpdrdt
    real :: xvcn
    integer :: nvcn
    integer :: j,k,n
    integer :: levp1, levm1

    levp1 = levs + 1
    levm1 = levs - 1
 
    cons0   = 0.d0  
    cons0p5 = 0.5d0 
    cons1   = 1.d0  
    cons2   = 2.d0 
    rk= rd /cp

!-------------------------------------------------------
   
    sinra=sinlat
    coriol=cons2*omega*sinra
    sinra=sinra/rerth

!-------------------------------------------------------
 
    clog2=log(cons2)     ! constant
    delta=cvap/cp  ! check these cpv cpd (at const p for vapor and dry
    delta1=delta-cons1
    rk1 = rk + 1.e0
    rkr = 1.0/rk
 
 
    do j=1,lons_lat
        expq(j)=exp(qg(j))
    enddo
!.............................................................
! get vertical coordinate for vcnhyb going bot. to top.
    do k=1,levp1
        do j=1,lons_lat
            si(j,levs+2-k)= ak5(k)+bk5(k)*expq(j) !ak(k) bk(k) go top to bottom
        enddo
    enddo
 
 
    do  k=1,levs
        do j=1,lons_lat
            sl(j,k) = cons0p5*(si(j,k)+si(j,k+1))
        enddo
    enddo
!.............................................................
 
    do k=1,levp1
        do j=1,lons_lat
            pk5(j,k)=ak5(k) + bk5(k)*expq(j)
        enddo
    enddo
 
    do k=1,levs
        do j=1,lons_lat
            dpk(j,k)=    pk5(j,k+1) - pk5(j,k)
            rdel(j,k)=    cons1/dpk(j,k)            ! constant
            rdel2(j,k)=  cons0p5/dpk(j,k)            ! constant
        enddo
    enddo
 
    k=1
    do j=1,lons_lat
        alfa(j,1)=clog2                          ! constant
    enddo
 
    do j=1,lons_lat
        rlnp(j,1)= 99999.99
    enddo

    do  k=2,levs
        do j=1,lons_lat
            rlnp(j,k)= log( pk5(j,k+1)/pk5(j,k) )
            alfa(j,k)= cons1-( pk5(j,k)/dpk(j,k) )*rlnp(j,k)
        enddo
    enddo
 
    spdmax=0.
    do  k=1,levs
        do j=1,lons_lat
            ek(j,k)=(ug(j,levs+1-k)*ug(j,levs+1-k) &
                   + vg(j,levs+1-k)*vg(j,levs+1-k))*rcl
            if (ek(j,k) .gt. spdmax(levs+1-k))  spdmax(levs+1-k)=ek(j,k)
        enddo
    enddo

    do k=1,levs
        do j=1,lons_lat
            cg(j,k)=(ug(j,levs+1-k)*dlam(j)+vg(j,levs+1-k)*dphi(j))*rcl
        enddo
    enddo

    k=1
    do j=1,lons_lat
        db(j,1)=dg(j,levs)*dpk(j,1)
        cb(j,1)=cg(j,1)*dbk(1)
    enddo
 
    do k=1,levm1
        do j=1,lons_lat
            db(j,k+1)=db(j,k)+dg(j,levs-k)*dpk(j,k+1)
            cb(j,k+1)=cb(j,k)+cg(j,k+1)*dbk(k+1)
        enddo
    enddo
 
    do j=1,lons_lat
        dqdt(j)= -db(j,levs)/expq(j)-cb(j,levs)
        dot(j,    1)=cons0                    !constant
        dot(j,levp1)=cons0                    !constant
    enddo

    do k=1,levs-1
        do j=1,lons_lat
            dot(j,k+1)=-expq(j)*(bk5(k+1)*dqdt(j)+cb(j,k)) -db(j,k)
        enddo
    enddo

    do k=1,levp1
        do j=1,lons_lat
            dotinv(j,k)=dot(j,levp1+1-k)
        enddo
    enddo
 
! variables are in bottom to top order  !!!!!!!!!!!!!!!!!
! do horizontal advection.
 
    k=1
    do j=1,lons_lat
        dudt(j,levs+1-k)=-ug(j,levs+1-k)*dudl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dudf(j,levs+1-k)
        dvdt(j,levs+1-k)=-ug(j,levs+1-k)*dvdl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dvdf(j,levs+1-k)
        dtdt(j,levs+1-k)=-ug(j,levs+1-k)*dtdl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dtdf(j,levs+1-k)
    enddo
 
    k=levs
    do j=1,lons_lat
        dudt(j,levs+1-k)=-ug(j,levs+1-k)*dudl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dudf(j,levs+1-k)
        dvdt(j,levs+1-k)=-ug(j,levs+1-k)*dvdl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dvdf(j,levs+1-k)
        dtdt(j,levs+1-k)=-ug(j,levs+1-k)*dtdl(j,levs+1-k) &
                         -vg(j,levs+1-k)*dtdf(j,levs+1-k)
    enddo
 
    do k=2,levm1
        do j=1,lons_lat
            dudt(j,levs+1-k)=-ug(j,levs+1-k)*dudl(j,levs+1-k) &
                             -vg(j,levs+1-k)*dudf(j,levs+1-k)
            dvdt(j,levs+1-k)=-ug(j,levs+1-k)*dvdl(j,levs+1-k) &
                             -vg(j,levs+1-k)*dvdf(j,levs+1-k)
            dtdt(j,levs+1-k)=-ug(j,levs+1-k)*dtdl(j,levs+1-k) &
                             -vg(j,levs+1-k)*dtdf(j,levs+1-k)
        enddo
    enddo
! add coriolis,deformation note  coriolis sign for s.hemi
    do k=1,levs
        do j=1,lons_lat
            dudt(j,levs+1-k)=dudt(j,levs+1-k)+vg(j,levs+1-k)*coriol
            dvdt(j,levs+1-k)=dvdt(j,levs+1-k)-ug(j,levs+1-k)*coriol &
                            -sinra*ek(j,k)
        enddo
    enddo
!.................................................................
! calculate pressure force:
    k=1
    do j=1,lons_lat
       cofb(j,k)=-rdel(j,k)*(                 alfa(j,k)*dbk(k))
    enddo
 
    do k=2,levs
        do j=1,lons_lat
            cofb(j,k)=-rdel(j,k)*(bk5(k)*rlnp(j,k)+alfa(j,k)*dbk(k))
        enddo
    enddo
 
    do k=1,levs
        do j=1,lons_lat
            uprs(j,k)=cofb(j,k)*rd*tg(j,levs+1-k)*expq(j)*dlam(j)
            vprs(j,k)=cofb(j,k)*rd*tg(j,levs+1-k)*expq(j)*dphi(j)
        enddo
    enddo

    do k=1,levs
        do j=1,lons_lat
            cofa(j,k)=-rdel(j,k)*(  &
                bk5(k+1)*pk5(j,k)/pk5(j,k+1) - bk5(k) &
                +rlnp(j,k)*( bk5(k)-pk5(j,k)*dbk(k)*rdel(j,k)))
        enddo
    enddo
!.................................................................
    do k=1,levs
        do j=1,lons_lat
            px1u(j,k)=cons0              ! grid topography =0 for testing
            px1v(j,k)=cons0              ! grid tpopgraphy =0 for testing
        enddo
    enddo
! see programming notes for looping in  calculating  px2 and px3
    do j=1,lons_lat
        px2(j,levs)=cons0                             ! constant
        px2(j,levs-1)= -rd*( bk5(levp1)/pk5(j,levp1)-bk5(levs)/pk5(j,levs))*tg(j,1)
    enddo
 
 
    do k=2,levs-1
        do j=1,lons_lat
            px2(j,levs-k)=px2(j,levp1-k) &
            -rd*(bk5(levs+2-k)/pk5(j,levs+2-k)-bk5(levp1-k)/pk5(j,levp1-k))* &
             tg(j,k)
        enddo
    enddo
 
    do k=1,levs
        do j=1,lons_lat
            px2u(j,k)=px2(j,k)*expq(j)*dlam(j)
            px2v(j,k)=px2(j,k)*expq(j)*dphi(j)
        enddo
    enddo

    do j=1,lons_lat
       px3u(j,levs)=cons0 ! constant
       px3v(j,levs)=cons0 ! constant
    enddo
 
    do j=1,lons_lat
       px3u(j,levs-1)=-rd*rlnp(j,levs)*dtdl(j,1)
       px3v(j,levs-1)=-rd*rlnp(j,levs)*dtdf(j,1)
    enddo
 
    do k=2,levs-1
        do j=1,lons_lat
            px3u(j,levs-k)=px3u(j,levp1-k)-rd*rlnp(j,levp1-k)*dtdl(j,k)
            px3v(j,levs-k)=px3v(j,levp1-k)-rd*rlnp(j,levp1-k)*dtdf(j,k)
        enddo
    enddo

    do k=1,levs
        do j=1,lons_lat
            px3u(j,k)=px3u(j,k)/rcl
            px3v(j,k)=px3v(j,k)/rcl
        enddo
    enddo
    do k=1,levs
        do j=1,lons_lat
            px4u(j,k)=-rd*alfa(j,k)*dtdl(j,levs+1-k)/rcl
            px4v(j,k)=-rd*alfa(j,k)*dtdf(j,levs+1-k)/rcl
        enddo
    enddo
 
    do k=1,levs
    do j=1,lons_lat
        px5u(j,k)=-cofa(j,k)*rd*tg(j,levs+1-k)*expq(j)*dlam(j)
        px5v(j,k)=-cofa(j,k)*rd*tg(j,levs+1-k)*expq(j)*dphi(j)
    enddo
    enddo
 
    do k=1,levs
    do j=1,lons_lat
        uphi(j,k)=px1u(j,k)+px2u(j,k)+px3u(j,k)+px4u(j,k)+px5u(j,k)
        vphi(j,k)=px1v(j,k)+px2v(j,k)+px3v(j,k)+px4v(j,k)+px5v(j,k)
    enddo
    enddo
    do k=1,levs
    do j=1,lons_lat
        dudt(j,levs+1-k)=dudt(j,levs+1-k)+uphi(j,k)+uprs(j,k)
        dvdt(j,levs+1-k)=dvdt(j,levs+1-k)+vphi(j,k)+vprs(j,k)
    enddo
    enddo
!............................................................
    do k=1,levs
    do j=1,lons_lat
        worka(j,k)=rk*tg(j,levs+1-k)/(cons1+(delta1)* &
        rqg(j,levs+1-k,1)) * rdel(j,k)
    enddo
    enddo
 
    k=1
    do j=1,lons_lat
       workb(j,1)= alfa(j,1)*( dg(j,levs)*dpk(j,1)+expq(j)*cb(j,1)*dbk(1) )
    enddo
 
    do k=2,levs
    do j=1,lons_lat
        workb(j,k)=rlnp(j,k)*( db(j,k-1)+expq(j)*cb(j,k-1) ) &
        +alfa(j,k)*( dg(j,levs+1-k)*dpk(j,k)+expq(j)*cg(j,k)*dbk(k) )
    enddo
    enddo
 
    k=1
    do j=1,lons_lat
       workc(j,1)=expq(j)*cg(j,1)*dbk(1)
    enddo
 
    do k=2,levs
    do j=1,lons_lat
        workc(j,k)=expq(j)*cg(j,k)*( dbk(k)+ck(k)*rlnp(j,k)*rdel(j,k) )
    enddo
    enddo
 
    do k=1,levs
    do j=1,lons_lat
       dtdt(j,levs+1-k)= dtdt(j,levs+1-k)+worka(j,k)*( -workb(j,k) + workc(j,k))
    enddo
    enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    do n=1,ntrac
        k=1
        do j=1,lons_lat
            drdt(j,levs+1-k,n)=-ug(j,levs+1-k)*drdl(j,levs+1-k,n) &
                          -vg(j,levs+1-k)*drdf(j,levs+1-k,n)
        enddo
        k=levs
        do j=1,lons_lat
            drdt(j,levs+1-k,n)=-ug(j,levs+1-k)*drdl(j,levs+1-k,n) &
                          -vg(j,levs+1-k)*drdf(j,levs+1-k,n)
        enddo
        do k=2,levm1
            do j=1,lons_lat
                drdt(j,levs+1-k,n)=-ug(j,levs+1-k)*drdl(j,levs+1-k,n) &
                          -vg(j,levs+1-k)*drdf(j,levs+1-k,n)
            enddo
        enddo
    enddo 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
! do vertical advection
 
      k=1
      do j=1,lons_lat
        zadv(j,levs+1-k,1)= -rdel2(j,k)*dot(j,k+1)*( ug(j,levs-k)-ug(j,levs+1-k))
 
        zadv(j,levs+1-k,2)= -rdel2(j,k)*dot(j,k+1)*( vg(j,levs-k)-vg(j,levs+1-k))
 
        zadv(j,levs+1-k,3)= -rdel2(j,k)*dot(j,k+1)*( tg(j,levs-k)-tg(j,levs+1-k))
      enddo
 
      k=levs
      do j=1,lons_lat
        zadv(j,levs+1-k,1)= -rdel2(j,k)*dot(j,k)*( ug(j,levs+1-k)-ug(j,levs+2-k) )
 
        zadv(j,levs+1-k,2)= -rdel2(j,k)*dot(j,k)*( vg(j,levs+1-k)-vg(j,levs+2-k) )
 
        zadv(j,levs+1-k,3)= -rdel2(j,k)*dot(j,k)*( tg(j,levs+1-k)-tg(j,levs+2-k) )
      enddo
 
      do k=2,levm1
      do j=1,lons_lat
        zadv(j,levs+1-k,1)= &
         -rdel2(j,k)*( dot(j,k+1)*( ug(j,levs  -k)-ug(j,levs+1-k) ) + &
          dot(j,k  )*( ug(j,levs+1-k)-ug(j,levs+2-k) ) )
 
        zadv(j,levs+1-k,2)= &
         -rdel2(j,k)*( dot(j,k+1)*( vg(j,levs  -k)-vg(j,levs+1-k) ) + &
          dot(j,k  )*( vg(j,levs+1-k)-vg(j,levs+2-k) ) )
 
        zadv(j,levs+1-k,3)= &
         -rdel2(j,k)*( dot(j,k+1)*( tg(j,levs  -k)-tg(j,levs+1-k) ) + &
          dot(j,k  )*( tg(j,levs+1-k)-tg(j,levs+2-k) ) )
 
      enddo
      enddo

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Fanglin Yang, June 2007
! 1. use Total Variation Diminishing (TVD) flux-limited scheme
!    for vertical advection of tracers ( J. Thuburn, QJRMS, 1993, 469-487)
!    Vertical advection  dQ/dt = AA = -W*dQ/dP = -[d(Q*W)/dP-Q*dW/dP]
!    let BB=d(Q*W)/dP and CC=-Q*dW/dP, AA=-(BB+CC), then Q(n+1)=Q(n)+AA*dt
! 2. The current scheme is central in space and central in time.  To use
!    the TVD scheme for vertical advection, the time differencing must be
!    forward in time otherwise it is unstable.  However, for horizonatl
!    advection which is center in space, the forward-in-time scheme is
!    always unstable.  To overcome this conflict, the vertical adevtion 
!    from time step n-1 is used to get mean advection at current time
!    step.  Then, the central-in-time scheme is applied to both the 
!    vertical and horizontal advections of tracers.
!
     do n=1,ntrac
         k=1
         do j=1,lons_lat
             zadv(j,levs+1-k,3+n)= -rdel2(j,k)*dot(j,k+1)*( rqg(j,levs-k,n)-rqg(j,levs+1-k,n) )
         enddo
 
         k=levs
         do j=1,lons_lat
             zadv(j,levs+1-k,3+n)= -rdel2(j,k)*dot(j,k)*( rqg(j,levs+1-k,n)-rqg(j,levs+2-k,n) )
         enddo
 
         do k=2,levm1
             do j=1,lons_lat
                 zadv(j,levs+1-k,3+n)=-rdel2(j,k)* &
                  ( dot(j,k+1)*( rqg(j,levs  -k,n)-rqg(j,levs+1-k,n) ) + &
                  dot(j,k  )*( rqg(j,levs+1-k,n)-rqg(j,levs+2-k,n) ) )
             enddo
         enddo
     enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call vcnhyb(lons_lat,levs,3+ntrac,deltim, &
                 si,sl,dotinv,zadv,nvcn,xvcn)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
! add vertical filterd advection
     do k=1,levs
     do j=1,lons_lat
      dudt(j,k)=dudt(j,k)+zadv(j,k,1)
      dvdt(j,k)=dvdt(j,k)+zadv(j,k,2)
      dtdt(j,k)=dtdt(j,k)+zadv(j,k,3)
     enddo
     enddo
     do  n=1,ntrac
      do k=1,levs
      do j=1,lons_lat
       drdt(j,k,n)=drdt(j,k,n)+zadv(j,k,3+n)
      enddo
      enddo
     enddo
 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! this multiplication must be on  completed tendencies.
     do k=1,levs
     do j=1,lons_lat
       dudt(j,levs+1-k)=dudt(j,levs+1-k)*rcl
       dvdt(j,levs+1-k)=dvdt(j,levs+1-k)*rcl
     enddo
     enddo

!     print *,' leave gfidi_hyb_fd ' 		! hmhj

     return
end subroutine gfidi_hyb



!--------------------------------------------------------------------------------    
subroutine vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!--------------------------------------------------------------------------------   
!c                .      .    .                                       .
!c subprogram:    vcnhyb      vertical advection instability filter
!c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
!c
!c abstract: filters vertical advection tendencies
!c   in the dynamics tendency equation in order to ensure stability
!c   when the vertical velocity exceeds the cfl criterion.
!c   the vertical velocity in this case is sigmadot.
!c   for simple second-order centered eulerian advection,
!c   filtering is needed when vcn=zdot*dt/dz>1.
!c   the maximum eigenvalue of the linear advection equation
!c   with second-order implicit filtering on the tendencies
!c   is less than one for all resolvable wavenumbers (i.e. stable)
!c   if the nondimensional filter parameter is nu=(vcn**2-1)/4.
!c
!c program history log:
!c   97-07-30  iredell
!c
!c usage:    call vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!c
!c   input argument list:
!c     im       - integer number of gridpoints to filter
!c     km       - integer number of vertical levels
!c     nm       - integer number of fields
!c     dt       - real timestep in seconds
!c     zint     - real (im,km+1) interface vertical coordinate values
!c     zmid     - real (im,km) midlayer vertical coordinate values
!c     zdot     - real (im,km+1) vertical coordinate velocity
!c     zadv     - real (im,km,nm) vertical advection tendencies
!c
!c   output argument list:
!c     zadv     - real (im,km,nm) vertical advection tendencies
!c     nvcn     - integer number of points requiring filtering
!c     xvcn     - real maximum vertical courant number
!c
!c   subprograms called:
!c     tridim_hyb   - tridiagonal matrix solver
!c
    implicit none
    integer,intent(in):: im,km,nm
    real,intent(in):: dt,zint(im,km+1),zmid(im,km),zdot(im,km+1)
    real,intent(inout):: zadv(im,km,nm)
    integer,intent(out):: nvcn
    real,intent(out):: xvcn
    integer :: i,j,k,n,ivcn(im)
    logical :: lvcn(im)
    real :: zdm,zda,zdb,vcn(im,km-1)
    real :: rnu,cm(im,km),cu(im,km-1),cl(im,km-1)
    real :: rr(im,km,nm)
!c  compute vertical courant number
!c  increase by 10% for safety
      nvcn=0
      xvcn=0.
      lvcn=.false.
      do k=1,km-1
        do i=1,im
          zdm=zmid(i,k)-zmid(i,k+1)
          vcn(i,k)=abs(zdot(i,k+1)*dt/zdm)*1.1
          lvcn(i)=lvcn(i).or.vcn(i,k).gt.1
          xvcn=max(xvcn,vcn(i,k))
        enddo
      enddo
!c  determine points requiring filtering
      if(xvcn.gt.1) then
        do i=1,im
          if(lvcn(i)) then
            ivcn(nvcn+1)=i
            nvcn=nvcn+1
          endif
        enddo
!c  compute tridiagonal matrim
        do j=1,nvcn
          cm(j,1)=1
        enddo
        do k=1,km-1
          do j=1,nvcn
            i=ivcn(j)
            if(vcn(i,k).gt.1) then
             zdm=zmid(i,k)-zmid(i,k+1)
             zda=zint(i,k+1)-zint(i,k+2)
             zdb=zint(i,k)-zint(i,k+1)
              rnu=(vcn(i,k)**2-1)/4
              cu(j,k)=-rnu*zdm/zdb
              cl(j,k)=-rnu*zdm/zda
              cm(j,k)=cm(j,k)-cu(j,k)
              cm(j,k+1)=1-cl(j,k)
            else
              cu(j,k)=0
              cl(j,k)=0
              cm(j,k+1)=1
            endif
          enddo
        enddo
!c  fill fields to be filtered
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              rr(j,k,n)=zadv(i,k,n)
            enddo
          enddo
        enddo
!c  solve tridiagonal system
        call tridim_hyb(nvcn,im,km,km,nm,cl,cm,cu,rr,cu,rr)
!c  replace filtered fields
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              zadv(i,k,n)=rr(j,k,n)
            enddo
          enddo
        enddo
      endif
end subroutine vcnhyb

!--------------------------------------------------------------------------------   
subroutine tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!--------------------------------------------------------------------------------   
!c                .      .    .                                       .
!c subprogram:    tridim_hyb      solves tridiagonal matrix problems.
!c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
!c
!c abstract: this routine solves multiple tridiagonal matrix problems
!c   with multiple right-hand-side and solution vectors for every matrix.
!c   the solutions are found by eliminating off-diagonal coefficients,
!c   marching first foreward then backward along the matrix diagonal.
!c   the computations are vectorized around the number of matrices.
!c   no checks are made for zeroes on the diagonal or singularity.
!c
!c program history log:
!c   97-07-30  iredell
!c
!c usage:    call tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!c
!c   input argument list:
!c     l        - integer number of tridiagonal matrices
!c     lx       - integer first dimension (lx>=l)
!c     n        - integer order of the matrices
!c     nx       - integer second dimension (nx>=n)
!c     m        - integer number of vectors for every matrix
!c     cl       - real (lx,2:n) lower diagonal matrix elements
!c     cm       - real (lx,n) main diagonal matrix elements
!c     cu       - real (lx,n-1) upper diagonal matrix elements
!c                (may be equivalent to au if no longer needed)
!c     r        - real (lx,nx,m) right-hand-side vector elements
!c                (may be equivalent to a if no longer needed)
!c
!c   output argument list:
!c     au       - real (lx,n-1) work array
!c     a        - real (lx,nx,m) solution vector elements
!c
!c attributes:
!c   language: fortran 77.
!c   machine:  cray.
!c
    implicit none
    integer, intent(in) :: l, n, lx, nx, m
    real :: cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m), &
            au(lx,n-1),a(lx,nx,m)
    integer :: i, j, k
    real :: fk

!c  march up
      do i=1,l
        fk=1./cm(i,1)
        au(i,1)=fk*cu(i,1)
      enddo
      do j=1,m
        do i=1,l
          fk=1./cm(i,1)
          a(i,1,j)=fk*r(i,1,j)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)=fk*cu(i,k)
        enddo
        do j=1,m
          do i=1,l
            fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
            a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
          enddo
        enddo
      enddo
!c  march down
      do j=1,m
        do i=1,l
          fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
          a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
        enddo
      enddo
      do k=n-1,1,-1
        do j=1,m
          do i=1,l
            a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
          enddo
        enddo
      enddo
end subroutine tridim_hyb

end module gfidi_mod
