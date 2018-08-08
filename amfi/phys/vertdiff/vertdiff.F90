module vertdiff_mod

use moninp_mod, only : moninp

implicit none
private

public :: do_vertical_diffusion

real, allocatable, dimension(:,:,:) :: dqdt1, qgrs1

logical :: initialized = .false.

contains

subroutine init_vertical_diffusion(levs,ntrac,imax)
    integer :: levs, ntrac, imax
    
    if (initialized) return

    allocate(dqdt1(levs,ntrac,imax), qgrs1(levs,ntrac,imax))

    initialized = .true.

    return
end subroutine init_vertical_diffusion

!--------------------------------------------------------------------------------   
subroutine do_vertical_diffusion (imax,levs,ntrac,dvdt,dudt,dtdt,dqdt,ugrs,vgrs,tgrs,qgrs, &
                                  prsik,rb,ffmm,ffhh,qss,hflx,evap,stress,wind,kpbl, &
                                  prsi,del,prsl,prslk,phii,phil,dtp, &
                                  dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,xkzm)
!--------------------------------------------------------------------------------  
    implicit none
    integer, intent(in) :: imax, levs, ntrac 
    real, intent(in) :: dtp
    real, dimension(levs,imax,ntrac), intent(in) :: qgrs
    real, dimension(levs,imax),   intent(in) :: ugrs, vgrs, tgrs, del, prsl, prslk, phil
    real, dimension(levs+1,imax), intent(in) :: prsik, prsi, phii
    real, dimension(imax),        intent(in) :: rb, ffmm, ffhh, qss, hflx, evap, stress, wind

    real, dimension(levs,imax,ntrac), intent(out) :: dqdt
    real, dimension(levs,imax), intent(out) :: dvdt, dudt, dtdt
    real, dimension(imax), intent(out) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, hpbl, gamt, gamq, xkzm
    integer, dimension(imax), intent(out) :: kpbl

    integer :: i

    call init_vertical_diffusion(levs,ntrac,imax)        

    do i = 1, imax
        qgrs1(1:levs,1:ntrac,i) = qgrs(1:levs,i,1:ntrac)
    enddo

    dqdt1 = 0.

    do i = 1, imax
        call moninp(levs, ntrac, dvdt(:,i), dudt(:,i), dtdt(:,i), dqdt1(:,:,i), ugrs(:,i), vgrs(:,i), &
                    tgrs(:,i), qgrs1(:,:,i), prsik(1,i), rb(i), ffmm(i), ffhh(i), qss(i), &
                    hflx(i), evap(i), stress(i), wind(i), kpbl(i), prsi(:,i), del(:,i), &
                    prsl(:,i), prslk(:,i), phii(:,i), phil(:,i), dtp, dusfc1(i), &
                    dvsfc1(i), dtsfc1(i), dqsfc1(i), hpbl(i), gamt(i), gamq(i), xkzm(i))
    enddo

    do i = 1, imax
        dqdt(1:levs,i,1:ntrac) = dqdt1(1:levs,1:ntrac,i)
    enddo

    return
end subroutine do_vertical_diffusion

end module vertdiff_mod
