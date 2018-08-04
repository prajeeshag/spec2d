module vertdiff_mod

use moninp_mod, only : moninp

implicit none
private

public :: do_vertical_diffusion

contains


!--------------------------------------------------------------------------------   
subroutine do_vertical_diffusion (imax,levs,ntrac,dvdt,dudt,dtdt,dqdt,ugrs,vgrs,tgrs,qgrs, &
                                  prsik,rb,ffmm,ffhh,qss,hflx,evap,stress,wind,kpbl, &
                                  prsi,del,prsl,prslk,phii,phil,dtp, &
                                  dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,xkzm)
!--------------------------------------------------------------------------------  
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

    real, dimension(levs,ntrac) :: dqdt1
    real, dimension(levs,ntrac) :: qgrs1
    integer :: i                                                

    do i = 1, imax
        qgrs1 = qgrs(:,i,:)
        call moninp(levs, ntrac, dvdt(:,i), dudt(:,i), dtdt(:,i), dqdt1, ugrs(:,i), vgrs(:,i), &
                    tgrs(:,i), qgrs1, prsik(:,i), rb(i), ffmm(i), ffhh(i), qss(i), &
                    hflx(i), evap(i), stress(i), wind(i), kpbl(i), prsi(:,i), del(:,i), &
                    prsl(:,i), prslk(:,i), phii(:,i), phil(:,i), dtp, dusfc1(i), &
                    dvsfc1(i), dtsfc1(i), dqsfc1(i), hpbl(i), gamt(i), gamq(i), xkzm(i))
        dqdt(:,i,:) = dqdt1
    enddo

end subroutine do_vertical_diffusion

end module vertdiff_mod
