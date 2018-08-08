module gwdrag_conv_mod

use gwdc_mod, only : gwdc

use constants_mod, only : GRAV, CP_AIR, RDGAS, RVGAS

implicit none
private

public :: gwdrag_conv

real, parameter :: FVIRT = RVGAS/RDGAS-1.

contains

!--------------------------------------------------------------------------------
subroutine gwdrag_conv(km, imax, u, v, t, q, plyr, delp, cuhr, plvl, ktop, &
                       kbot, kuo, fu, fv, dlen, tauctx, taucty)
!--------------------------------------------------------------------------------
    integer, intent(in) :: km,imax
    real, dimension(km,imax), intent(in) :: u, v, t, q, plyr, delp, cuhr
    real, dimension(km+1,imax), intent(in) :: plvl
    integer, dimension(imax), intent(in) :: ktop, kbot, kuo
    real, dimension(km,imax), intent(out) :: fu, fv
    real, dimension(imax), intent(in) :: dlen
    real, dimension(imax), intent(out) :: tauctx, taucty
    
    real, dimension(imax) :: qmax, cumabs
    real, dimension(km,imax) :: cumchr
    integer :: i, k, k1

    qmax(:)     = 0.0
    cumabs(:)   = 0.0
    cumchr(:,:) = 0.0

    do k = 1, km
        do i = 1, imax
            if ( k>=kbot(i) .and. k<=ktop(i) ) then
                qmax(i)     = max( qmax(i), cuhr(k,i) )
                cumabs(i)   = cuhr(k,i) + cumabs(i)
            endif
        enddo
    enddo

    do i = 1, imax
        do k = kbot(i), ktop(i)
            do k1 = kbot(i), k
                cumchr(k,i) = cuhr(k1,i) + cumchr(k,i)
            enddo
            cumchr(k,i) = cumchr(k,i) / cumabs(i)
         enddo
    enddo

    do i = 1, imax
        call gwdc(km, u(:,i), v(:,i), t(:,i), q(:,i), plyr(:,i), plvl(:,i), &
             delp(:,i), qmax(i), cumchr(:,i), ktop(i), kbot(i), kuo(i), &
             fu(:,i), fv(:,i), GRAV, CP_AIR, RDGAS, FVIRT, dlen(i), tauctx(i), &
             taucty(i))
    enddo

    return
end subroutine gwdrag_conv

end module gwdrag_conv_mod
