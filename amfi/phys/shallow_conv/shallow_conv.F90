module shallow_conv_mod

use shalcvt3_mod, only : shalcvt3

implicit none
private

public :: shallow_conv 

contains

!--------------------------------------------------------------------------------   
subroutine shallow_conv(imax, km, dt, del, prsi, prsl, prslk, kuo, q, t)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: imax, km
    real,    intent(in) :: dt
    real, dimension(km,imax), intent(in) :: del, prsl, prslk
    real, dimension(km+1,imax), intent(in) :: prsi
    integer, dimension(imax), intent(in) :: kuo
    real, dimension(km,imax), intent(inout) :: q, t

    real, dimension(imax) :: dpshc
    integer, dimension(imax) :: kmsc
    integer :: i, k

    dpshc = 0.3 * prsi(1,:)
    kmsc = 1 

    do i = 1, imax
        do k = 2, km
            if ( prsi(1,i)-prsi(k,i) <= dpshc(i) ) then
                kmsc(i) = k
            else
                exit
            endif
        enddo
    enddo 

    do i = 1, imax 
        call shalcvt3(kmsc(i), dt, del(:,i), prsi(:,i), prsl(:,i), prslk(:,i), &
                      kuo(i), q(:,i), t(:,i))
    enddo
    
    return
end subroutine shallow_conv

end module shallow_conv_mod 
