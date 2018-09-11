module aqua_planet_mod

use mpp_mod, only : mpp_error, FATAL, NOTE
use constants_mod, only : PI, KELVIN

implicit none

private

public :: aquape_init_temp, get_aquape_sst, aquape_o3

include "o3_dat"

include "temp_dat"

contains

subroutine aquape_init_temp(latin, prsl, ta)
    real, intent(in) :: latin, prsl(:) ! pressure levels in hPa(millibar)
    real, intent(out) :: ta(:)
    
    integer :: j1, j2, i
    real :: ta1(nlev_temp), w1, w2, w12

    real :: inlat(1)

    if (abs(latin)>90.) call mpp_error(FATAL,'init_aquqpe_temp: abs(nlatin)>90.')

    if (size(prsl)/=size(ta)) call mpp_error(FATAL, 'init_aquqpe_temp: size(prsl)/=size(ta)')

    inlat = abs(latin)
   
    do i = 1, size(ta1) 
        call interp_vert(temp(:,i), ta1(i:i), lat_temp, inlat)
    end do

    call interp_vert(ta1, ta, level_temp, prsl)

end subroutine aquape_init_temp


subroutine aquape_o3(latin, prsl, ozone)
    real, intent(in) :: latin, prsl(:) ! pressure levels in hPa(millibar)
    real, intent(out) :: ozone(:)
    
    integer :: j1, j2, i
    real :: ozone1(nlev_o3), w1, w2, w12

    real :: inlat(1)

    if (abs(latin)>90.) call mpp_error(FATAL,'init_aquqpe_temp: abs(nlatin)>90.')

    if (size(prsl)/=size(ozone)) call mpp_error(FATAL, 'init_aquqpe_temp: size(prsl)/=size(ozone)')

    inlat = abs(latin)
   
    do i = 1, size(ozone1) 
        call interp_vert(o3(:,i), ozone1(i:i), lat_o3, inlat)
    end do

    call interp_vert(ozone1, ozone, level_o3, prsl)

end subroutine aquape_o3


subroutine get_aquape_sst(sst, lat)
    real, dimension(:,:), intent(out) :: sst(:,:)
    real, dimension(:,:), intent(in) :: lat(:,:)
    integer :: i, j

    if (any(shape(sst)-shape(lat)/=0)) call mpp_error(FATAL,'get_aquape_sst: shape mismatch shape(sst)/=shape(lat)')

    sst = 0.
    do i = 1, size(sst,2)    
        do j = 1, size(sst,1)
            if (abs(lat(j,i)*PI/180.)>=PI/3.) cycle
            sst(j,i) = 27.*(1-sin((3.*lat(j,i)*PI/180.)/2.))
        end do
    end do

    sst = sst + KELVIN

end subroutine get_aquape_sst

!--------------------------------------------------------------------------------   
subroutine interp_vert (fldin,fldout,axin,axout)
!--------------------------------------------------------------------------------   
    real, dimension(:), intent(in) :: fldin, axin, axout
    real, dimension(:), intent(out) :: fldout

    integer :: i1, i2, ni, j
    real :: w1, w2, tmp(size(axin)), w12

    ni = size(axin)

    do j = 1, size(axout)
        tmp = axin-axout(j)
        i1 = 0; i2 = 0
        if (any(tmp>=0)) i1 = minloc(tmp, 1, mask=tmp>=0.)
        if (any(tmp<=0)) i2 = maxloc(tmp, 1, mask=tmp<=0.)

        if(i1<=0.and.i2<=0) call mpp_error(FATAL,'interp_vert: both i1 and i2 <= 0')

        if (i1<=0) i1=i2
        if (i2<=0) i2=i1

        w1 = abs(axout(j)-axin(i2))
        w2 = abs(axout(j)-axin(i1))

        if (w1==0.or.w2==0.) then
            w1 = 1.; w2 = 1.
        endif

        w12 = w1 + w2

        w1 = w1/w12; w2=w2/w12

        fldout(j) = fldin(i1)*w1+fldin(i2)*w2
    end do

    return

end subroutine interp_vert

end module aqua_planet_mod

