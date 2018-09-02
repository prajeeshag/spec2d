module cu_conv_mod

use fms_mod, only : open_namelist_file, close_file, mpp_error, FATAL, NOTE, uppercase

use mersenne_twister, only : random_setseed, random_number

use sascnv_mod, only : sascnv

use sascnvn_mod, only : sascnvn

use constants_mod, only : KELVIN

implicit none
private

public :: cu_conv, init_cu_conv

character (len=8) :: conv_scheme='SAS'
logical :: cloudice=.false.

character (len=8) :: avail_conv_scheme(2) = ['SAS','SASNEW']

integer :: seed0, iseed
integer :: is, ie, js, je

logical :: initialized = .false., debug=.false.

contains

!--------------------------------------------------------------------------------   
subroutine init_cu_conv(isc,iec,jsc,jec)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: isc,iec,jsc,jec
    integer :: unit, stat
    real :: wrk(1)
    namelist/cu_conv_nml/conv_scheme, cloudice, debug

    unit = open_namelist_file()
    read(unit,nml=cu_conv_nml)
    call close_file(unit)

    is = isc
    ie = iec
    js = jsc
    je = jec

    conv_scheme = uppercase(conv_scheme) 
    
    if (.not.any(conv_scheme==avail_conv_scheme)) &
        call mpp_error(FATAL,'cu_conv_mod: conv_scheme '//trim(conv_scheme) &
                       //' not available')

    call mpp_error(NOTE,'Using '//trim(conv_scheme)//' cumulus convection scheme')

    seed0 = 2009+12+12+0
    call random_setseed(seed0)
    call random_number(wrk)
    seed0 = seed0 + nint(wrk(1)*1000.0)
    iseed = mod(100.0*sqrt(1.*3600),1.0e9) + 1 + seed0

    initialized = .true.
    
end subroutine init_cu_conv

!--------------------------------------------------------------------------------   
subroutine cu_conv (dtp, del, prsl, pgr, phil, clw, cli, q0, t0, u0, v0, cldwrk, &
                    rain, kbot, ktop, kcnv, slmsk, vvel)
!--------------------------------------------------------------------------------   
    real, intent(in) :: dtp
    real, dimension(:,:,:), intent(in) :: del, prsl, phil, vvel 
    real, dimension(:,:), intent(in) :: slmsk, pgr
    real, dimension(:,:,:), intent(inout) :: clw, cli, q0, t0, u0, v0
    integer, dimension(:,:), intent(out) :: kbot, ktop, kcnv
    real, dimension(:,:), intent(out) :: rain, cldwrk

    integer :: imax, km, i, j, k
    character(len=512) :: msg

    if(.not.initialized) call mpp_error(FATAL,'cu_conv_mod: module not initialized!')

    imax = size(del,2)*size(del,3)
    km = size(del,1)

    if (debug) then
        do i = 1, size(t0,3)
            do j = 1, size(t0,2)
                do k = 1, size(t0,1)
                    if (t0(k,j,i)>=KELVIN+100..or.t0(k,j,i)<=KELVIN-160.) then
                        write(msg,'(A,3(I4,1x),F15.10)') &
                                'cu_conv: temperature out of range: ', k, js+j-1, is+i-1, t0(k,j,i)
                        call mpp_error(FATAL,trim(msg))
                    endif
                end do
            end do
        end do
    end if

    call sascnv_drv(imax, km, dtp, del, prsl, pgr, phil, &
                    clw, cli, q0, t0, u0, v0, cldwrk, &
                    rain, kbot, ktop, kcnv, slmsk, vvel)

end subroutine cu_conv

!--------------------------------------------------------------------------------   
subroutine sascnv_drv(imax, km, dtp, del, prsl, pgr, phil, &
                      clw, cli, q0, t0, u0, v0, cldwrk, &
                      rain, kbot, ktop, kcnv, slmsk, vvel)
!--------------------------------------------------------------------------------   

    integer, intent(in) :: imax, km
    real, intent(in) :: dtp
    real, dimension(km,imax), intent(in) :: del, prsl, phil, vvel 
    real, dimension(imax), intent(in) :: slmsk, pgr
    real, dimension(km,imax), intent(inout) :: clw, cli, q0, t0, u0, v0
    integer, dimension(imax), intent(out) :: kbot, ktop, kcnv
    real, dimension(imax), intent(out) :: rain, cldwrk

    real, dimension(imax) :: xkt2
    real, dimension(km,2) :: clw1
    real, dimension(km) :: ud_mf, dd_mf, dt_mf
    integer :: i
    integer, parameter :: ncld = 1

    clw1 = -999.0
     
    if (trim(conv_scheme)=='SAS') then
        call random_setseed(iseed)
        call random_number(xkt2)
        do i = 1, imax
            clw1(:,1) = clw(:,i)
            if (cloudice) clw1(:,2) = cli(:,i)
            call sascnv(km, dtp, del(:,i), prsl(:,i), pgr(i), phil(:,i), &
                   clw1, q0(:,i), t0(:,i), u0(:,i), v0(:,i), cldwrk(i), &
                   rain(i), kbot(i), ktop(i), kcnv(i), slmsk(i), &
                   vvel(:,i), xkt2(i), ncld)
            clw(:,i) = clw1(:,1)
            if (cloudice) cli(:,i) = clw1(:,2)
        enddo
    else
        do i = 1, imax
            clw1(:,1) = clw(:,i)
            if (cloudice) clw1(:,2) = cli(:,i)
            call sascnvn(km, dtp, del(:,i), prsl(:,i), pgr(i), phil(:,i), &
                   clw1, q0(:,i), t0(:,i), u0(:,i), v0(:,i), cldwrk(i), &
                   rain(i), kbot(i), ktop(i), kcnv(i), slmsk(i), &
                   vvel(:,i), ncld, ud_mf, dd_mf, dt_mf)
            clw(:,i) = clw1(:,1)
            if (cloudice) cli(:,i) = clw1(:,2)
        enddo
    endif

end subroutine sascnv_drv

end module cu_conv_mod
