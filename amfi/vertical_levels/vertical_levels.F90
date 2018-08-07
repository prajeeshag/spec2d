module vertical_levels_mod

use mpp_mod, only : mpp_error, FATAL, NOTE, WARNING, mpp_pe, mpp_root_pe
use fms_io_mod, only : read_data
use constants_mod, only : RDGAS, CP_AIR
use omegtes_mod, only : omegtes

implicit none
private

public :: init_vertical_levels, get_ak_bk, get_pressure_at_levels, &
          get_vertical_vel

real, parameter :: rk = RDGAS/CP_AIR
real, parameter :: rk1 = rk + 1.0, rkr = 1.0/rk, r100=100.0, pt01=0.01

real, allocatable :: ak(:), bk(:)
real, allocatable :: ck(:), dbk(:)
real, allocatable :: bkl(:), si(:)
real, allocatable :: sl(:)

integer :: nlevs

logical :: initialized=.false.

contains

!--------------------------------------------------------------------------------   
subroutine  init_vertical_levels(nlevs_in)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: nlevs_in
    integer :: k
    real :: psurff = 101.3

    nlevs = nlevs_in

    call set_ak_bk(nlevs_in)

    if (mpp_pe()==mpp_root_pe()) then
        print *, 'ak=', ak
        print *, 'bk=', bk
    endif

    allocate(dbk(nlevs),bkl(nlevs),ck(nlevs))
    allocate(si(nlevs+1),sl(nlevs))

    do k = 1, nlevs
        dbk(k) = bk(k+1)-bk(k)
        bkl(k) = (bk(k+1)+bk(k))*0.5
        ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
    enddo

    do k = 1, nlevs+1
        si(nlevs+2-k) = ak(k)/psurff + bk(k)
    enddo

    do k = 1, nlevs
        sl(k) = 0.5*(si(k)+si(k+1))
    enddo

    initialized = .true.
    return
end subroutine init_vertical_levels

!--------------------------------------------------------------------------------   
subroutine get_vertical_vel(pgr,dphi,dlam,div,u,v,vvel)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(in), dimension(:,:) :: pgr, dphi, dlam
    real, intent(in), dimension(:,:,:) :: div, u, v
    real, intent(out), dimension(:,:,:) :: vvel

    integer :: km, imax

    if(.not.initialized) call mpp_error('get_vertical_vel','module not initialized',fatal) 

    km = size(vvel,1)
    imax = size(pgr,2) * size(pgr,1)

    call omegtes_drv(imax, km, pgr, dphi, dlam, div, u, v, vvel)
    
    return

end subroutine get_vertical_vel

!--------------------------------------------------------------------------------   
subroutine omegtes_drv(imax,km,pgr, dphi, dlam, div, u, v, vvel)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: km, imax
    real, intent(in), dimension(imax) :: pgr, dphi, dlam
    real, intent(in), dimension(km,imax) :: div, u, v
    real, intent(out), dimension(km,imax) :: vvel
    
    integer :: i 
    
    do i = 1, imax
        call omegtes(km, ak, bk, dbk, ck, pgr(i), dphi(i), &
                 dlam(i), div(:,i), u(:,i), v(:,i), vvel(:,i))
    enddo
end subroutine omegtes_drv


!--------------------------------------------------------------------------------   
subroutine get_pressure_at_levels(pgr,prsi,prsl,prsik,prslk)
!--------------------------------------------------------------------------------   
    implicit none
    real, intent(in)  :: pgr(:,:)
    real, intent(out), optional :: prsi(:,:,:), prsik(:,:,:)
    real, intent(out), optional :: prsl(:,:,:), prslk(:,:,:)

    real, dimension(nlevs+1,size(pgr,1),size(pgr,2)) :: prsi1, prsi1k
    real, dimension(nlevs  ,size(pgr,1),size(pgr,2)) :: prsl1, prsl1k
    real :: tem(size(pgr,1),size(pgr,2))
    integer :: k

    if(.not.initialized) call mpp_error('vertical_levels_mod','module not initialized',fatal)

    if(.not.present(prsi).and. &
       .not.present(prsl).and. &
       .not.present(prsik).and. &
       .not.present(prslk)) return

    do k=1,nlevs+1
        prsi1(nlevs+2-k,:,:) = ak(k) + bk(k) * pgr(:,:)
    enddo

    if (present(prsi)) prsi = prsi1

    if(.not.present(prsl).and. &
       .not.present(prsik).and. &
       .not.present(prslk)) return
    
    prsi1k(:,:,:) = (prsi1(:,:,:)*pt01) ** rk

    if (present(prsik)) prsik = prsi1k
   
    if(.not.present(prsl).and. &
       .not.present(prslk)) return
    
    do k=1,nlevs
        tem(:,:) = rk1 * (prsi1(k,:,:) - prsi1(k+1,:,:))
        prsl1k(k,:,:) = (prsi1k(k,:,:)*prsi1(k,:,:)-prsi1k(k+1,:,:)*prsi1(k+1,:,:))/tem
        prsl1(k,:,:)    = r100 * prsl1k(k,:,:) ** rkr
    enddo

    if(present(prsl)) prsl = prsl1
    if(present(prslk)) prslk = prsl1k
    
    return
end subroutine get_pressure_at_levels


!--------------------------------------------------------------------------------   
subroutine get_ak_bk(ak_out,bk_out,ck_out,dbk_out,bkl_out,si_out,sl_out)
!--------------------------------------------------------------------------------   
    real, intent(out), optional :: ak_out(:), bk_out(:)
    real, intent(out), optional :: ck_out(:), dbk_out(:)
    real, intent(out), optional :: bkl_out(:), si_out(:)
    real, intent(out), optional :: sl_out(:)
   
    if(.not.initialized) call mpp_error('get_ak_bk','module not initialized',fatal) 

    if(present(ak_out)) ak_out(:)=ak 
    if(present(bk_out)) bk_out(:)=bk 
    if(present(ck_out)) ck_out(:)=ck 
    if(present(dbk_out)) dbk_out(:)=dbk 
    if(present(bkl_out)) bkl_out(:)=bkl 
    if(present(si_out)) si_out(:)=si 
    if(present(sl_out)) sl_out(:)=sl 

    return
end subroutine


!--------------------------------------------------------------------------------   
subroutine set_ak_bk(nlevs_in)
!-------------------------------------------------------------------------------- 
    integer, intent(in) :: nlevs_in

    select case(nlevs_in)
    case(64)
        allocate(ak(65),bk(65)) 
        ak = &
         [0.00000000000,   0.06424700165,   0.13778999329,   0.22195799255, &
         0.31826599121,   0.42843399048,   0.55442401123,   0.69845697021, &
         0.86305798340,   1.05107995605,   1.26575195313,   1.51071105957, &
         1.79005102539,   2.10836596680,   2.47078808594,   2.88303808594, &
         3.35145996094,   3.88305200195,   4.48549316406,   5.16714599609, &
         5.93704980469,   6.80487402344,   7.77714990234,   8.83253710938, &
         9.93661425781,  11.05485253906,  12.15293652344,  13.19706542969, &
        14.15431640625,  14.99307421875,  15.68348925781,  16.19796679688, &
        16.51173632813,  16.61160351562,  16.50314453125,  16.19731542969, &
        15.70889257813,  15.05634179687,  14.26143457031,  13.34867089844, &
        12.34449023437,  11.27634765625,  10.17171191406,   9.05705078125, &
         7.95690820312,   6.89311718750,   5.88420605469,   4.94502880859, &
         4.08661401367,   3.31621704102,   2.63755297852,   2.05114990234, &
         1.55478894043,   1.14398803711,   0.81248901367,   0.55271997070, &
         0.35622299194,   0.21401499939,   0.11689900208,   0.05571200180, &
         0.02151600075,   0.00574100018,   0.00057499999,   0.00000000000, &
         0.00000000000]

        bk = &
          [0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000, &
          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000, &
          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000, &
          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000, &
          0.00000000000,   0.00000000000,   0.00000000000,   0.00000000000, &
          0.00000000000,   0.00000000000,   0.00003697000,   0.00043106001, &
          0.00163591001,   0.00410671020,   0.00829401985,   0.01463712007, &
          0.02355588041,   0.03544161841,   0.05064684153,   0.06947457790, &
          0.09216690809,   0.11881218851,   0.14926877618,   0.18329623342, &
          0.22057019174,   0.26068544388,   0.30316412449,   0.34746849537, &
          0.39301824570,   0.43921080232,   0.48544332385,   0.53113484383, &
          0.57574665546,   0.61879962683,   0.65988701582,   0.69868290424, &
          0.73494523764,   0.76851469278,   0.79930973053,   0.82731884718, &
          0.85259068012,   0.87522357702,   0.89535498619,   0.91315102577, &
          0.92879730463,   0.94249105453,   0.95443409681,   0.96482759714, &
          0.97386759520,   0.98174226284,   0.98862659931,   0.99467116594, &
          1.00000000000]
    case default
        call mpp_error('set_ak_bk', 'levs not supported', FATAL)
    end select

end subroutine set_ak_bk

end module vertical_levels_mod

