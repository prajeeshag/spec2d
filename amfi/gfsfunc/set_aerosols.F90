module set_aerosols_mod

use fms_mod, only : open_namelist_file, close_file
use module_radsw_parameters, only : NBDSW
use module_radlw_parameters, only : NBDLW
use data_override_mod, only : data_override
use mpp_mod, only : handle_error=>mpp_error, fatal, warning, note
use physcons_mod, only : con_pi, con_rd, con_fvirt, con_g,con_rog,&
                     con_t0c, con_c, con_boltz, con_plnk
use time_manager_mod, only: time_type
 
implicit none
private

public :: set_aerosols, init_set_aerosols

integer, parameter, public :: NF_AESW = 3 
integer, parameter, public :: NF_AELW = 3

real, allocatable :: aod_tropo_sw_in(:,:,:), asy_tropo_sw_in(:,:,:), ssa_tropo_sw_in(:,:,:)
real, allocatable :: extcoef_sw_in(:,:,:)
real, allocatable :: z_tropo_in(:)
real :: z1_in=-1.0, rdz_clim
real :: dz_in_tropo = 0.0, z0_in_tropo =0.0

integer :: nlevs_of_tropo_in = 0

integer :: js, je, is, ie, nlay

logical :: dz_in_equal=.true.

logical :: initialized=.false., use_this_aerosol=.true.

namelist /radiation_aerosol_nml/ dz_in_equal, nlevs_of_tropo_in, z1_in, z0_in_tropo, &
                                 dz_in_tropo, use_this_aerosol

contains


!--------------------------------------------------------------------------------   
subroutine init_set_aerosols(isc,iec,jsc,jec,nlev)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: isc, iec, jsc, jec, nlev
    integer :: i, unit
    character (len=32) :: name

    if (initialized) return

    is = isc; ie = iec
    js = jsc; je = jec
    nlay = nlev

    unit = open_namelist_file()
    read(unit,nml=radiation_aerosol_nml)
    call close_file(unit)

    if (.not.use_this_aerosol) then
      initialized = .true.
      return
    endif
    
    if (nlevs_of_tropo_in <=0) call handle_error(fatal, 'nlevs_of_tropo_in if aerosol data should be greater than 0')

    allocate(aod_tropo_sw_in(NBDSW,js:je,is:ie)) 
    allocate(ssa_tropo_sw_in(NBDSW,js:je,is:ie)) 
    allocate(asy_tropo_sw_in(NBDSW,js:je,is:ie)) 
    allocate(extcoef_sw_in(nlevs_of_tropo_in,js:je,is:ie)) 

    allocate(z_tropo_in(nlevs_of_tropo_in))

    rdz_clim = 1.0 / dz_in_tropo

      do i=1,nlevs_of_tropo_in
        z_tropo_in(i) = z0_in_tropo + (i-1) * dz_in_tropo
      end do

    initialized = .true.

end subroutine init_set_aerosols


subroutine get_aerosols(Time)
    type(time_type), intent(in) :: Time
    logical :: overriden

    if(.not.initialized) call handle_error(FATAL,'set_aerosol_mod: module not initilialized')

    if (.not.use_this_aerosol) return
 
    overriden = .true.
     
    call data_override('ATM','aod_sw',aod_tropo_sw_in,Time,overriden)
    if (.not.overriden) call handle_error(fatal, 'get_aerosols: Could not read aod_sw.') 
 
    call data_override('ATM','ssa_sw',ssa_tropo_sw_in,Time,overriden)  
    if (.not.overriden) call handle_error(fatal, 'get_aerosols: Could not read ssa_sw.')
 
    call data_override('ATM','asy_sw',asy_tropo_sw_in,Time,overriden)  
    if (.not.overriden) call handle_error(fatal, 'get_aerosols: Could not read asy_sw.')
 
    call data_override('ATM','extcoef_sw',extcoef_sw_in,Time,overriden)  
    if (.not.overriden) call handle_error(fatal, 'get_aerosols: Could not read extcoef_sw.')

end subroutine get_aerosols


!--------------------------------------------------------------------------------   
subroutine set_aerosols (Time, prsi, prsl, tlay, geopl, oro, aerosw, aerolw)
!--------------------------------------------------------------------------------   
    implicit none
    type(time_type), intent(in) :: Time
    real, dimension(nlay+1,js:je,is:ie), intent(in) :: prsi
    real, dimension(nlay,js:je,is:ie), intent(in) :: prsl, tlay
    real, dimension(nlay,js:je,is:ie), intent(in) :: geopl
    real, dimension(js:je,is:ie),   intent(in) :: oro

    real, dimension(nlay,NBDSW,NF_AESW,js:je,is:ie), intent(out) :: aerosw
    real, dimension(nlay,NBDLW,NF_AELW,js:je,is:ie), intent(out) :: aerolw

    integer :: i, i1, i2, j1, j2, k, m, m1, kp,jk,j,l
    integer ::jwl, jl, isc, iec

    integer, dimension(js:je,is:ie):: kindex

    real :: zdeltag(nlay,js:je,is:ie), min_index

    real :: zh(nlay,js:je,is:ie), delp(nlay,js:je,is:ie)

    real, dimension(nlay,js:je,is:ie) :: ext_coef_sw
    real, dimension(js:je,is:ie) :: ext_coef_sum 

    if(.not.initialized) call handle_error(FATAL,'set_aerosol_mod: module not initilialized')

    if (.not.use_this_aerosol) return

    aerosw = 0.
    aerolw = 0.

    call get_aerosols(Time)

    ! (i) calculate altitude above NN and layer thickness in 
    !     gfs for altitude profiles

    do k=1,nlay
      delp(k,:,:)=prsi(k,:,:)-prsi(k+1,:,:)
    enddo

    zdeltag = (delp*tlay/prsl)*con_rog

    DO k=1,nlay
      zh(k,:,:) = geopl(k,:,:)/con_g + oro(:,:)
    END DO

    ! (ii) Vertical indexing of extinction coeff in for gfs grid (following
    ! Echam) 

    ext_coef_sw=0.0

    do jk=1,nlay
      kindex=int(zh(jk,:,:)*rdz_clim+0.50)
      where(kindex<1) kindex=1
      do i = is, ie
          do j = js, je
              if (kindex(j,i) > 0 .and. kindex(j,i) <= nlevs_of_tropo_in) then
                  ext_coef_sw(jk,j,i)= extcoef_sw_in(kindex(j,i),j,i)
              end if
          end do
      end do
    end do

    ! normalize height profile for all modes(anthropogenic and natural)
    ! Needed even if exticntion data is in model levels

    ext_coef_sum=0.0

    ext_coef_sw = zdeltag*ext_coef_sw
    DO jk=1,nlay
      ext_coef_sum=ext_coef_sum+ext_coef_sw(jk,:,:)
    END DO

    where (ext_coef_sum <= 0.0)
      ext_coef_sum = 1.0
    end where

    do jk=1,nlay
      ext_coef_sw(jk,:,:)=ext_coef_sw(jk,:,:)/ext_coef_sum
    enddo
    
    ! calculate optical properties
    ! aerosol optical depth 
    do jk=1,nlay
      do jwl=1,nbdsw
          aerosw(jk,jwl,1,:,:) = aod_tropo_sw_in(jwl,:,:) * ext_coef_sw(jk,:,:)
      end do

      aerosw(jk,1:nbdsw,2,:,:) = ssa_tropo_sw_in(1:nbdsw,:,:)
      aerosw(jk,1:nbdsw,3,:,:)  = asy_tropo_sw_in(1:nbdsw,:,:)

      where (aerosw(:,1:nbdsw,1,:,:) == 0.0)
        aerosw(:,1:nbdsw,2,:,:) = 1.0
        aerosw(:,1:nbdsw,3,:,:) = 0.0
      end where
    enddo

  end subroutine set_aerosols

end module set_aerosols_mod
