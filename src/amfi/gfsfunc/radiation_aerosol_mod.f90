module radiation_aerosol_mod

  use module_radsw_main, only : NBDSW
  use module_radlw_main, only : NBDLW
  use data_override_mod, only : data_override
  use mpp_mod, only : handle_error=>mpp_error, fatal, warning, note, mpp_pe, &
    mpp_root_pe
  use physcons_mod, only : con_pi, con_rd, con_fvirt, con_g,con_rog,&
                       con_t0c, con_c, con_boltz, con_plnk
  use diag_manager_mod, only: diag_axis_init, reg_df=>register_diag_field, send_data
  use time_manager_mod, only: time_type
  use mpp_domains_mod, only : domain2D, mpp_get_compute_domain
  use fms_mod, only : open_namelist_file, close_file

   
  implicit none
  private

  integer, parameter, public :: NF_AESW = 3     ! num of output fields for sw rad
  integer, parameter, public :: NF_AELW = 3     ! num of output fields for lw rad

  real, allocatable :: aod_tropo_sw_in(:,:,:), asy_tropo_sw_in(:,:,:), ssa_tropo_sw_in(:,:,:)

  real, allocatable :: aod_anthrop_sw_in(:,:,:,:), asy_anthrop_sw_in(:,:,:,:), ssa_anthrop_sw_in(:,:,:,:)

  real, allocatable :: aod_strato_sw_in(:,:,:,:), asy_strato_sw_in(:,:,:,:), ssa_strato_sw_in(:,:,:,:)

  real, allocatable :: extcoef_sw_in(:,:,:)

  integer :: nlevs_of_tropo_in = 35, nlevs_of_strato_in = 70, nlevs_of_anthrop_in = 64
  integer :: id_aero
  logical :: dz_in_equal=.true.
  real :: z1_in=-1.0, rdz_clim
  real :: dz_in_strato = 500.0, z0_in_strato = 5000.0, dz_in_tropo = 500.0, z0_in_tropo = 250.0
  character (len=32) :: vert_interp_opt='nearest'
  
  logical :: initialized=.false., use_this_aerosol=.true., strato=.true., tropo = .true., anthrop = .true., natural = .true.
  
  integer :: id_aod_tropo_sw_in=0, id_ssa_tropo_sw_in=0, id_asy_tropo_sw_in=0, id_extcoef_sw_in=0
  integer :: id_zhaero=0, id_oro=0, id_ext_coef_sw=0

  integer :: id_aod_anthrop_sw_in1=0, id_aod_anthrop_sw_in2=0, id_aod_anthrop_sw_in3=0
  integer :: id_aod_anthrop_sw_in4=0, id_aod_anthrop_sw_in5=0, id_aod_anthrop_sw_in6=0
  integer :: id_aod_anthrop_sw_in7=0, id_aod_anthrop_sw_in8=0, id_aod_anthrop_sw_in9=0
  integer :: id_aod_anthrop_sw_in10=0, id_aod_anthrop_sw_in11=0, id_aod_anthrop_sw_in12=0
  integer :: id_aod_anthrop_sw_in13=0, id_aod_anthrop_sw_in14=0

  integer :: id_ssa_anthrop_sw_in1=0, id_ssa_anthrop_sw_in2=0, id_ssa_anthrop_sw_in3=0
  integer :: id_ssa_anthrop_sw_in4=0, id_ssa_anthrop_sw_in5=0, id_ssa_anthrop_sw_in6=0
  integer :: id_ssa_anthrop_sw_in7=0, id_ssa_anthrop_sw_in8=0, id_ssa_anthrop_sw_in9=0
  integer :: id_ssa_anthrop_sw_in10=0, id_ssa_anthrop_sw_in11=0, id_ssa_anthrop_sw_in12=0
  integer :: id_ssa_anthrop_sw_in13=0, id_ssa_anthrop_sw_in14=0
!
  integer :: id_asy_anthrop_sw_in1=0, id_asy_anthrop_sw_in2=0, id_asy_anthrop_sw_in3=0
  integer :: id_asy_anthrop_sw_in4=0, id_asy_anthrop_sw_in5=0, id_asy_anthrop_sw_in6=0
  integer :: id_asy_anthrop_sw_in7=0, id_asy_anthrop_sw_in8=0, id_asy_anthrop_sw_in9=0
  integer :: id_asy_anthrop_sw_in10=0, id_asy_anthrop_sw_in11=0, id_asy_anthrop_sw_in12=0
  integer :: id_asy_anthrop_sw_in13=0, id_asy_anthrop_sw_in14=0
!
  integer :: id_aod_strato_sw_in1=0, id_aod_strato_sw_in2=0, id_aod_strato_sw_in3=0
  integer :: id_aod_strato_sw_in4=0, id_aod_strato_sw_in5=0, id_aod_strato_sw_in6=0
  integer :: id_aod_strato_sw_in7=0, id_aod_strato_sw_in8=0, id_aod_strato_sw_in9=0
  integer :: id_aod_strato_sw_in10=0, id_aod_strato_sw_in11=0, id_aod_strato_sw_in12=0
  integer :: id_aod_strato_sw_in13=0, id_aod_strato_sw_in14=0

  integer :: is, ie, js, je, ilen, jlen, nlev
  
  real, allocatable:: z_strato_in(:), z_tropo_in(:)
  real, allocatable, dimension(:,:) :: oro

  namelist /radiation_aerosol_nml/ dz_in_equal, nlevs_of_tropo_in, z1_in, vert_interp_opt, &
    nlevs_of_strato_in, z0_in_strato, dz_in_strato, z0_in_tropo, dz_in_tropo, &
    use_this_aerosol, strato, tropo, nlevs_of_anthrop_in, anthrop, natural

  public :: set_aerosols, init_aerosol_mod 
  
  contains

    subroutine init_aerosol_mod(Time,domain_in,nlev_in,axes_in)
      type(time_type), intent(in) :: Time
      type(domain2D), intent(in) :: domain_in 
      integer, intent(in) :: nlev_in
      integer, intent(in) :: axes_in(:)
      integer :: i, unit
      character (len=32) :: name

      if (initialized) return

      unit = open_namelist_file()
      read(unit,nml=radiation_aerosol_nml)
      if (mpp_pe()==mpp_root_pe()) write(*,nml=radiation_aerosol_nml)
      call close_file(unit)

      if (.not.use_this_aerosol) then
        initialized = .true.
        return
      endif

      nlev = nlev_in

      call mpp_get_compute_domain(domain_in,js,je,is,ie)
      jlen = je - js + 1
      ilen = ie - is + 1

      if (nlevs_of_tropo_in <=0) call handle_error(fatal, 'nlevs_of_tropo_in if aerosol data should be greater than 0')
!     if (dz_in_equal) then
!       if (dz_in<=0.0) call handle_error(fatal, 'if dz_in_equal=.true. then dz hould be greter than 0.0')
!       if (z1_in<0.0.and.trim(vert_interp_opt)/='nearest') &
!       call handle_error(fatal, 'if dz_in_equal=.true. then z1_in should be greater than or equal to 0.0')
!     endif

      allocate(aod_tropo_sw_in(NBDSW, js:je, is:ie)) 
      allocate(ssa_tropo_sw_in(NBDSW, js:je, is:ie)) 
      allocate(asy_tropo_sw_in(NBDSW, js:je, is:ie)) 
      allocate(extcoef_sw_in(nlevs_of_tropo_in, js:je, is:ie)) 

      allocate(aod_anthrop_sw_in(nlev, js:je, is:ie, NBDSW))
      allocate(ssa_anthrop_sw_in(nlev, js:je, is:ie, NBDSW))
      allocate(asy_anthrop_sw_in(nlev, js:je, is:ie, NBDSW))

      allocate(aod_strato_sw_in(nlevs_of_strato_in, js:je, is:ie, NBDSW))
!     allocate(ssa_strato_sw_in(js:je, is:ie, nlevs_of_strato_in,NBDSW))
!     allocate(asy_strato_sw_in(js:je, is:ie, nlevs_of_strato_in,NBDSW))
      allocate(z_strato_in(nlevs_of_strato_in))
      allocate(z_tropo_in(nlevs_of_tropo_in))

      rdz_clim = 1.0 / dz_in_tropo     !lay(k+1)-lay(k)

    ! Calculate the heights of stratospheric aerosol file 
    do i=1,nlevs_of_strato_in
      z_strato_in(i) = z0_in_strato + (i-1) * dz_in_strato
    end do

    ! Calculate the heights of tropospheric aerosol file 
    do i=1,nlevs_of_tropo_in
      z_tropo_in(i) = z0_in_tropo + (i-1) * dz_in_tropo
    end do

    if(mpp_pe()==mpp_root_pe()) print *, 'z_strato_in=', z_strato_in
 
 !   do i = 1, NBDSW
  !    write(name,'(A,I2.2)') 'aod',i
   ! enddo

    !id_extcoef_sw_in =  register_var('extcoef_sw_in', 'Extinction Coef', ' ', nlevs_of_tropo_in )
    !id_aod_tropo_sw_in = register_var('aod_tropo_sw_in', 'AOD input', ' ', NBDSW) 
    !id_zhaero = register_var('zhaero', 'zhaero', ' ', levs) 
    !id_oro = register_var('oro','oro',' ')
    !id_ext_coef_sw = register_var('ext_coef_sw','ext_coef_sw',' ',levs)
  
    initialized = .true.

  end subroutine init_aerosol_mod


  
  subroutine get_aerosols(Time)
    !This should be called before grrad, Either in dotstep or in the beging of gloopr.
    ! Get the data from aerosol input files given in atm_data_table.
    type(time_type), intent(in) :: Time
    logical :: ov
    character(len=32) :: varnm
    integer :: i


    if (.not.use_this_aerosol) return
   
    call data_override('ATM','aod_sw',aod_tropo_sw_in,Time,ov,kxy=1)
    if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read aod_sw.') 
    !call send_data(id_aod_tropo_sw_in, aod_tropo_sw_in)
   
    call data_override('ATM','ssa_sw',ssa_tropo_sw_in,Time,ov,kxy=1)  
    if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read ssa_sw.')
    !call send_data(id_ssa_tropo_sw_in, ssa_tropo_sw_in)
   
    call data_override('ATM','asy_sw',asy_tropo_sw_in,Time,ov,kxy=1)  
    if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read asy_sw.')
    !call send_data(id_asy_tropo_sw_in, asy_tropo_sw_in)
   
    call data_override('ATM','extcoef_sw',extcoef_sw_in,Time,ov,kxy=1)  
    if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read extcoef_sw.')
    !call send_data(id_extcoef_sw_in, extcoef_sw_in)


    do i = 1, NBDSW
      write(varnm,*) i
      varnm = 'aod_anthrop_sw'//trim(adjustl(varnm))
      call data_override('ATM',trim(varnm),aod_anthrop_sw_in(:,:,:,i),Time,ov,kxy=1)
      if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read '//trim(varnm))
      !call send_data(id_aod_anthrop_sw_in1, aod_anthrop_sw_in(:,:,:,1))
    end do

    do i = 1, NBDSW
      write(varnm,*) i
      varnm = 'ssa_anthrop_sw'//trim(adjustl(varnm))
      call data_override('ATM',trim(varnm),ssa_anthrop_sw_in(:,:,:,i),Time,ov,kxy=1)
      if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read '//trim(varnm))
      !call send_data(id_ssa_anthrop_sw_in1, ssa_anthrop_sw_in(:,:,:,1))
    end do

    do i = 1, NBDSW
      write(varnm,*) i
      varnm = 'asy_anthrop_sw'//trim(adjustl(varnm))
      call data_override('ATM',trim(varnm),asy_anthrop_sw_in(:,:,:,i),Time,ov,kxy=1)
      if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read '//trim(varnm))
      !call send_data(id_asy_anthrop_sw_in1, asy_anthrop_sw_in(:,:,:,1))
    end do

    if (strato) then
      do i = 1, NBDSW
        write(varnm,*) i
        varnm = 'aod_strato_sw'//trim(adjustl(varnm))
        call data_override('ATM',trim(varnm),aod_strato_sw_in(:,:,:,i),Time,ov,kxy=1)
        if (.not.ov) call handle_error(fatal, 'get_aerosols: Could not read '//trim(varnm))
        !call send_data(id_aod_strato_sw_in1, aod_strato_sw_in(:,:,:,1))
      end do
    endif

  end subroutine get_aerosols


  subroutine set_aerosols (Time, prsi, prsl, tlay, geopl, oro, aerosw, aerolw)
    implicit none
    !  ---  inputs:
    type(time_type), intent(in) :: Time
    real, dimension(1:,js:,is:), intent(in) :: prsi, prsl, tlay
    real, dimension(1:,js:,is:), intent(in) :: geopl
    real, dimension(js:,is:), intent(in) :: oro

    !  ---  outputs:
    real, dimension(1:nlev,1:NBDSW,1:NF_AESW,js:je,is:ie), intent(out) :: aerosw
    real, dimension(1:nlev,1:NBDLW,1:NF_AELW,js:je,is:ie), intent(out) :: aerolw 

    !  ---  locals:

    integer :: i, i1, i2, j1, j2, k, m, m1, kp,jk,j,l
    integer :: jwl, jl, isc, iec,nlev_start

    integer, dimension(js:je, is:ie):: kindex

    real :: zdeltag(nlev,js:je,is:ie), min_index,zhif(nlev+1,js:je,is:ie)

    real :: zh(nlev,js:je,is:ie), delp(nlev,js:je,is:ie), wt(nlev,js:je,is:ie), donm(nlev,js:je,is:ie)

    real, dimension(nlev, js:je, is:ie):: ext_coef_sw
    real, dimension(js:je,is:ie):: ext_coef_sum 
    real, dimension(nlev,js:je,is:ie,NBDSW):: aod_strato_sw, ssa_strato_sw, asy_strato_sw
    real, dimension(nlev,js:je,is:ie,NBDSW):: aod_anthrop_sw, ssa_anthrop_sw, asy_anthrop_sw
    real, dimension(nlev,js:je,is:ie,NBDSW):: aod_natural_sw, ssa_natural_sw, asy_natural_sw

    !===>  ...  begin here

    aerosw = 0.
    aerolw = 0.

    if (.not.use_this_aerosol) return

    ! (i) calculate altitude above NN and layer thickness in 
    !     gfs for altitude profiles

    call get_aerosols(Time)

    do k = 1, nlev
      delp(k,js:je,is:ie) = prsi(k,js:je,is:ie) - prsi(k+1,js:je,is:ie)
    enddo

    zdeltag(1:nlev,js:je,is:ie) = (delp(1:nlev,js:je,is:ie) * & 
      tlay(1:nlev,js:je,is:ie) / &
      prsl(1:nlev,js:je,is:ie))*con_rog

    do k = 1, nlev
      zh(k,js:je,is:ie) = (geopl(k,js:je,is:ie)/con_g) + oro(js:je,is:ie)
    end do

    !interface height
    zhif(0,js:je,is:ie)=0.0
    do k=2,nlev
      zhif(k,js:je,is:ie) = ( zh(k+1,js:je,is:ie) + zh(k,js:je,is:ie) ) * 0.5
    end do

    ! put the stratospheric aerosols to the nearest model
    ! levels
    aod_strato_sw = 0.0
    do i = is, ie
      do j = js, je
        nlev_start=1
        do l=1, nlev-1
          do k = nlev_start, nlevs_of_strato_in
            if (z_strato_in(k).gt.zhif(l+1,j,i)) exit
            if ((z_strato_in(k) > zhif(l,j,i)).and.(z_strato_in(k) <= zhif(l+1,j,i))) then
              if ((z_strato_in(k) > zhif(l,j,i)).and.(z_strato_in(k) <= zh(l,j,i))) then
                aod_strato_sw(l,j,i,1:NBDSW) = aod_strato_sw(l,j,i,1:NBDSW)+ & 
                  ((z_strato_in(k)-zhif(l,j,i)) / (zh(l,j,i)-zhif(l,j,i))) & 
                  * aod_strato_sw_in(k,j,i,1:NBDSW)
              else
                aod_strato_sw(l,j,i,1:NBDSW)   = aod_strato_sw(l,j,i,1:NBDSW) + &
                  ((zhif(l+1,j,i)-z_strato_in(k)) / (zhif(l+1,j,i) - &
                  zh(l,j,i))) * aod_strato_sw_in(k,j,i,1:NBDSW)
              endif
            endif
          end do
          nlev_start = k
        end do
      end do
    end do

    !-----Multiply extinction coefficient with zh
    do j = js, je
      do i=is, ie
        do jwl=1,nbdsw
          do jk=1,nlev
            aod_strato_sw(jk,j,i,NBDSW) = (aod_strato_sw(jk,j,i,NBDSW)) *  zdeltag(jk,j,i)  * 0.001
          end do
        end do
      end do
    end do

    ! Put the stratospheric aerosols to zero wherever tropospheric aerosols are
    ! present
  
    do l=1,nlev
      do j = js, je
        do i=is, ie
          if (zh(l,j,i) < z_tropo_in(nlevs_of_tropo_in).or. &
            (zh(l,j,i) > z_strato_in(nlevs_of_strato_in))) then
            aod_strato_sw(l,j,i,1:NBDSW) = 0.0
          end if
        end do
      end do
    end do

    ! Put the anthropogenic aerosols to the model heights 
    !do i = is, ie
    !  do j = js, je
    !    do jwl = 1, nbdsw
    !      do l = 1, nlev
             aod_anthrop_sw = aod_anthrop_sw_in
             ssa_anthrop_sw = ssa_anthrop_sw_in
             asy_anthrop_sw = asy_anthrop_sw_in
             aod_anthrop_sw = nint(aod_anthrop_sw * 1000000.0) * 1E-6
    !      end do
    !    end do
    !  end do
    !end do

    ! (ii) Vertical indexing of extinction coeff in for gfs grid (following Echam) 

    ext_coef_sw(1:nlev, js:je, is:ie)=0.0

    do jk=1,nlev
      kindex(js:je,is:ie)=int(zh(jk,js:je,is:ie)*rdz_clim+0.50)
      where(kindex<1) kindex=1
      do j=js,je
        do i=is,ie
          if (kindex(j,i) > 0 .and. kindex(j,i) <= nlevs_of_tropo_in) then
            ext_coef_sw(jk,j,i)= extcoef_sw_in(kindex(j,i),j,i)
          end if
        end do
      end do
    end do

    ! normalize height profile for all modes(anthropogenic and natural)
    ! Needed even if exticntion data is in model levels

      ext_coef_sum(js:je, is:ie) = 0.0

   !   do jk=1,nlay
   !     ext_coef_sum(1:im)=ext_coef_sum(1:im)+                                  &
   !  &                   ext_coef_sw(1:im,jk)*zdeltag(1:im,jk)
   !   enddo

   !   where (ext_coef_sum(1:im) <= 0.0)
   !     ext_coef_sum(1:im)=1.0
   !   end where

      do jk=1,nlev
        ! ext_coef_sw(1:im,jk)=zdeltag(1:im,jk)*ext_coef_sw(1:im,jk)/ext_coef_sum(1:im)
        ext_coef_sw(jk,js:je,is:ie)=zdeltag(jk,js:je,is:ie)*ext_coef_sw(jk,js:je,is:ie)
        ext_coef_sum(js:je,is:ie)=ext_coef_sum(js:je,is:ie)+ext_coef_sw(jk,js:je,is:ie)
      end do

      where (ext_coef_sum(js:je,is:ie) <= 0.0)
        ext_coef_sum(js:je,is:ie)=1.0
      end where

      do jk=1,nlev
        ext_coef_sw(jk,js:je,is:ie)=ext_coef_sw(jk,js:je,is:ie)/ext_coef_sum(js:je,is:ie)
      enddo
      
      aerosw(1:nlev,1:NBDSW,1:NF_AESW,js:je,is:ie) = 0.0
   
! calculate optical properties
! aerosol optical depth 
      do jk=1,nlev
        do jwl=1,nbdsw
          if(use_this_aerosol) then
            if(tropo) then
              do j=js, je
                do i=is, ie
                  if (zh(jk,j,i) < z_tropo_in(nlevs_of_tropo_in)) then
                    if(natural) then
                      aod_natural_sw(jk,j,i,jwl) = aod_tropo_sw_in(jwl,j,i) * ext_coef_sw(jk,j,i)
                      aerosw(jk,jwl,1,j,i) = aerosw(jk,jwl,1,j,i) + aod_natural_sw(jk,j,i,jwl)
                    end if ! end if for natural

                    if(anthrop) then
                      aerosw(jk,jwl,1,j,i) = aerosw(jk,jwl,1,j,i) + aod_anthrop_sw(jk,j,i,jwl)
                    end if ! end if for anthrop

                  end if ! end if to make sure tropospheric aerosol is below 17 km
                end do ! end do for im loop
              end do ! end do for im loop
            end if ! end if for tropo

            if (strato) then
              do j=js,je
                do i=is,ie
                  if (zh(jk,j,i) > z_tropo_in(nlevs_of_tropo_in).and.(zh(jk,j,i) < z_strato_in(nlevs_of_strato_in))) then
                      aerosw(jk,jwl,1,j,i) = aerosw(jk,jwl,1,j,i) + aod_strato_sw(jk,j,i,jwl)
                  end if ! end if to make sure that the stratospheric aerosol goes only in stratosphere 
                end do ! end do for im
              end do ! end do for im
            end if ! end if for strato
          end if ! end if for use_this_aerosol redundant, as control will not enter this subroutine if false

          if(tropo) then  
            do j=js,je
              do i=is,ie
                if (zh(jk,j,i) < z_tropo_in(nlevs_of_tropo_in)) then
                  ssa_natural_sw(jk,j,i,jwl) = ssa_tropo_sw_in(jwl,j,i)
                  asy_natural_sw(jk,j,i,jwl)  = asy_tropo_sw_in(jwl,j,i)

                  if(natural.and.anthrop.and.(.not.(aerosw(jk,jwl,1,j,i)==0.0))) then
                    aerosw(jk,jwl,2,j,i) = ((ssa_natural_sw(jk,j,i,jwl) * &
                      aod_natural_sw(jk,j,i,jwl)) + (ssa_anthrop_sw(jk,j,i,jwl) * &
                      aod_anthrop_sw(jk,j,i,jwl))) / (aod_natural_sw(jk,j,i,jwl)+ &
                      aod_anthrop_sw(jk,j,i,jwl))

                    aerosw(jk,jwl,3,j,i) = ((ssa_natural_sw(jk,j,i,jwl) * &
                      aod_natural_sw(jk,j,i,jwl) * asy_natural_sw(jk,j,i,jwl)) + &
                      (ssa_anthrop_sw(jk,j,i,jwl) * aod_anthrop_sw(jk,j,i,jwl) * &
                       asy_anthrop_sw(jk,j,i,jwl))) / ((ssa_natural_sw(jk,j,i,jwl) * &
                       aod_natural_sw(jk,j,i,jwl)) +   (ssa_anthrop_sw(jk,j,i,jwl) * &
                       aod_anthrop_sw(jk,j,i,jwl)))
                  end if ! end if for natural and anthrop 

                  if(anthrop.and.(.not.(natural))) then
                    aerosw(jk,jwl,2,j,i) = ssa_anthrop_sw(jk,j,i,jwl)
                    aerosw(jk,jwl,3,j,i) = asy_anthrop_sw(jk,j,i,jwl)
                  end if ! end if for anthrop and no natural

                  if(natural.and.(.not.(anthrop))) then
                    aerosw(jk,jwl,2,j,i) = ssa_natural_sw(jk,j,i,jwl)
                    aerosw(jk,jwl,3,j,i) = asy_natural_sw(jk,j,i,jwl)
                  end if ! end if for natural and no anthrop
                end if ! end if for limiting tropospheric aerosols to below 17 km
              end do ! end do for im
            end do ! end do for im
          end if ! end if for tropo
        end do! end do for nbdsw

        if(strato) then
          do j=js,je
            do i=is,ie
              if (zh(jk,j,i) > z_tropo_in(nlevs_of_tropo_in).and.(zh(jk,j,i) < z_strato_in(nlevs_of_strato_in))) then
                aerosw(jk,1:nbdsw,2,j,i)  = 0.99 
                aerosw(jk,1:nbdsw,3,j,i)  = 0.6366 
              end if ! end if for ssa and asy in stratosphere  
            end do !end do for im 
          end do !end do for im 
        end if ! end if for strato

        where (aerosw(jk,1:nbdsw,1,js:je,is:ie) == 0.0)
          aerosw(jk,1:nbdsw,2,js:je,is:ie) = 1.0
          aerosw(jk,1:nbdsw,3,js:je,is:ie) = 0.0
        end where

      enddo ! end do for jk/nlay

      aerolw(:,:,:,:,:) = 0.0

  end subroutine set_aerosols

end module radiation_aerosol_mod 
