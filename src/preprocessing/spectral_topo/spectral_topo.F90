
#ifndef test_gauss_legendre
program main

  use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error, &
                      mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_max, &
                      mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                      mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather, &
                      mpp_declare_pelist, mpp_set_current_pelist, &
                      mpp_get_current_pelist, mpp_sum
  
  use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain, &
                              mpp_global_field, mpp_get_global_domain
  
  use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, &
                      stdlog, stdout, stderr
  use fms_io_mod, only : fms_io_exit, restart_file_type, register_restart_field, &
          file_exist, field_exist, field_size
  
  use transforms_mod, only : get_spherical_wave, get_lonsP, compute_ucos_vcos, compute_vor_div, &
                             spherical_to_grid, grid_to_spherical, init_transforms, get_latsF, &
                             register_spec_restart, save_spec_restart, get_latsP, &
                             restore_spec_restart, end_transforms, save_wisdom, get_fourier_wave
  
  use ocpack_mod, only : oc_ny, oc_nx, oc_nfour, ocpack_typeP, npack=>oc_npack, get_ocpackP, init_ocpack, &
                         end_ocpack
  
  implicit none

  integer :: layout(2) = [1,1]
  logical :: reduced=.True., packed=.True.

  type(domain2d) :: domain_g1, domain_g2

  integer, allocatable :: yextent1(:), xextent1(:), sph_wave1(:,:), f_wave1(:,:)
  !integer :: num_lat1=8000, trunc1=7998
  integer :: num_lat1=4000, trunc1=3998
  !integer :: num_lat1=2000, trunc1=1998
  !integer :: num_lat1=94, trunc1=62
  integer :: jsp1, jep1, isp1, iep1, ilenp1, jlenp1
  integer :: ocnx1, ocny1, nwaves_oe1

  integer, allocatable :: yextent2(:), xextent2(:), sph_wave2(:,:), f_wave2(:,:)
  integer :: num_lat2=94, trunc2=62
  integer :: jsp2, jep2, isp2, iep2, ilenp2, jlenp2
  integer :: ocnx2, ocny2, nwaves_oe2

  real, allocatable :: topog1(:,:), topog2(:,:), gtopo1(:,:,:), gtopo2(:,:,:)
  complex, allocatable :: stopo1(:,:,:), stopo2(:,:,:)
  integer :: i,j,l,m,j1

  call mpp_init()
  call fms_init()

  layout = 0 
  if (any(layout==0)) then
      layout = [mpp_npes()/50,50]
  elseif (layout(1)*layout(2)/=mpp_npes()) then
      call mpp_error(FATAL,'init_atmos: product of layout should be equal to npes')
  endif

  print *, layout

  allocate(yextent1(layout(1)))
  allocate(xextent1(layout(2)))

  allocate(yextent2(layout(1)))
  allocate(xextent2(layout(2)))

  call init_ocpack(num_lat1, trunc1, layout, yextent=yextent1, xextent=xextent1,&
                  isreduced=reduced, ispacked=packed)

  ocnx1 = oc_nx()
  ocny1 = oc_ny()

  if (mpp_pe()==mpp_root_pe()) print *, ocnx1, ocny1, xextent1, yextent1

  if (sum(yextent1)/=ocny1) then
      print *, 'ocny=', ocny1, 'sum(yextent)=', sum(yextent1), 'yextent=', yextent1
  endif
  if (sum(xextent1)/=ocnx1) then
      print *, 'ocnx=', ocnx1, 'sum(xextent)=', sum(xextent1), 'xextent=', xextent1
  endif

  call mpp_define_domains([1,ocny1,1,ocnx1], layout, domain_g1, &
                          xextent=yextent1, yextent=xextent1, kxy=1)

  call mpp_get_compute_domain(domain_g1, jsp1, jep1, isp1, iep1)
  ilenp1 = iep1 - isp1 + 1
  jlenp1 = jep1 - jsp1 + 1

  call init_transforms(domain_g1,trunc1,nwaves_oe1)

  print *, 'nwaves_oe = ', nwaves_oe1

  allocate(sph_wave1(nwaves_oe1,2))
  allocate(f_wave1(nwaves_oe1,2))

  call get_spherical_wave(sph_wave1)
  call get_fourier_wave(f_wave1)

  allocate(topog1(ocny1,ocnx1))
  allocate(stopo1(1,nwaves_oe1,2))
  allocate(gtopo1(1,jsp1:jep1,isp1:iep1))

  !call read_data('data/topo_16016x8000_Octa.nc','topo',topog1)
  !gtopo1(1,:,:) = topog1(jsp1:jep1,isp1:iep1)
  !call grid_to_spherical(gtopo1,stopo1,do_trunc=.true.)

  !call end_transforms()
  !call end_ocpack()
!
!
!  call init_ocpack(num_lat2, trunc2, layout, yextent=yextent2, xextent=xextent2,&
!                  isreduced=reduced, ispacked=packed)
!
!  ocnx2 = oc_nx()
!  ocny2 = oc_ny()
!  print *, ocnx2, ocny2
!
!  if (sum(yextent2)/=ocny2) then
!      print *, 'ocny=', ocny2, 'sum(yextent)=', sum(yextent2), 'yextent=', yextent2
!  endif
!  if (sum(xextent2)/=ocnx2) then
!      print *, 'ocnx=', ocnx2, 'sum(xextent)=', sum(xextent2), 'xextent=', xextent2
!  endif
!
!  call mpp_define_domains([1,ocny2,1,ocnx2], layout, domain_g2, &
!                          xextent=yextent2, yextent=xextent2, kxy=1)
!
!  call mpp_get_compute_domain(domain_g2, jsp2, jep2, isp2, iep2)
!  ilenp2 = iep2 - isp2 + 1
!  jlenp2 = jep2 - jsp2 + 1
!
!  call init_transforms(domain_g2,trunc2,nwaves_oe2)
!
!  allocate(sph_wave2(nwaves_oe2,2))
!  allocate(f_wave2(nwaves_oe2,2))
!
!  allocate(topog2(ocny2,ocnx2))
!  allocate(gtopo2(1,jsp2:jep2,isp2:iep2))
!  allocate(stopo2(1,nwaves_oe2,2))
!
!  call get_spherical_wave(sph_wave2)
!  call get_fourier_wave(f_wave2)
!
!  do i = 1, 2
!    do j = 1, nwaves_oe2
!      l = sph_wave2(j,i)
!      m = f_wave1(j,i)
!      j1 = maxloc(sph_wave1(:,i), 1, mask=(sph_wave1(:,i)==l.and.f_wave1(:,i)==m))
!      print *, j1, count(sph_wave1(:,i)==l.and.f_wave1(:,i)==m)
!      stopo2(1,j,i) = stopo1(1,j1,i)
!    end do
!  end do
!
!  call spherical_to_grid(stopo2, grid=gtopo2)
!
!  call write_data("out.nc","topo",gtopo2(1,:,:))

  call fms_io_exit()
  call mpp_exit()

end program main

#else

program main
  use gauss_and_legendre_mod
  use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error, &
                      mpp_npes, mpp_get_current_pelist, mpp_pe, mpp_max, &
                      mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                      mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather, &
                      mpp_declare_pelist, mpp_set_current_pelist, &
                      mpp_get_current_pelist, mpp_sum
  
  use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain, &
                              mpp_global_field, mpp_get_global_domain
  
  use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init, &
                      stdlog, stdout, stderr
  use fms_io_mod, only : fms_io_exit, restart_file_type, register_restart_field, &
          file_exist, field_exist, field_size

  implicit none

  integer :: nlat = 8000, nlat2, trunc=7998
  integer :: num_fourier, num_spherical
  real, allocatable :: sin_hem(:), wts_hem(:)
  real, allocatable :: legendre(:,:,:)
  integer :: fourier_inc=1, jhg(1)=[1]

  call mpp_init()
  call fms_init()

  nlat2 = nlat/2
  num_fourier = trunc
  num_spherical = trunc + 1
  if (mod(num_spherical,2)==0) num_spherical = num_spherical + 1

  allocate(legendre(0:num_fourier, 0:num_spherical, size(jhg)))
  allocate(sin_hem(nlat2),wts_hem(nlat2))
  call compute_gaussian(sin_hem, wts_hem, nlat2)
  call compute_legendre(legendre, num_fourier, fourier_inc,  &
                      num_spherical, sin_hem, nlat2, jhg, .true.)
  call fms_io_exit()
  call mpp_exit()

end program main
#endif
