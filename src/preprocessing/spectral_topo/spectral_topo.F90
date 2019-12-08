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
                             restore_spec_restart, end_transforms, save_wisdom
  
  use ocpack_mod, only : oc_ny, oc_nx, oc_nfour, ocpack_typeP, npack=>oc_npack, get_ocpackP, init_ocpack
  
  implicit none

  integer :: layout(2) = [1,1]
  integer, allocatable :: yextent(:), xextent(:), sph_wave(:,:)
  logical :: reduced=.True., packed=.True.
  integer :: num_lat=94, trunc=62
  type(domain2d) :: domain_g
  integer :: jsp, jep, isp, iep, ilenp, jlenp
  integer :: ocnx, ocny, nwaves_oe

  call mpp_init()
  call fms_init()
 
  if (any(layout==0)) then
      layout = [1,mpp_npes()]
  elseif (layout(1)*layout(2)/=mpp_npes()) then
      call mpp_error(FATAL,'init_atmos: product of layout should be equal to npes')
  endif

  allocate(yextent(layout(1)))
  allocate(xextent(layout(2)))

  call init_ocpack(num_lat, trunc, layout, yextent=yextent, xextent=xextent,&
                  isreduced=reduced, ispacked=packed)

  ocnx = oc_nx()
  ocny = oc_ny()

  if (sum(yextent)/=ocny) then
      print *, 'ocny=', ocny, 'sum(yextent)=', sum(yextent), 'yextent=', yextent
  endif
  if (sum(xextent)/=ocnx) then
      print *, 'ocnx=', ocnx, 'sum(xextent)=', sum(xextent), 'xextent=', xextent
  endif

  call mpp_define_domains([1,ocny,1,ocnx], layout, domain_g, &
                          xextent=yextent, yextent=xextent, kxy=1)

  call mpp_get_compute_domain(domain_g, jsp, jep, isp, iep)
  ilenp = iep - isp + 1
  jlenp = jep - jsp + 1

  call init_transforms(domain_g,trunc,nwaves_oe)

  allocate(sph_wave(nwaves_oe,2))

  call get_spherical_wave(sph_wave)

  print *, sph_wave

  !allocate(topog(ocny,ocnx))

  !allocate(stopo(1,nwaves_oe,2))

end program main

