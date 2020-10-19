program main
  use topog_regularization_mod
  use constants_mod, only : RADIUS
  use transforms_mod, only : transforms_init 
  use fms_io_mod, only : field_size, read_data, fms_io_init, write_data, &
                         fms_io_exit
  use fms_mod, only : fms_init, fms_end

  implicit none

  integer, parameter :: trunc=62
  integer :: lat_max, num_lon
  real, allocatable :: topo_uf(:,:), topo_re(:,:)
  real :: lambda, actual_fraction_smoothed
  character(len=256) :: fname="topo.nc",fldnm="topo"
  integer :: siz(4)
  real :: fact=0.95

  call fms_init()
  call fms_io_init()
  
  write(*,*) "Enter : fname, fldnm, fact"
  read(*,*) fname, fldnm, fact

  call field_size(trim(fname),trim(fldnm),siz)

  print *, siz

  num_lon = siz(1)
  lat_max = siz(2)

  allocate(topo_uf(num_lon,lat_max))
  allocate(topo_re(num_lon,lat_max))

  topo_uf = 0.

  call read_data(trim(fname),trim(fldnm),topo_uf)

  where(topo_uf<0.) topo_uf=0.

  call transforms_init(radius, lat_max, num_lon, trunc, 1, trunc, .false.)
  call compute_lambda(fact, topo_uf<=0., topo_uf, lambda, actual_fraction_smoothed) 
  call regularize(lambda, topo_uf<=0., topo_uf, topo_re, actual_fraction_smoothed)
  
  call write_data('topo_re', 'topo_re', topo_re) 

  call fms_io_exit()
  call fms_end()
  

end program main
