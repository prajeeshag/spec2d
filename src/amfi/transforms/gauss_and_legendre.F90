module gauss_and_legendre_mod

use       fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                         write_version_number, NOTE, read_data, write_data

use constants_mod, only: pi

use netcdf

!-----------------------------------------------------------------------
!   computes Gaussian latitudes and associated Legendre polynomials
!
!-----------------------------------------------------------------------

implicit none
private

character(len=128), parameter :: version = '$Id: gauss_and_legendre.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=128), parameter :: tagname ='$Name: tikal $'

logical :: entry_to_logfile_done=.false.

character (len=512) :: fnm

public :: compute_legendre, compute_gaussian


contains

!-----------------------------------------------------------------------
subroutine compute_legendre(legendre, num_fourier, fourier_inc,  &
                num_spherical, sin_lat, n_lat, jhg, write_coef)
!-----------------------------------------------------------------------

integer, intent (in) :: num_fourier, fourier_inc, num_spherical, n_lat
real,    intent (in), dimension(n_lat) :: sin_lat
integer, intent(in) :: jhg(:)
logical, optional :: write_coef

!real, intent (out), dimension(0:num_fourier, 0:num_spherical,n_lat)  :: legendre
real, intent (out), dimension(0:, 0:, :)  :: legendre

integer :: j, m, n , fourier_max, j1, j2, ios
real, dimension(0:num_fourier*fourier_inc, 0:num_spherical) :: poly, eps, m2, l2
real, dimension(0:num_fourier*fourier_inc) :: b
real, dimension(n_lat) :: cos_lat
logical :: file_exist

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done=.true.
endif

fourier_max = num_fourier*fourier_inc

write(fnm,'("pnm_",I4.4,"_",I4.4,"_",I4.4,"_",I4.4,".nc")') &
        num_fourier, num_spherical, n_lat, 1
inquire(file='legendre/'//trim(fnm),exist=file_exist)

write(fnm,'("pnm_",I4.4,"_",I4.4,"_",I4.4)') &
        num_fourier, num_spherical, n_lat

if (file_exist) then
  call error_mesg('compute_legendre','Reading from Legendre poly file legendre/'//trim(fnm), NOTE)
  do j = 1, size(jhg)
    j1 = jhg(j)
    call read_legendre(legendre(:,:,j),j1)
    !write(fnm,'("pnm_",I4.4,"_",I4.4,"_",I4.4,"_",I4.4,".nc")') &
    !          num_fourier, num_spherical, n_lat, j1
    !open(10,file='legendre/'//trim(fnm),form='unformatted',status='old',action='read',iostat=ios)
    !if (ios/=0) call error_mesg('compute_legendre','Error opening file legendre/'//trim(fnm), FATAL)
    !read(10,iostat=ios) legendre(:,:,j)
    !if (ios/=0) call error_mesg('compute_legendre','Error reading file legendre/'//trim(fnm), FATAL)
    !close(10)
  end do
return
endif

do j=1,n_lat
  cos_lat(j)=sqrt(1-sin_lat(j)*sin_lat(j))
end do

do n=0,num_spherical
  do m=0,fourier_max
    m2(m,n)   = m**2
    l2(m,n)   = (m+n)**2
  end do
end do

eps= sqrt((l2 -m2)/(4.0*l2 - 1.0))

do m=1,fourier_max
  b(m) = SQRT(0.5*(2.0*float(m) + 1.0)/float(m))
end do

poly(0,0) = sqrt(0.5)
j2=0
do j=1,n_lat

  if (mpp_pe()==mpp_root_pe()) print *, "compute_legendre j=", j

  do m=1,fourier_max
    poly(m,0) = b(m)*cos_lat(j)*poly(m-1,0)
  end do

  do m=0,fourier_max
     poly(m,1) =  sin_lat(j)*poly(m,0)/eps(m,1)
  end do

  do n=2,num_spherical
    do m=0,fourier_max
      poly(m,n) = (sin_lat(j)*poly(m,n-1) - eps(m,n-1)*poly(m,n-2))/eps(m,n)
    end do
  end do

  
  if (present(write_coef) .and. write_coef) then
    !write(fnm,'("pnm_",I4.4,"_",I4.4,"_",I4.4,"_",I4.4,".nc")') &
    !    num_fourier, num_spherical, n_lat, j
    if (mpp_pe()==mpp_root_pe()) then
      !open(10,file=trim(fnm),form='unformatted',iostat=ios)
      !if (ios/=0) call error_mesg('compute_legendre','Error opening file '//trim(fnm), FATAL)
      !write(10,iostat=ios) poly
      !if (ios/=0) call error_mesg('compute_legendre','Error writting file '//trim(fnm), FATAL)
      !close(10)
      call write_legendre(poly,j)
    endif
  else
    if (any(jhg==j)) then
      j1 = maxloc(jhg,1,jhg==j)
      do n = 0, num_spherical
        do m = 0,num_fourier
          legendre(m,n,j1) = poly(m*fourier_inc,n)
        end do
      end do
      j2 = j2 + 1
      if (j2==size(jhg)) return
    endif
  endif

end do

return
end subroutine compute_legendre

!----------------------------------------------------------------------
subroutine compute_gaussian(sin_hem, wts_hem, n_hem)
!----------------------------------------------------------------------
!
!     reference:
!       press, h. william, et. al., numerical recipes (fortran version),
!       cambridge, england: cambridge university press (1990)
! 
!------------------------------------------------------------------------

integer, intent (in) :: n_hem
real, intent (out), dimension(n_hem) :: sin_hem, wts_hem

real :: converg
integer :: itermax
integer :: i, iter, j, n, nprec
real :: pp, p1, p2, p3, z, z1

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done=.true.
endif

! must use a more relaxed convergence criteria on the
! workstations than that for the cray T90
! fez code is commented out

!if(kind(converg).eq.8) then
!  converg = 1.0E-15
!else if(kind(converg).eq.4) then
!  converg = 1.0E-7
!else
!  call error_mesg('compute_gaussian','dont know what value to use for converg', FATAL)
!end if

! The 2 lines of code below will yeild a different result than the fez code
! when kind(converg)=4. converg is 1.0E-6 instead of 1.0E-7
! This should be investigated further, but it's OK for now because it yeilds
! the same result on the HPCS. -- pjp
nprec = precision(converg)
converg = .1**nprec


itermax = 10

n=2*n_hem
do i=1,n_hem
  z = cos(pi*(i - 0.25)/(n + 0.5))
  do iter=1,itermax
     p1 = 1.0
     p2 = 0.0

     do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
     end do

     pp = n*(z*p1 - p2)/(z*z - 1.0E+00)
     z1 = z
     z  = z1 - p1/pp
     if(ABS(z - z1) .LT. converg) go to 10
  end do
  call error_mesg('compute_gaussian','abscissas failed to converge in itermax iterations', FATAL)

  10  continue

  sin_hem (i)     = z
  wts_hem (i)     = 2.0/((1.0 - z*z)*pp*pp)

end do

return
end subroutine compute_gaussian
!----------------------------------------------------------------------


subroutine write_legendre(pnm,l)

  real, intent(in) :: pnm(:,:)
  integer, intent(in) :: l
  character(len=512) :: filename
  integer :: stat, ncid, xid, yid, varid

  write(filename,'(A,"_",I4.4,".nc")') trim(fnm), l

  stat = nf90_create(trim(filename),OR(nf90_clobber,nf90_netcdf4),ncid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_def_dim(ncid,'x',size(pnm,1),xid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_def_dim(ncid,'y',size(pnm,2),yid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_def_var(ncid,'pnm',nf90_double,[xid,yid],varid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_def_var_deflate(ncid,varid,0,0,0)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_enddef(ncid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)
  
  stat = nf90_put_var(ncid,varid,pnm)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_close(ncid)
  if (stat/=nf90_noerr) call error_mesg('write_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

end subroutine write_legendre



subroutine read_legendre(pnm,l)

  real, intent(out) :: pnm(:,:)
  integer, intent(in) :: l
  character(len=512) :: filename
  integer :: stat, ncid, xid, yid, varid

  write(filename,'(A,"_",I4.4,".nc")') 'legendre/'//trim(fnm), l

  stat = nf90_open(trim(filename),0,ncid)
  if (stat/=nf90_noerr) call error_mesg('read_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_inq_varid(ncid,'pnm',varid)
  if (stat/=nf90_noerr) call error_mesg('read_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_get_var(ncid,varid,pnm)
  if (stat/=nf90_noerr) call error_mesg('read_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

  stat = nf90_close(ncid)
  if (stat/=nf90_noerr) call error_mesg('read_legendre', &
    trim(filename)//": "//trim(nf90_strerror(stat)), FATAL)

end subroutine read_legendre


!subroutine handle_err(status)
!   integer, intent ( in) :: status
!    if(status /= nf90_noerr) then
!       print *, trim(nf90_strerror(status))
!    stop “Stopped”
!         end if
!end subroutine handle_err

end module gauss_and_legendre_mod


