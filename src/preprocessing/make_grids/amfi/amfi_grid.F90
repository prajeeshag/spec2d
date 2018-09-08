program amfi_grid

use netcdf
use mpp_mod, only : mpp_init, mpp_error, FATAL, NOTE
use ocpack_mod, only : init_ocpack, oc_nx, oc_ny, ocpack_typeP, oc_maxlon, &
        oc_npack, get_ocpackP, oc_nlat
use gauss_and_legendre_mod, only : compute_gaussian
use constants_mod, only : PI

implicit none

real, allocatable :: lonbe(:,:), lonbw(:,:), lonc(:,:)
real, allocatable :: latbs(:,:), latbn(:,:), latc(:,:)

type(ocpack_typeP), allocatable :: ocP(:,:)
integer :: stat, num_lat, xn
real, parameter :: eps1 = 1e-10

call mpp_init()

num_lat = 0

write(*,*) 'Enter num_lat: '
read(*,*) num_lat

if (num_lat<4.or.mod(num_lat,2)/=0) call mpp_error(FATAL,'num_lat should be >=4, and a even number')

call init_ocpack(num_lat,num_lat/2,[1,1])

call make_grid()

contains 

subroutine make_grid()
    real :: dlonf, dlon, dlat, area, area1, sumwts
    integer :: nlon, is, ie, i, j, k, n, g, i1, i2
    integer :: ocnx, ocny, nlat, ierr, ncid, ocnx_id, ocny_id, lonbw_id, &
               lonbe_id, lonc_id, latbs_id, latbn_id, latc_id
    real, allocatable :: sin_hem(:), wts_hem(:), wts_lat(:), sin_lat(:), latb(:)

    ocnx = oc_nx(); ocny = oc_ny(); nlat = oc_nlat()

    allocate(lonbe(ocny,ocnx), lonbw(ocny,ocnx), &
             latbs(ocny,ocnx), latbn(ocny,ocnx), &
             lonc(ocny,ocnx), latc(ocny,ocnx))

    allocate(latb(nlat+1))

    allocate(sin_hem(nlat/2),wts_hem(nlat/2),wts_lat(nlat),sin_lat(nlat))
 
    allocate(ocP(oc_npack(),oc_ny()))

    call get_ocpackP(ocP)

    call compute_gaussian(sin_hem, wts_hem, nlat/2)

    wts_lat(1:nlat/2) = wts_hem
    wts_lat(nlat:nlat/2+1:-1) = wts_hem

    sin_lat(1:nlat/2) = -sin_hem
    sin_lat(nlat:nlat/2+1:-1) = sin_hem

    latb(1) = -0.5*PI
    sumwts = 0.
    do j = 1, nlat-1
        sumwts = sumwts + wts_lat(j)
        latb(j+1) = asin(sumwts-1.)
    end do
    latb(nlat+1) = 0.5*PI

    latb = latb*180./PI
    
    do n = 1, oc_npack()
        do j = 1, oc_ny()
            is = ocP(n,j)%is
            ie = ocP(n,j)%ie
            nlon = ocP(n,j)%ilen
            dlon = 360./nlon
            g = ocP(n,j)%g
            do i = is, ie
                lonbw(j,i) = -dlon/2. + (i-is)*dlon
                lonbe(j,i) = dlon/2. + (i-is)*dlon
                lonc(j,i) = (i-is)*dlon
                latbs(j,i) = latb(g)
                latbn(j,i) = latb(g+1)
                latc(j,i) = asin(sin_lat(g))*180./PI
            end do
        end do
    end do
    
    ierr = nf90_create('amfi_grid.nc',NF90_CLOBBER,ncid)

    ierr = nf90_def_dim(ncid,'ocnx',ocnx,ocnx_id)
    call handle_err(ierr)
    ierr = nf90_def_dim(ncid,'ocny',ocny,ocny_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'lonbw',NF90_DOUBLE,[ocny_id,ocnx_id],lonbw_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'lonbe',NF90_DOUBLE,[ocny_id,ocnx_id],lonbe_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'lonc',NF90_DOUBLE,[ocny_id,ocnx_id],lonc_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'latbs',NF90_DOUBLE,[ocny_id,ocnx_id],latbs_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'latbn',NF90_DOUBLE,[ocny_id,ocnx_id],latbn_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'latc',NF90_DOUBLE,[ocny_id,ocnx_id],latc_id)
    call handle_err(ierr)

    ierr = nf90_put_att(ncid,NF90_GLOBAL,'maxlon',oc_maxlon())
    call handle_err(ierr)
    ierr = nf90_put_att(ncid,NF90_GLOBAL,'ocny',oc_ny())
    call handle_err(ierr)
    ierr = nf90_put_att(ncid,NF90_GLOBAL,'ocnx',oc_nx())
    call handle_err(ierr)
    ierr = nf90_put_att(ncid,NF90_GLOBAL,'nlat',oc_nlat())
    call handle_err(ierr)
    ierr = nf90_enddef(ncid)
    call handle_err(ierr)

    ierr = nf90_put_var(ncid,lonbw_id,lonbw) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,lonbe_id,lonbe) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,lonc_id,lonc) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,latbs_id,latbs) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,latbn_id,latbn) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,latc_id,latc) 
    call handle_err(ierr)

    ierr = nf90_close(ncid)
    call handle_err(ierr)

end subroutine make_grid

subroutine handle_err(status)
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
        call mpp_error(FATAL,trim(nf90_strerror(status)))
    end if
end subroutine handle_err

end program amfi_grid
