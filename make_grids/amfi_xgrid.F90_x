program amfi_xgrid

use iso_c_binding

use netcdf
use mpp_mod, only : mpp_init, mpp_error, FATAL, NOTE
use mkxgrid_mod, only : vtx, clipin, poly_area
use ocpack_mod, only : init_ocpack, oc_nx, oc_ny, ocpack_typeP, oc_maxlon, &
        oc_npack, get_ocpackP, oc_nlat

implicit none

interface 
    integer(C_INT) function mppncc(fnm) bind(C,name="mppncc")
        import
        character(kind=C_CHAR) :: fnm(*)
    end function mppncc
end interface

integer, allocatable :: xi(:,:), xj(:,:)
real, allocatable :: xf(:,:)

type(ocpack_typeP), allocatable :: ocP(:,:)
integer :: stat, num_lat, xn
real, parameter :: eps1 = 1e-10

xn = 0

call mpp_init()

num_lat = 94

call make_xgrid_oc2rg(num_lat)

contains 

subroutine make_xgrid_oc2rg(nlat)
    integer, intent(in) :: nlat
    integer(C_INT) :: num_in = 4, n_out
    type(vtx) :: v1(2), v2(2), v_out(4)
    type(vtx), allocatable :: vtxoc(:,:,:), vtxrg(:,:,:)
    real, allocatable :: latb(:), lonfb(:)
    real :: dlonf, dlon, dlat, area, area1
    integer :: nlon, is, ie, i, j, k, n, g, i1, i2
    integer :: ierr, ncid, np_id, xn_id, xi_id, xj_id, xf_id

    if (allocated(xf)) return

    call init_ocpack(nlat,nlat/2,1)
    
    allocate(latb(nlat+1))
    allocate(xi(2,nlat*oc_maxlon()*2))
    allocate(xj(2,nlat*oc_maxlon()*2))
    allocate(xf(2,nlat*oc_maxlon()*2))

    allocate(ocP(oc_npack(),oc_ny()))
    allocate(vtxoc(2,oc_ny(),oc_nx()), vtxrg(2,nlat,oc_maxlon()))
    
    dlat=180./nlat
    
    latb(1) = -90.
    do j = 2, nlat+1
        latb(j) = latb(j-1) + dlat
    end do
    latb(nlat+1)=90.
    
    call get_ocpackP(ocP)
    
    do n = 1, oc_npack()
        do j = 1, oc_ny()
            is = ocP(n,j)%is
            ie = ocP(n,j)%ie
            nlon = ocP(n,j)%ilen
            dlon = 360./nlon
            g = ocP(n,j)%g
            do i = is, ie
                vtxoc(1,j,i)%y = latb(g)
                vtxoc(2,j,i)%y = latb(g+1)
                vtxoc(1,j,i)%x = -dlon/2. + (i-is)*dlon
                vtxoc(2,j,i)%x = dlon/2. + (i-is)*dlon
            end do
        end do
    end do
    
    dlon = 360./oc_maxlon()
    do i = 1, oc_maxlon()
        do j = 1, nlat
            vtxrg(1,j,i)%y = latb(j)        
            vtxrg(2,j,i)%y = latb(j+1)        
            vtxrg(1,j,i)%x = -dlon/2. + (i-1)*dlon
            vtxrg(2,j,i)%x = dlon/2. + (i-1)*dlon
        enddo
    enddo        
    
    xn = 0         
    do n = 1, oc_npack()
        do j = 1, oc_ny()
            is = ocP(n,j)%is
            ie = ocP(n,j)%ie
            g = ocP(n,j)%g
            do i1 = 1, oc_maxlon()
                do i2 = is, ie 
                    n_out = clipin(vtxrg(:,g,i1),vtxoc(:,j,i2),v_out)
                    if (n_out<1) cycle
                    area=poly_area(v_out, n_out)
                    if (area<eps1) cycle
                    xn = xn + 1
                    xj(1,xn) = j
                    xi(1,xn) = i2
                    xj(2,xn) = g
                    xi(2,xn) = i1

                    v_out(1) = vtxoc(1,j,i2)
                    v_out(3) = vtxoc(2,j,i2)
                    v_out(2)%y = v_out(1)%y
                    v_out(2)%x = v_out(3)%x
                    v_out(4)%y = v_out(3)%y
                    v_out(4)%x = v_out(1)%x
                    area1 = poly_area(v_out, 4)
                    xf(1,xn) = area/area1

                    v_out(1) = vtxrg(1,g,i1)
                    v_out(3) = vtxrg(2,g,i1)
                    v_out(2)%y = v_out(1)%y
                    v_out(2)%x = v_out(3)%x
                    v_out(4)%y = v_out(3)%y
                    v_out(4)%x = v_out(1)%x
                    area1 = poly_area(v_out, 4)
                    xf(2,xn) = area/area1
                end do
            end do
        end do
    end do

    ierr = nf90_create('amfi_xgrd.nc',NF90_CLOBBER,ncid)

    ierr = nf90_def_dim(ncid,'np',2,np_id)
    call handle_err(ierr)
    ierr = nf90_def_dim(ncid,'xn',xn,xn_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'xi',NF90_INT,[np_id,xn_id],xi_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'xj',NF90_INT,[np_id,xn_id],xj_id)
    call handle_err(ierr)
    ierr = nf90_def_var(ncid,'xf',NF90_DOUBLE,[np_id,xn_id],xf_id)
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

    ierr = nf90_put_var(ncid,xi_id,xi(:,1:xn)) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,xj_id,xj(:,1:xn)) 
    call handle_err(ierr)
    ierr = nf90_put_var(ncid,xf_id,xf(:,1:xn)) 
    call handle_err(ierr)
    ierr = nf90_close(ncid)
    call handle_err(ierr)

end subroutine make_xgrid_oc2rg

subroutine handle_err(status)
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
        call mpp_error(FATAL,trim(nf90_strerror(status)))
    end if
end subroutine handle_err

end program amfi_xgrid
