program amfi_xgrid

  use iso_c_binding

  use netcdf
  use mpp_mod, only : mpp_init, mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_exit
  use mpp_io_mod, only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_io_init, mpp_get_axes
  use mpp_io_mod, only : mpp_open, mpp_get_info, mpp_get_fields, mpp_get_atts, fieldtype, axistype
  use mpp_io_mod, only : mpp_read, atttype, MPP_OVERWR, mpp_modify_meta, mpp_copy_meta, mpp_write
  use mpp_io_mod, only : mpp_close, mpp_write_meta
  use fms_mod, only : fms_init
  use fms_io_mod, only : read_data, field_size, fms_io_init, fms_io_exit, write_data
  use mkxgrid_mod, only : vtx, clipin, poly_area
  use gauss_and_legendre_mod, only : compute_gaussian
  use constants_mod, only : PI
  
  implicit none
  
  integer, allocatable :: pxi(:), pxj(:)
  integer, allocatable :: rxi(:), rxj(:)
  real, allocatable :: rxf(:)
  real, allocatable :: pxf(:)
  real, allocatable :: xxf(:)
  real, allocatable :: lat(:), lon(:)
  real, allocatable :: latb(:), lonb(:)

  integer :: nlat, nlon, ntime, nvar, natt, ndim
  integer :: ounit, iunit, yid, xid
  type(axistype), allocatable :: oaxes(:)
  type(fieldtype), allocatable :: ofields(:)
  type(fieldtype), allocatable :: fields(:)
  type(axistype) :: time_axis
  character(len=32) :: time_axis_name="NOTIME"
  character(len=128) :: lon_nm='lon', lat_nm='lat'
  
  integer :: stat, xn, i, np=10
  real, parameter :: eps1 = 1e-10
  real, allocatable :: data2d(:,:), rdata2d(:,:)
  real, allocatable :: data3d(:,:,:), rdata3d(:,:,:)
  
  character(len=512) :: file_nm="in.nc", gridf="amfi_grid.nc"
  character(len=128) :: var_nm(50)
  integer, target :: kk, jj, ii
  integer, pointer :: kp, jp, ip

  integer :: ocnx, ocny

  logical :: debug=.false.

  
  namelist/snf_nml/file_nm, var_nm, lon_nm, lat_nm, np, debug
  
  call mpp_init()
  call mpp_io_init()
  call fms_init()
  call fms_io_init()
  
  do i = 1, size(var_nm)
    var_nm(i)="NONAME"
  end do

  read(*,nml=snf_nml)

  if (mpp_pe()==mpp_root_pe()) then
    write(*,nml=snf_nml)
  endif

  call initialize_input_grid()

  call make_xgrid_oc2rg()

  if (open_files_copy_meta()) then
    call process_vars()
  else
    call mpp_error(NOTE,"Nothing to do!!!!!")
  endif 

  call fms_io_exit()
  call mpp_exit()
  
  contains

    subroutine initialize_input_grid()
      real, allocatable :: wts_hem(:), wts_lat(:), sin_hem(:)
      integer :: siz(4), i, j
      real :: dlon, sumwts=0., dlat

      call field_size(trim(file_nm),trim(lon_nm),siz)
      nlon = siz(1)
      call mpp_error(NOTE,"read nlon")

      call field_size(trim(file_nm),trim(lat_nm),siz)
      nlat = siz(1)
      call mpp_error(NOTE,"read nlat")

      allocate(lat(nlat))
      allocate(latb(nlat+1))

      allocate(lon(nlon))
      allocate(lonb(nlon+1))

      call read_data(trim(file_nm),trim(lat_nm),lat)
      call mpp_error(NOTE,"read lat")
      call read_data(trim(file_nm),trim(lon_nm),lon)
      call mpp_error(NOTE,"read lon")

      if (ifgaussian(lat)) then
        call mpp_error(NOTE,"Input data is in gaussian latitudes.")
        if (mod(nlat,2)/=0) call mpp_error(FATAL,"Input nlat is not an even number")
        allocate(sin_hem(nlat/2),wts_hem(nlat/2))
        allocate(wts_lat(nlat))
        call compute_gaussian(sin_hem, wts_hem, nlat/2)
        wts_lat(1:nlat/2) = wts_hem
        wts_lat(nlat:nlat/2+1:-1) = wts_hem

        latb(1) = -0.5*PI
        sumwts = 0.
        do j = 1, nlat-1
            sumwts = sumwts + wts_lat(j)
            latb(j+1) = asin(sumwts-1.)
        end do
        latb(nlat+1) = 0.5*PI
        latb = latb*180./PI
        deallocate(sin_hem, wts_lat, wts_hem)
      else
        call mpp_error(NOTE,"Input data is in regular latitudes.")
        dlat = abs(lat(2)-lat(1))
        latb(1) = -90.
        do j = 1, nlat-1
          latb(j+1) = latb(j) + dlat
        end do
        latb(nlat+1) = 90.
      endif
      if (.not.ifs2n(lat)) latb=latb*-1.

      dlon = abs(lon(2)-lon(1))
      lonb(1) = lon(1) - dlon*0.5
      do j = 1, nlon
        lonb(j+1) = lonb(j) + dlon
      end do

    end subroutine initialize_input_grid
  

    function open_files_copy_meta()
      type(axistype), allocatable :: axes(:)
      type(atttype), allocatable :: gatt(:)
      character(len=128) :: nm, cart
      integer :: n, siz(4), ndimv
      logical :: donothing, open_files_copy_meta

      open_files_copy_meta=.true.

      call mpp_open(iunit,trim(file_nm),action=MPP_RDONLY, &
                    form=MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)

      call mpp_get_info(iunit, ndim, nvar, natt, ntime)

      allocate(fields(nvar))
      allocate(axes(ndim))
      allocate(gatt(natt))

      call mpp_get_axes(iunit, axes, time_axis=time_axis)

      if (ntime > 0) then
        call mpp_get_atts(time_axis, name=time_axis_name)
      endif

      call mpp_get_fields(iunit,fields)

      donothing=.true.
      do n = 1, nvar
        call mpp_get_atts(fields(n),axes=axes,ndim=ndimv)
        if (.not.haslatlon(axes(1:ndimv))) cycle
        donothing=.false.
        exit
      end do

      if (donothing) then
        open_files_copy_meta=.false.
        return
      endif

      call mpp_open(ounit,"o_"//trim(file_nm),action=MPP_OVERWR, &
                    form=MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)

      !do n = 1, natt
      !  call mpp_copy_meta(ounit,gatt(n))
      !end do

      allocate(oaxes(ndim))
      call mpp_get_axes(iunit, oaxes)
      call modify_axes(oaxes)

      do n = 1, ndim
        call mpp_copy_meta(ounit,oaxes(n))
      end do

      call modify_fields(fields)

      do n = 1, ndim
        if(trim(get_name(oaxes(n)))/=trim(time_axis_name)) call mpp_write(ounit,oaxes(n))
      end do

      !call mpp_close(ounit)

    end function open_files_copy_meta


    subroutine process_vars()
      integer :: n, ndimv
      type(axistype) :: axes(ndim)
      character(len=64) :: nm

      do n = 1, nvar
        call mpp_get_atts(fields(n),name=nm,ndim=ndimv,axes=axes)
        if (.not.haslatlon(axes(1:ndimv))) then
        
        else  
          print *, trim(nm)
          call init_rearrange(axes(1:ndimv))
          if (notTXY(axes(1:ndimv))<1) then
            call do_interpolation2d(n,hastime(axes(1:ndimv)))
          endif
        endif
      end do

      call mpp_close(ounit)
      call mpp_close(iunit)

    end subroutine process_vars


    subroutine modify_fields(flds)
      type(fieldtype), intent(in) :: flds(:)
      integer :: i, j, ndimv, ndimv1, idk
      type(axistype) :: axes(ndim), axes1(ndim)
      character(len=128) :: nm, units, longname
      real :: minv, maxv, scalev, add, missing

      print *, 'ndim=', ndim

      allocate(ofields(nvar))
      do i = 1, nvar
        call mpp_get_atts(flds(i), name=nm, units=units, longname=longname, min=minv, &
          max=maxv, missing=missing, ndim=ndimv, axes=axes, scale=scalev, add=add)


        if (.not.haslatlon(axes(1:ndimv))) then
          ofields(i) = flds(i)
          cycle
        end if

        if (hastime(axes(1:ndimv))) then
          axes1(ndimv) = time_axis
        endif

        idk=notTXY(axes(1:ndimv))

        if (idk>0) then
          axes1(1) = axes(idk)
          axes1(2) = oaxes(yid)
          axes1(3) = oaxes(xid)
        else
          axes1(1) = oaxes(yid)
          axes1(2) = oaxes(xid)
        endif

        call mpp_write_meta(ounit, ofields(i), axes=axes1(1:ndimv), name=nm, units=units, &
          longname=longname, standard_name=longname, min=minv, max=maxv, missing=missing, &
          scale=scalev, add=add)

     end do

    end subroutine modify_fields


    function notTXY(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: notTXY, i

      notTXY=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lat_nm) .or. &
            trim(get_name(axes(i)))==trim(lon_nm) .or. &
            trim(get_name(axes(i)))==trim(time_axis_name)) cycle
        notTXY=i
        exit
      end do
    end function notTXY

    function id_Y(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: id_Y, i

      id_Y=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lat_nm)) then
          id_Y=i
          exit
        endif
      end do
    end function id_Y

    function id_X(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: id_X, i

      id_X=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lon_nm)) then
          id_X=i
          exit
        endif
      end do
    end function id_X


    subroutine modify_axes(axes)
      type(axistype), intent(inout) :: axes(:)
      real, allocatable :: tmp(:)

      integer :: i, j

      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lat_nm)) then
          allocate(tmp(ocny))
          tmp=[(j, j=1, ocny)]
          call mpp_modify_meta(axes(i), 'y', 'nounit', 'y', 'Y', tmp)
          yid=i
        elseif (trim(get_name(axes(i)))==trim(lon_nm)) then
          allocate(tmp(ocnx))
          tmp=[(j, j=1, ocnx)]
          call mpp_modify_meta(axes(i), 'x', 'nounit', 'x', 'X', tmp)
          xid=i
        end if
        if (allocated(tmp)) deallocate(tmp)
      end do

    end subroutine modify_axes


    subroutine do_interpolation2d(ni,istime)
      integer, intent(in) :: ni
      logical, intent(in) :: istime
      integer :: siz(4), n

      call mpp_get_atts(fields(ni),siz=siz)

      if (allocated(data2d)) deallocate(data2d)

      allocate(data2d(siz(1),siz(2)))
      if (.not.allocated(rdata2d)) allocate(rdata2d(ocny,ocnx))

      if (istime) then
        call mpp_read(iunit, fields(ni), data2d, 1)
        rdata2d = 0.
        do n = 1, xn
          ip = rxi(n)
          jp = rxj(n)
          if (debug) print *, ii, jj, kk
          rdata2d(pxj(n),pxi(n)) = rdata2d(pxj(n),pxi(n)) + &
            data2d(ii,jj) * pxf(n)
        end do
        call mpp_write(ounit,ofields(ni),rdata2d,1.)
      else
        call mpp_read(iunit, fields(ni), data2d)
        rdata2d = 0.
        do n = 1, xn
          ip = rxi(n)
          jp = rxj(n)
          rdata2d(pxj(n),pxi(n)) = rdata2d(pxj(n),pxi(n)) + &
            data2d(ii,jj) * pxf(n)
        end do
        call mpp_write(ounit,ofields(ni),rdata2d)
      end if

    end subroutine do_interpolation2d



    subroutine init_rearrange(axes)
      type(axistype), intent(in) :: axes(:)
      integer :: i, ndim
      character(len=128) :: nm

      if (trim(get_name(axes(1)))==trim(lon_nm)) then
        ip => ii
        if (trim(get_name(axes(2)))==trim(lat_nm)) then
          jp => jj
          kp => kk
        else
          jp => kk
          kp => jj
        endif
      elseif (trim(get_name(axes(2)))==trim(lon_nm)) then
        ip => jj
        if (trim(get_name(axes(1)))==trim(lat_nm)) then
          jp => ii
          kp => kk
        else
          jp => kk
          kp => ii
        endif
      elseif (trim(get_name(axes(3)))==trim(lon_nm)) then
        ip => kk
        if (trim(get_name(axes(1)))==trim(lat_nm)) then
          jp = ii
          kp = jj
        else
          jp = jj
          kp = ii
        endif 
      endif

    end subroutine init_rearrange

    function get_name(axes)
      type(axistype), intent(in) :: axes
      character(len=128) :: get_name

      call mpp_get_atts(axes,name=get_name)
      return
    end function get_name

    function hastime(axes)
      type(axistype) :: axes(:)
      logical :: hastime
      integer :: i
      character(len=128) :: nm
      
      hastime = .false.

      do i = 1, size(axes,1)
        call mpp_get_atts(axes(i), name=nm)
        if (trim(nm)==trim(time_axis_name)) then
          hastime = .true.
          return
        endif
      end do
    end function hastime

    function haslatlon(axes)
      type(axistype) :: axes(:)
      logical :: haslatlon, haslat, haslon
      integer :: i
      character(len=128) :: nm

      haslat=.false.
      haslon=.false.
      do i = 1, size(axes,1)
        call mpp_get_atts(axes(i),name=nm)
        if ((.not.haslat) .and. trim(nm)==trim(lat_nm)) haslat=.true.
        if ((.not.haslon) .and. trim(nm)==trim(lon_nm)) haslon=.true.
      end do

      haslatlon = haslat.and.haslon

    end function haslatlon

    function ifs2n(x)
      real, intent(in) :: x(:)
      logical :: ifs2n

      ifs2n = .true.

      ifs2n=(x(2)-x(1))>0.
      return

    end function ifs2n

    function ifgaussian(x)
      real, intent(in) :: x(:)
      logical :: ifgaussian

      real :: dx1, dx2

      ifgaussian = .false.

      dx1 = abs(x(2)-x(1))
      dx2 = abs(x(3)-x(2))

      ifgaussian=abs(dx1-dx2)>eps1
      return

    end function ifgaussian

    subroutine make_xgrid_oc2rg()
      integer(C_INT) :: num_in = 4, n_out
      type(vtx) :: v1(2), v2(2), v_out(4)
      type(vtx), allocatable :: vtxoc(:,:,:), vtxrg(:,:,:)
      real, allocatable :: lonbw(:,:), lonbe(:,:), latbn(:,:), latbs(:,:)
      integer, allocatable :: xi(:,:), xj(:,:)
      real, allocatable :: xf(:,:)
    
      real :: area, area1
      integer :: is, ie, i, j, k, n, g, i1, i2, xnsize, siz(4)
    
      xnsize=size(lonb,1)*size(latb,1)*np
    
      call field_size(gridf,'latbs',siz)
      xnsize = max(xnsize,siz(1)*siz(2)*np)
      print *, 'siz=',siz(:)
   
      ocny=siz(1)
      ocnx=siz(2)

      allocate(latbs(siz(1),siz(2)))
      allocate(latbn(siz(1),siz(2)))
      allocate(lonbe(siz(1),siz(2)))
      allocate(lonbw(siz(1),siz(2)))
      allocate(vtxoc(2,siz(1),siz(2)))
    
      allocate(vtxrg(2,size(latb,1)-1,size(lonb,1)-1))
    
      allocate(xi(2,xnsize))
      allocate(xj(2,xnsize))
      allocate(xf(3,xnsize))

      call read_data(gridf,'latbs',latbs)
      call read_data(gridf,'latbn',latbn)
      call read_data(gridf,'lonbw',lonbw)
      call read_data(gridf,'lonbe',lonbe)
    
      do j = 1, size(latbs,1)
        do i = 1, size(latbs,2)
          vtxoc(1,j,i)%y = latbs(j,i)
          vtxoc(2,j,i)%y = latbn(j,i)
          vtxoc(1,j,i)%x = lonbw(j,i)
          vtxoc(2,j,i)%x = lonbe(j,i) 
        end do
      end do
      
      do j = 1, size(latb,1)-1
        do i = 1, size(lonb,1)-1
          vtxrg(1,j,i)%y = latb(j)        
          vtxrg(2,j,i)%y = latb(j+1)        
          vtxrg(1,j,i)%x = lonb(i)
          vtxrg(2,j,i)%x = lonb(i+1)
        enddo
      enddo        
      
      xn = 0         
      do j = 1, size(latbs,1)
        do i2 = 1, size(latbs,2)
          do g = 1, size(latb,1)-1
            do i1 = 1, size(lonb,1)-1
              n_out = clipin(vtxrg(:,g,i1),vtxoc(:,j,i2),v_out)
              if (n_out<1) cycle
              area=poly_area(v_out, n_out)
              if (area<eps1) cycle
              xn = xn + 1
              if (xn > xnsize) then
                print *, xn, xnsize
                call mpp_error(FATAL,"xn > xnsize")
              endif
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
              xf(3,xn) = area
            end do
          end do
        end do
      end do
    
      deallocate(latbs)
      deallocate(latbn)
      deallocate(lonbe)
      deallocate(lonbw)
      deallocate(vtxoc)
      deallocate(vtxrg)
    
      allocate(pxi(xn))
      allocate(pxj(xn))
      allocate(pxf(xn))
      allocate(rxi(xn))
      allocate(rxj(xn))
      allocate(rxf(xn))
      allocate(xxf(xn))
    
      pxi=xi(1,1:xn)
      pxj=xj(1,1:xn)
      pxf=xf(1,1:xn)
      rxi=xi(2,1:xn)
      rxj=xj(2,1:xn)
      rxf=xf(2,1:xn)
      xxf=xf(3,1:xn)

      print *, maxval(pxi), maxval(pxj), maxval(rxi), maxval(rxj)
    
      deallocate(xi)
      deallocate(xj)
      deallocate(xf)
    
    end subroutine make_xgrid_oc2rg

    subroutine handle_err(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
        call mpp_error(FATAL,trim(nf90_strerror(status)))
      end if
    end subroutine handle_err

end program amfi_xgrid
