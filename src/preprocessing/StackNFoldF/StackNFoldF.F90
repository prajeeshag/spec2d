program amfi_xgrid

  use iso_c_binding

  use netcdf
  use mpp_mod, only : mpp_init, mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_exit, mpp_npes, mpp_broadcast
  use mpp_mod, only : mpp_set_current_pelist, mpp_declare_pelist, WARNING
  use mpp_io_mod, only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_io_init, mpp_get_axes
  use mpp_io_mod, only : mpp_open, mpp_get_info, mpp_get_fields, mpp_get_atts, fieldtype, axistype
  use mpp_io_mod, only : mpp_read, atttype, MPP_OVERWR, mpp_modify_meta, mpp_copy_meta, mpp_write
  use mpp_io_mod, only : mpp_close, mpp_write_meta, mpp_get_times, mpp_get_axis_data
  use fms_io_mod, only : field_size, read_data, fms_io_init, fms_io_exit, fms_io_init
  use fms_io_mod, only : open_namelist_file, close_file
  use mkxgrid_mod, only : vtx, clipin, poly_area
  use gauss_and_legendre_mod, only : compute_gaussian
  use constants_mod, only : PI
  
  implicit none
  
  character(len=128), parameter :: gridf="amfi_grid.nc", gridf1="grid_spec.nc"
  character(len=13), dimension(6) :: lat_units=['degrees_north',& 
                                                'degree_north ',& 
                                                'degree_N     ',& 
                                                'degrees_N    ',& 
                                                'degreeN      ',& 
                                                'degreesN     ']

  character(len=13), dimension(6) :: lon_units=['degrees_east ',& 
                                                'degree_east  ',& 
                                                'degree_E     ',& 
                                                'degrees_E    ',& 
                                                'degreeE      ',& 
                                                'degreesE     ']

  character(len=10), dimension(2) :: lat_names=['latitude  ', &
                                                'lat       ']

  character(len=10), dimension(2) :: lon_names=['longitude ', &
                                                'lon       ']

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
  real, allocatable :: times(:)
  character(len=32) :: time_axis_name="NOTIME"
  character(len=128) :: lon_nm='NONAME', lat_nm='NONAME'
  
  integer :: stat, xn, i, np=10
  real, parameter :: eps1 = 1e-10
  logical, allocatable, dimension(:,:) :: lmask
  integer, allocatable, dimension(:,:) :: max_area_loc
  
  character(len=512) :: in_file="in.nc"
  character(len=512) :: out_file="snf_out.nc"

  character(len=128) :: var_nm(50)
  integer, target :: kk, jj, ii
  integer, pointer :: kp, jp, ip

  integer :: ocnx, ocny, search_frac=4
  integer, allocatable :: pes(:)

  logical :: debug=.false., typedata=.false., landdata=.false., valid_pe=.true.
  integer :: unit

  
  namelist/snf_nml/in_file, out_file, lon_nm, lat_nm, np, debug, typedata, landdata, search_frac
  
  call mpp_init()
  call mpp_io_init()
  call fms_io_init()
  
  do i = 1, size(var_nm)
    var_nm(i)="NONAME"
  end do

  if (mpp_pe()==mpp_root_pe()) then
    read(*,nml=snf_nml)
    write(*,nml=snf_nml)
  endif

  call broadcast_char(in_file)
  call broadcast_char(out_file)
  call broadcast_char(lon_nm)
  call broadcast_char(lat_nm)

  call mpp_broadcast(np,mpp_root_pe())
  call mpp_broadcast(debug,mpp_root_pe())
  call mpp_broadcast(typedata,mpp_root_pe())
  call mpp_broadcast(landdata,mpp_root_pe())
  call mpp_broadcast(search_frac,mpp_root_pe())

  call initialize_input_grid()

  if (valid_pe) then
    call make_xgrid_oc2rg()

    if (open_files_copy_meta()) then
      call process_vars()
    else
      call mpp_error(WARNING,"Nothing to do!!!!!")
    endif
  else
      call mpp_error(WARNING,"Nothing to do for this pe!!!!!")
  endif

  call fms_io_exit()
  call mpp_exit()
  
  contains

    subroutine initialize_input_grid()
      real, allocatable :: wts_hem(:), wts_lat(:), sin_hem(:)
      integer :: siz(4), i, j
      real :: dlon, sumwts=0., dlat
      type(axistype), allocatable :: axes(:)
      integer :: ntime1

      call mpp_open(iunit,trim(in_file),action=MPP_RDONLY, &
                    form=MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)

      call mpp_get_info(iunit, ndim, nvar, natt, ntime)

      ntime1 = max(ntime,1)

      if (ntime1<mpp_npes()) then
        allocate(pes(ntime1))
        do i = 0, ntime1-1
          pes(i+1) = i
        end do

        call mpp_declare_pelist(pes)
        
        valid_pe = .false.
        if (any(mpp_pe()==pes)) then
           call mpp_set_current_pelist(pes,no_sync=.true.)
           valid_pe=.true.
        endif
      else
        valid_pe = .true.
        allocate(pes(mpp_npes()))
        do i = 0, mpp_npes()-1
          pes(i+1) = i
        end do
      end if

      if (.not.valid_pe) return

      allocate(fields(nvar))
      allocate(axes(ndim))

      call mpp_get_axes(iunit, axes, time_axis=time_axis)

      call set_lonlat(axes)

      allocate(latb(nlat+1))
      allocate(lonb(nlon+1))

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


    subroutine set_lonlat(axes)
      type(axistype), intent(in) :: axes(:)
      
      integer :: i, siz
      character(len=13) :: units
      character(len=128) :: name

      if (trim(lon_nm)=='NONAME') then
        do i = 1, size(axes,1)
          call mpp_get_atts(axes(i),units=units,name=name,len=siz)
          units = trim(adjustl(units))
          name = trim(adjustl(name))
          if (any(units==lon_units) .or. any(name(1:10)==lon_names) ) then
            lon_nm=name
            nlon=siz
            if (allocated(lon)) call mpp_error(FATAL,"Cannot have more than one longitude axis in file")
            allocate(lon(nlon))
            call mpp_get_axis_data(axes(i),lon)
          endif
        end do
      else
        do i = 1, size(axes,1)
          call mpp_get_atts(axes(i),units=units,name=name,len=siz)
          name = trim(adjustl(name))
          if (trim(name)==trim(lon_nm)) then
            nlon=siz
            if (allocated(lon)) call mpp_error(FATAL,"Cannot have more than one longitude axis in file")
            allocate(lon(nlon))
            call mpp_get_axis_data(axes(i),lon)
            exit
          endif
        end do
      end if

      if (trim(lon_nm)=='NONAME') call mpp_error(FATAL,"Could'nt find longitude")

      if (trim(lat_nm)=='NONAME') then
        do i = 1, size(axes,1)
          call mpp_get_atts(axes(i),units=units,name=name,len=siz)
          units = trim(adjustl(units))
          name = trim(adjustl(name))
          if (any(units==lat_units) .or. any(name(1:10)==lat_names) ) then
            lat_nm=name
            nlat=siz
            if (allocated(lat)) call mpp_error(FATAL,"Cannot have more than one latitude axis in file")
            allocate(lat(nlat))
            call mpp_get_axis_data(axes(i),lat)
          endif
        end do
      else
        do i = 1, size(axes,1)
          call mpp_get_atts(axes(i),units=units,name=name,len=siz)
          name = trim(adjustl(name))
          if (trim(name)==trim(lat_nm)) then
            nlat=siz
            if (allocated(lat)) call mpp_error(FATAL,"Cannot have more than one latitude axis in file")
            allocate(lat(nlat))
            call mpp_get_axis_data(axes(i),lat)
            exit
          endif
        end do
      end if

      if (trim(lat_nm)=='NONAME') call mpp_error(FATAL,"Could'nt find latitude")

    end subroutine set_lonlat
  

    function open_files_copy_meta()
      type(axistype) :: axes(ndim)
      type(atttype), allocatable :: gatt(:)
      character(len=128) :: nm, cart
      integer :: n, siz(4), ndimv
      logical :: donothing, open_files_copy_meta

      open_files_copy_meta=.true.

      if (ntime > 0) then
        call mpp_get_atts(time_axis, name=time_axis_name)
        allocate(times(ntime))
        call mpp_get_times(iunit, times)
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

      call mpp_open(ounit,trim(out_file),action=MPP_OVERWR, &
                    form=MPP_NETCDF,threading=MPP_MULTI,&
                    nc4parallel=.true.)

      call mpp_write_meta(ounit,'in_file',cval=trim(in_file))
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
      integer :: id_k, t, ti, id_x, id_y
      integer :: tstart, tstep, tend, siz(4), siz_r(4)
      logical :: hstime, hslatlon, istype
      character(len=32) :: units

      tstart = mpp_pe()+1
      tstep = mpp_npes()
      tend = ntime + (tstep-mod(ntime,tstep))

      if (mpp_pe()==mpp_root_pe()) then
        do n = 1, nvar
          call mpp_get_atts(fields(n), name=nm, ndim=ndimv, axes=axes, &
            units=units, siz=siz)

          istype = typedata.or.trim(adjustl(units))=='types'

          hstime = hastime(axes(1:ndimv))

          if (hstime) cycle

          if (.not.haslatlon(axes(1:ndimv))) then
          
          else  
            call mpp_error(NOTE,'Interpolating .... '//trim(nm))
            call init_rearrange(axes(1:ndimv))

            id_k=id_Kaxis(axes(1:ndimv))
            id_x=id_Xaxis(axes(1:ndimv))
            id_y=id_Yaxis(axes(1:ndimv))

            siz_r(1)=siz(id_x)
            siz_r(2)=siz(id_y)

            ti = -1

            if (id_k<1) then
              call do_interpolation2d(fields(n), ofields(n), ti, siz, siz_r)
            else
              siz_r(3)=siz(id_k)
              call do_interpolation3d(fields(n), ofields(n), ti, siz, siz_r)
            endif

          endif
        end do
      endif

      do t = tstart, tend, tstep
        ti = t
        if (ti > ntime) ti=ti-tstep
        do n = 1, nvar
          call mpp_get_atts(fields(n), name=nm, ndim=ndimv, axes=axes, &
            units=units, siz=siz)

          istype = typedata.or.trim(adjustl(units))=='types'

          hstime = hastime(axes(1:ndimv))

          if (.not.hstime) cycle

          if (.not.haslatlon(axes(1:ndimv))) then
          
          else  
            call mpp_error(NOTE,'Interpolating .... '//trim(nm))
            call init_rearrange(axes(1:ndimv))

            id_k=id_Kaxis(axes(1:ndimv))
            id_x=id_Xaxis(axes(1:ndimv))
            id_y=id_Yaxis(axes(1:ndimv))

            siz_r(1)=siz(id_x)
            siz_r(2)=siz(id_y)

            if (id_k<1) then
              call do_interpolation2d(fields(n), ofields(n), ti, siz, siz_r)
            else
              siz_r(3)=siz(id_k)
              call do_interpolation3d(fields(n), ofields(n), ti, siz, siz_r)
            endif

          endif
        end do
      end do

      call mpp_close(ounit)
      call mpp_close(iunit)

    end subroutine process_vars


    subroutine do_interpolation2d(field,ofield,itime,siz,siz_r)

      type(fieldtype), intent(in) :: field, ofield
      integer, intent(in) :: itime, siz(:), siz_r(:)

      real, dimension(siz(1),siz(2)) :: idata
      real, dimension(siz_r(1),siz_r(2)) :: idatar
      real, dimension(ocny,ocnx) :: odata
      real :: missing
      character(len=32) :: units
      logical :: istype


      call mpp_get_atts(field, units=units, missing=missing)

      if (itime>0) then
        print *, 'pe, time =', mpp_pe(), itime
        call mpp_read(iunit, field, idata, itime)
      else
        call mpp_read(iunit, field, idata)
      endif

      istype = typedata.or.trim(adjustl(units))=='types'

      if (istype) then
        where(idata<1) idata=missing   !everything < 1 is considered as missing
        call interpolate2d_type(idata,odata,missing) 
      else
        call interpolate2d(idata,odata,missing)
      endif

      if(any(odata==missing)) then
        call rearrange2d(idata,idatar)
        call fill_miss2d(odata,idatar,missing,lmask)
        if (istype) where(odata==missing) odata=0. !zero is the missing value for type data
      end if

      call write_data2d(ofield, itime, odata)

    end subroutine do_interpolation2d


    subroutine do_interpolation3d(field,ofield,itime,siz,siz_r)
      type(fieldtype), intent(in) :: field, ofield
      integer, intent(in) :: itime, siz(:), siz_r(:)

      real, dimension(siz(1),siz(2),siz(3)) :: idata
      real, dimension(siz_r(1),siz_r(2),siz_r(3)) :: idatar
      real, dimension(siz_r(3),ocny,ocnx) :: odata
      real :: missing
      character(len=32) :: units
      logical :: istype


      call mpp_get_atts(field, units=units, missing=missing)

      if (itime>0) then
        print *, 'pe, time =', mpp_pe(), itime
        call mpp_read(iunit, field, idata, itime)
      else
        call mpp_read(iunit, field, idata)
      endif

      istype = typedata.or.trim(adjustl(units))=='types'

      if (istype) then
        where (idata<1) idata = missing   !everything < 1 is considered as missing
        call interpolate3d_type(idata,odata,missing) 
      else
        call interpolate3d(idata,odata,missing)
      end if

      if(any(odata==missing)) then
        call rearrange3d(idata,idatar)
        call fill_miss3d(odata,idatar,missing,lmask)
        if (istype) where(odata==missing) odata=0. !zero is the missing value for type data
      end if

      call write_data3d(ofield, itime, odata)

    end subroutine do_interpolation3d


    subroutine interpolate2d_type(idata, odata, missing)
      real, intent(in) :: idata(:,:), missing
      real, intent(out) :: odata(:,:)

      !local
      integer :: n, t, k, ntypes, maxtype
      integer, allocatable :: itypes(:)
      real, dimension(size(idata,1),size(idata,2)) :: iidata
      real, dimension(size(odata,1),size(odata,2)) :: iodata, rodata

      iidata = idata

      ntypes=0
      where(iidata==missing) iidata=0.

      do while (any(iidata/=0.))
        ntypes = ntypes + 1
        maxtype = maxval(iidata)
        where(iidata==maxtype) iidata=0. 
      end do

      if (ntypes==0) then
        odata=0.
        return
      endif

      allocate(itypes(ntypes))

      iidata = idata
      where(iidata==missing) iidata=0.
      ntypes=0
      do while (any(iidata/=0.))
        ntypes = ntypes + 1
        maxtype = maxval(iidata)
        itypes(ntypes) = maxtype
        where(iidata==maxtype) iidata=0. 
      end do

      iodata = missing
      rodata = 0.
      do n = 1, ntypes
        iidata = 0.
        maxtype = itypes(n)
        where(iidata==maxtype) iidata=1.
        call interpolate2d(iidata,odata,missing)
        where(odata>rodata) 
          iodata = maxtype
          rodata = odata
        endwhere
      end do

      odata = iodata
      deallocate(itypes)

    end subroutine interpolate2d_type


    subroutine interpolate3d_type(idata, odata, missing)
      real, intent(in) :: idata(:,:,:), missing
      real, intent(out) :: odata(:,:,:)

      !local
      integer :: n, t, k, ntypes, maxtype
      integer, allocatable :: itypes(:)
      real, dimension(size(idata,1),size(idata,2),size(idata,3)) :: iidata
      real, dimension(size(odata,1),size(odata,2),size(odata,3)) :: iodata, rodata

      iidata = idata

      ntypes=0
      where(iidata==missing) iidata=0.

      do while (any(iidata/=0.))
        ntypes = ntypes + 1
        maxtype = maxval(iidata)
        where(iidata==maxtype) iidata=0. 
      end do

      if (ntypes==0) then
        odata=0.
        return
      endif

      allocate(itypes(ntypes))

      iidata = idata
      where(iidata==missing) iidata=0.
      ntypes=0
      do while (any(iidata/=0.))
        ntypes = ntypes + 1
        maxtype = maxval(iidata)
        itypes(ntypes) = maxtype
        where(iidata==maxtype) iidata=0. 
      end do

      iodata = missing
      rodata = 0.
      do n = 1, ntypes
        iidata = 0.
        maxtype = itypes(n)
        where(iidata==maxtype) iidata=1.
        call interpolate3d(iidata,odata,missing)
        where(odata>rodata) 
          iodata = maxtype
          rodata = odata
        endwhere
      end do

      odata = iodata
      deallocate(itypes)

    end subroutine interpolate3d_type


    subroutine rearrange2d(idata,idatar)
      real, intent(in), dimension(:,:) :: idata
      real, intent(out), dimension(:,:) :: idatar

        do jj = 1, size(idata,2)
          do ii = 1, size(idata,1)
            idatar(ip,jp) = idata(ii,jj)
          end do
        end do

    end subroutine rearrange2d


    subroutine rearrange3d(idata,idatar)
      real, intent(in), dimension(:,:,:) :: idata
      real, intent(out), dimension(:,:,:) :: idatar

      do kk = 1, size(idata,3) 
        do jj = 1, size(idata,2)
          do ii = 1, size(idata,1)
            idatar(ip,jp,kp) = idata(ii,jj,kk)
          end do
        end do
      end do

    end subroutine rearrange3d


    subroutine interpolate2d(idata,odata,missing)
      real, intent(in), dimension(:,:) :: idata
      real, intent(out), dimension(:,:) :: odata
      real, intent(in) :: missing

      real, dimension(size(odata,1),size(odata,2)) :: denom
      integer :: n, j, i, j1, i1, swipe
      real :: rval
     
      odata = 0.
      denom = 0.
      do n = 1, xn
        ip = rxi(n)
        jp = rxj(n)
        if (debug) print *, ii, jj, kk
        if (idata(ii,jj)/=missing) then
          odata(pxj(n),pxi(n)) = odata(pxj(n),pxi(n)) + &
            idata(ii,jj) * pxf(n)
          denom(pxj(n),pxi(n)) = denom(pxj(n),pxi(n)) + pxf(n)
        endif
      end do

      where(denom>0.) 
        odata=odata/denom
      elsewhere
        odata=missing
      endwhere

    end subroutine interpolate2d


    subroutine interpolate3d(idata,odata,missing)
      real, intent(in), dimension(:,:,:) :: idata
      real, intent(out), dimension(:,:,:) :: odata
      real, intent(in) :: missing

      real, dimension(size(odata,1),size(odata,2),size(odata,3)) :: denom
      integer :: n, j, i, j1, i1, swipe, k
      real :: rval
     
      odata = 0.
      denom = 0.
      do n = 1, xn
        ip = rxi(n)
        jp = rxj(n)
        do k = 1, size(odata,1)
          kp = k
          if (debug) print *, ii, jj, kk
          if (idata(ii,jj,kk)/=missing) then
            odata(k,pxj(n),pxi(n)) = odata(k,pxj(n),pxi(n)) + &
              idata(ii,jj,kk) * pxf(n)
            denom(k,pxj(n),pxi(n)) = denom(k,pxj(n),pxi(n)) + pxf(n)
          endif
        end do
      end do

      where(denom>0.) 
        odata=odata/denom
      elsewhere
        odata=missing
      endwhere

    end subroutine interpolate3d


    subroutine fill_miss2d(odata,idatar,missing,mask)
      real, intent(inout), dimension(:,:) :: odata
      real, intent(in), dimension(:,:) :: idatar
      logical, intent(in), dimension(:,:) :: mask
      real, intent(in) :: missing

      integer :: i, j, j1, i1, swipe
      real :: rval

      do i = 1, size(odata,2)
        do j = 1, size(odata,1)
          if (odata(j,i)==missing.and.mask(j,i)) then
            j1 = rxj(max_area_loc(j,i))
            i1 = rxi(max_area_loc(j,i))
            rval = find_nn_val(i1,j1,idatar(:,:),missing,swipe) 
            if (rval==missing) then
              print *, "Could not find nearest neighbhor valid point after swipe: ", swipe
              call mpp_error(FATAL, "Could not find nearest neighbhor valid point, decrease search_frac")
            end if
            odata(j,i) = rval
          endif
        end do
      end do

    end subroutine fill_miss2d


    subroutine fill_miss3d(odata,idatar,missing,mask)
      real, intent(inout), dimension(:,:,:) :: odata
      real, intent(in), dimension(:,:,:) :: idatar
      logical, intent(in), dimension(:,:) :: mask
      real, intent(in) :: missing

      integer :: i, j, j1, i1, swipe, k
      real :: rval

      do k = 1, size(odata,1)
        do i = 1, size(odata,3)
          do j = 1, size(odata,2)
            if (odata(k,j,i)==missing.and.mask(j,i)) then
              j1 = rxj(max_area_loc(j,i))
              i1 = rxi(max_area_loc(j,i))
              rval = find_nn_val(i1,j1,idatar(:,:,k),missing,swipe) 
              if (rval==missing) then
                print *, "Could not find nearest neighbhor valid point after swipe: ", swipe
                call mpp_error(FATAL, "Could not find nearest neighbhor valid point, decrease search_frac")
              end if
              odata(k,j,i) = rval
            endif
          end do
        end do
      end do

    end subroutine fill_miss3d


    subroutine write_data2d(field, itime, dat)
      type(fieldtype), intent(in) :: field
      integer, intent(in) :: itime
      real, intent(in) :: dat(:,:)
      integer :: siz(4), n, t, k
      
      if (itime>0) then
        t = itime
        call mpp_write(ounit,field,dat,times(t))
      else
        call mpp_write(ounit,field,dat)
      endif

    end subroutine write_data2d


    subroutine write_data3d(field, itime, dat)
      type(fieldtype), intent(in) :: field
      integer, intent(in) :: itime
      real, intent(in) :: dat(:,:,:)
      integer :: siz(4), n, t, k
      
      if (itime>0) then
        t = itime
        call mpp_write(ounit,field,dat,times(t))
      else
        call mpp_write(ounit,field,dat)
      endif

    end subroutine write_data3d


    function find_nn_val(i_in,j_in,idat,missing,swipe)
      integer, intent(in) :: i_in, j_in
      real, intent(in) :: idat(:,:), missing
      integer, intent(out), optional :: swipe
      real :: find_nn_val
      integer :: is, js, ie, je
      integer :: ni, nj, nt
      integer :: e, w, s, n, t

      ni=size(idat,1)
      nj=size(idat,2)
      nt = 0

      find_nn_val=missing

      do nt = 1, max(nj,ni)/search_frac
        if(present(swipe)) swipe = nt
        is = i_in - nt; ie = i_in + nt
        js = j_in - nt; je = j_in + nt
        if ( js < 1 )  js = 1
        if ( je > nj ) je = nj
        if ( is < 1 )  is = ni + is
        if ( ie > ni ) ie = ie - ni

        do t = 0, nt
          e = i_in - t; w = i_in + t
          s = j_in - t; n = j_in + t
          if (e < 1) e = ni + e
          if (w > ni) w = w - ni
          if (n > nj) n = nj
          if (s < 1) s = 1

          ! bottom left
          if (idat(e,js)/=missing) then
            find_nn_val=idat(e,js)
            return
          endif

          !bottom right
          if (idat(w,js)/=missing) then
            find_nn_val=idat(w,js)
            return
          endif

          !top left
          if (idat(e,je)/=missing) then
            find_nn_val=idat(e,je)
            return
          endif

          !top right
          if (idat(w,je)/=missing) then
            find_nn_val=idat(w,je)
            return
          endif

          !left top
          if (idat(is,n)/=missing) then
            find_nn_val=idat(is,n)
            return
          endif

          !left bottom
          if (idat(is,s)/=missing) then
            find_nn_val=idat(is,s)
            return
          endif

          !right top
          if (idat(ie,n)/=missing) then
            find_nn_val=idat(ie,n)
            return
          endif

          !right bottom
          if (idat(ie,s)/=missing) then
            find_nn_val=idat(ie,s)
            return
          endif
        end do
      end do
    end function find_nn_val


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
          jp => ii
          kp => jj
        else
          jp => jj
          kp => ii
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
   
      ocny=siz(1)
      ocnx=siz(2)

      allocate(latbs(siz(1),siz(2)))
      allocate(latbn(siz(1),siz(2)))
      allocate(lonbe(siz(1),siz(2)))
      allocate(lonbw(siz(1),siz(2)))
      allocate(max_area_loc(siz(1),siz(2)))
      allocate(vtxoc(2,siz(1),siz(2)))
    
      allocate(vtxrg(2,size(latb,1)-1,size(lonb,1)-1))
    
      allocate(xi(2,xnsize))
      allocate(xj(2,xnsize))
      allocate(xf(3,xnsize))

      allocate(lmask(siz(1),siz(2)))

      lmask = .true.
      if (landdata) then 
        lmask=.false.
        call read_data(gridf1,'AREA_LND',latbs)
        where(latbs>0.) lmask=.true.
      endif

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
          if (ifs2n(latb)) then
            vtxrg(1,j,i)%y = latb(j)        
            vtxrg(2,j,i)%y = latb(j+1)
          else
            vtxrg(1,j,i)%y = latb(j+1)        
            vtxrg(2,j,i)%y = latb(j)
          endif
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

      max_area_loc = 0
      do n = 1, xn
        if (max_area_loc(pxj(n),pxi(n)) == 0) max_area_loc(pxj(n),pxi(n)) = n
        if ( pxf(n) > pxf(max_area_loc(pxj(n),pxi(n))) ) max_area_loc(pxj(n),pxi(n)) = n
      end do

      if (any(max_area_loc==0)) call mpp_error(FATAL,"max_area_loc == 0")

      !print *, maxval(pxi), maxval(pxj), maxval(rxi), maxval(rxj)
    
      deallocate(xi)
      deallocate(xj)
      deallocate(xf)
    
    end subroutine make_xgrid_oc2rg


    subroutine modify_fields(flds)
      type(fieldtype), intent(in) :: flds(:)
      integer :: i, j, ndimv, ndimv1, idk
      type(axistype) :: axes(ndim), axes1(ndim)
      character(len=128) :: nm, units, longname
      real :: minv, maxv, scalev, add, missing

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

        idk=id_Kaxis(axes(1:ndimv))

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


    function id_Kaxis(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: id_Kaxis, i

      id_Kaxis=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lat_nm) .or. &
            trim(get_name(axes(i)))==trim(lon_nm) .or. &
            trim(get_name(axes(i)))==trim(time_axis_name)) cycle
        id_Kaxis=i
        exit
      end do
    end function id_Kaxis


    function id_Yaxis(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: id_Yaxis, i

      id_Yaxis=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lat_nm)) then
          id_Yaxis=i
          exit
        endif
      end do
    end function id_Yaxis


    function id_Xaxis(axes)
      type(axistype), intent(inout) :: axes(:)
      integer :: id_Xaxis, i

      id_Xaxis=0
      do i = 1, size(axes,1)
        if (trim(get_name(axes(i)))==trim(lon_nm)) then
          id_Xaxis=i
          exit
        endif
      end do
    end function id_Xaxis


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



    subroutine broadcast_char(cdat)
      character(len=*), intent(inout) :: cdat
      character(len=len(cdat)) :: cdata(1)

      cdata=cdat

      call mpp_broadcast(cdata,len(cdat),mpp_root_pe())

      cdat=cdata(1)

    end subroutine broadcast_char

    subroutine handle_err(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
        call mpp_error(FATAL,trim(nf90_strerror(status)))
      end if
    end subroutine handle_err

end program amfi_xgrid
