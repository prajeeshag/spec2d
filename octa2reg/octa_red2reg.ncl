
nlon = 204
nlat = 94

if (mod(nlat,2).ne.0) then
	print("FATAL: NLAT should be a multiple of 2")
	exit
end if	

yunits=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", \
	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
xunits=(/"degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"/)

lonsperlat=new((/nlat/),integer)
nplon = nlon - (nlat/2-1)*4
lonsperlat(0) = nplon
lonsperlat(nlat-1) = nplon
do i = 1, nlat/2-1
	lonsperlat(i) = lonsperlat(i-1) + 4
	lonsperlat(nlat-1-i) = lonsperlat(i)
end do

nlat@double = True
nlon@double = True
lonf = lonGlobeF(nlon, "lon", "longitude", "degrees_E")
latf = latGau(nlat, "lat", "latitude", "degrees_N")

ocnx = nplon+max(lonsperlat)
ocny = nlat/2

is=new((/2,ocny/),integer)
ie=new((/2,ocny/),integer)
ilen=new((/2,ocny/),integer)
pack=new((/2,ocny/),integer)

is(0,:) = 0
ie(1,:) = ocnx-1

pack = -1
do i = 0, nlat/4-1
	pack(0,i*2) = i
	pack(0,i*2+1) = nlat-1-i
	pack(1,i*2) = nlat/2-i-1
	pack(1,i*2+1) = nlat/2+i
end do

if (mod(nlat,4).ne.0) then
	pack(0,nlat/2-1) = nlat/4
	pack(1,nlat/2-1) = nlat-nlat/4-1
end if

do i = 0, ocny-1
	ie(0,i) = is(0,i) + lonsperlat(pack(0,i))
	is(1,i) = ie(0,i) + 1
	ilen(0,i) = lonsperlat(pack(0,i))
	ilen(1,i) = lonsperlat(pack(1,i))
end do

do i = 0, nlat/2-1
	print("pack = " + pack(0,i) + " " +pack(1,i) +" "+is(1,i)+" "+ilen(1,i))
end do

exit

;--------------------------------------------------------------------------------	
function get_cart(ax:string,fid:file)
begin
;--------------------------------------------------------------------------------	
	var=fid->$ax$
	cart = ""
	atts = getvaratts(var)
	if (all(ismissing(atts))) then
		print("Cannot determine cartesian axis attribute for "+ax)
		exit
	end if	
	if (any(atts.eq."axis")) then
		cart = var@axis
		return(cart)
	else if (any(atts.eq."cartesian_axis")) then
		cart = var@cartesian_axis
		return(cart)
	else
		print("Cannot determine cartesian axis attribute for "+ax)
		exit
	end if
	end if
end

;--------------------------------------------------------------------------------	
function find_axis_nms(varnm:string,fid:file)
local j, xnm, ynm, vdmnm
begin
;--------------------------------------------------------------------------------	
	vdmnm=getfilevardims(fid,varnm)
	
	xnm=""
	ynm=""
	znm=""
	tnm=""
	do j = 0, dimsizes(vdmnm)-1
		cart2=get_cart(vdmnm(j),fid)
		if (cart2.eq."X".or.cart2.eq."x") then
			xnm=vdmnm(j)
		else if (cart2.eq."Y".or.cart2.eq."y") then
			ynm=vdmnm(j)	
		else if (cart2.eq."Z".or.cart2.eq."z") then
			znm=vdmnm(j)	
		else if (cart2.eq."t".or.cart2.eq."T") then
			tnm=vdmnm(j)	
		end if
		end if
		end if
		end if
	end do

	if (ynm.eq."".or.xnm.eq."") then
		print("Cannot find X or Y coordinate for variable: "+varnm)
		exit
	end if

	axnm = (/ynm,xnm/)
	if (znm.ne."") then
		axnm1 = axnm
		delete(axnm)
		axnm=new(dimsizes(axnm1)+1,typeof(axnm1))
		axnm(0) = znm
		axnm(1:) = axnm1
		delete(axnm1)
	end if
	
	if (tnm.ne."") then
		axnm1 = axnm
		delete(axnm)
		axnm=new(dimsizes(axnm1)+1,typeof(axnm1))
		axnm(0) = tnm
		axnm(1:) = axnm1
		delete(axnm1)
	end if

	return (axnm)
end

;-------------------------------------------------------------------------------
function stack_and_fold(dati:numeric)
begin
;--------------------------------------------------------------------------------	

	printVarSummary(dati)

return(True)

end


;--------------------------------------------------------------------------------	
procedure toOcta_and_write(fo:file,fi:file, vnm:string, xynm[*]:string, \
						   mxlon:integer)
local axnm, i, lonin, latin, nlatin, nlonin, dati
begin
;--------------------------------------------------------------------------------	

	latin=fi->$xynm(1)$
	nlatin=dimsizes(latin)

	lonin=fi->$xynm(0)$
	nlonin=dimsizes(lonin)

	if (dimsizes(xynm).eq.3) then
		dati = fi->$vnm$($xynm(0)$|:,$xynm(1)$|:,$xynm(2)$|:)
		xnm = xynm(2)
		ynm = xynm(1)
	else if (dimsizes(xynm).eq.4) then
		dati = fi->$vnm$($xynm(0)$|:,$xynm(1)$|:,$xynm(2)$|:,$xynm(3)$|:)
		xnm = xynm(3)
		ynm = xynm(2)
	else if (dimsizes(xynm).eq.2) then
		dati = fi->$vnm$($xynm(0)$|:,$xynm(1)$|:)
		xnm = xynm(1)
		ynm = xynm(0)
	end if
	end if
	end if

    datf = area_conserve_remap_Wrap(dati&$xnm$, dati&$ynm$, dati, lonf, latf, False)
	fo->$vnm$ = datf
	
	dato=stack_and_fold(datf)

	delete(dati)
	delete(datf)
end


varlist = (/"tp1"/)
fi = addfile("amfi_res.nc","r")

system("rm -rf out.nc")
fo = addfile("out.nc","c")

fvnms = getfilevarnames(fi)

do i = 0, dimsizes(fvnms)-1
	if (any(fvnms(i).eq.getfilevardims(fi,fvnms(i)))) then
		print("skipping dim var: " + fvnms(i))
		continue
	end if
	if (dimsizes(varlist).ne.0) then
		if (.not.any(fvnms(i).eq.varlist)) then
			print("skipping var: " + fvnms(i))
			continue
		end if
	end if
	axnm=find_axis_nms(fvnms(i),fi)
	toOcta_and_write(fo,fi,fvnms(i),axnm,nlon)	
	delete(axnm)
end do

