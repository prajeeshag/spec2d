#!/bin/bash 

usage() { echo "Usage: $0 -x NLON -y NLAT -i inputfile.nc -o outfile.nc [-v vars] [-p optionlist]" 1>&2; exit 1;}

while getopts 'x:y:i:o:v:p:' flag; do
    case "${flag}" in
    x) NLON="$OPTARG" ;;
    y) NLAT="$OPTARG" ;;
    i) infile="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
	v) valist="$OPTARG" ;;
	p) oplist="$OPTARG" ;;
    *)
		usage	
		;;
    esac
done


if [[ -z $NLON ]] || [[ -z $outfile ]] || \
       [[ -z $NLAT ]] || [[ -z $infile ]]; then
	usage
fi

varlist="\" \""
if ! [[ -z $valist ]]; then
varlist="\"$valist\""
fi

optlist="\" \""
if ! [[ -z $oplist ]]; then
optlist="\"$oplist\""
fi

tfile=$(mktemp)

echo $tfile

cat <<EOF > $tfile

;external DISRGRD "./discre.so"

begin

;********************************************************************************	
;GLOBAL SECTION
NLON = $NLON
NLAT = $NLAT

if (mod(NLAT,2).ne.0) then
	print("FATAL: NLAT should be a multiple of 2")
	exit
end if	

YUNITS=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", \
	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
XUNITS=(/"degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"/)

LONSPERLAT=new((/NLAT/),integer)
NPLON = NLON - (NLAT/2-1)*4
LONSPERLAT(0) = NPLON
LONSPERLAT(NLAT-1) = NPLON
do i = 1, NLAT/2-1
	LONSPERLAT(i) = LONSPERLAT(i-1) + 4
	LONSPERLAT(NLAT-1-i) = LONSPERLAT(i)
end do

NLAT@double = True
NLON@double = True
LONF = lonGlobeF(NLON, "lon", "longitude", "degrees_E")
LATF = latGau(NLAT, "lat", "latitude", "degrees_N")

OCNX = NPLON+max(LONSPERLAT)
OCNY = NLAT/2

IS=new((/2,OCNY/),integer)
IE=new((/2,OCNY/),integer)
ILEN=new((/2,OCNY/),integer)
PACK=new((/2,OCNY/),integer)

IS(0,:) = 0
IE(1,:) = OCNX-1

PACK = -1
do i = 0, NLAT/4-1
	PACK(0,i*2) = i
	PACK(0,i*2+1) = NLAT-1-i
	PACK(1,i*2) = NLAT/2-i-1
	PACK(1,i*2+1) = NLAT/2+i
end do

if (mod(NLAT,4).ne.0) then
	PACK(0,NLAT/2-1) = NLAT/4
	PACK(1,NLAT/2-1) = NLAT-NLAT/4-1
end if

do i = 0, OCNY-1
	IE(0,i) = IS(0,i) + LONSPERLAT(PACK(0,i)) - 1
	IS(1,i) = IE(0,i) + 1
	ILEN(0,i) = LONSPERLAT(PACK(0,i))
	ILEN(1,i) = LONSPERLAT(PACK(1,i))
end do

OCLON = new((/OCNY,OCNX/),typeof(LONF))
OCLAT = new((/OCNY,OCNX/),typeof(LATF))

do j = 0, 1
    do i = 0, OCNY-1
        lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
        OCLON(i,IS(j,i):IE(j,i)) = lonc
        OCLAT(i,IS(j,i):IE(j,i)) = LATF(PACK(j,i))
        delete(lonc)
    end do
end do

OCNXNY = OCNY*OCNX

OCLON1d = reshape(OCLON,(/OCNXNY/))
OCLAT1d = reshape(OCLAT,(/OCNXNY/))

varlist = rm_single_dims(str_split_csv($varlist,",",0))

optlist = rm_single_dims(str_split_csv($optlist,":",0))

NOPT = dimsizes(optlist)

OPTDIC = new((/2,NOPT/),string)

do i = 0, NOPT-1
	OPTDIC(:,i) = rm_single_dims(str_split_csv(optlist(i),"=",0))
end do

;GLOBAL SECTION END
;********************************************************************************	



;FUNCTION SECTION
;********************************************************************************	

;--------------------------------------------------------------------------------	
function get_option_val(nm:string,default:string,nv:integer)
local nm, val, i, nv, val1
begin
;--------------------------------------------------------------------------------

val1=""
do i = 0, NOPT-1
	if (nm .eq. OPTDIC(0,i)) then
		val1 = OPTDIC(1,i)
		break
	end if	
end do

if (val1.eq."") then
	return(default)
end if

vals = rm_single_dims(str_split_csv(val1,",",0))

if (nv.lt.0) then
	return(vals(0))
end if

if (dimsizes(vals).eq.1) then
	val = vals(0)
	return(val)
end if

if (nv .gt. dimsizes(vals)-1) then
	print("FATAL: number of option values for option "+nm+ \
		  " is less than the number of variable given in the list")
	exit
end if

return(vals(nv))

end

;--------------------------------------------------------------------------------	
function get_cart(ax:string,fid:file)
local i, ax, fid
begin
;--------------------------------------------------------------------------------	
	var=fid->\$ax\$
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
	end if
	end if

	if (any(atts.eq."calendar")) then
		cart = "T"
		return(cart)
	end if

	print("Cannot determine cartesian axis attribute for "+ax)
	exit
end

;--------------------------------------------------------------------------------	
function find_axis_nms(varnm:string,fid:file)
local j, xnm, ynm, vdmnm, i, varnm, fid, axnm
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
function stack_and_fold_unstrct(datf:numeric,xi:numeric,yi:numeric,ongrid:string)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, ongrid, ongrd
begin
;--------------------------------------------------------------------------------	

	siz = dimsizes(datf)
	ndim = dimsizes(siz)

	nx=siz(ndim-1)
	ny=siz(ndim-2)
	howmany = 1
	
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	osiz=new(ndim,integer)

	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do
	osiz(ndim-2) = OCNY
	osiz(ndim-1) = OCNX
	
	xxi = xi
	yyi = yi

	ongrd=stringtochar(ongrid)
	if (ongrd(0).eq."T") then
		delete(yyi)
		yyi = LATF
	end if

	if (ongrd(1).eq."T") then
		delete(xxi)
		xxi = LONF
	end if

	dati = linint2_points_Wrap(xxi, yyi, datf, True, OCLON1d, OCLAT1d, 0) 

	if (any(ismissing(dati))) then
		print("WARNING: packed field contains missing values!!!!")
	end if

	dato = reshape(dati,osiz)

	copy_VarCoords_2(dati,dato)
	copy_VarAtts(dati,dato)

	nnd = ndim - 1
	dato!nnd = "x"
	nnd = ndim - 2
	dato!nnd = "y"
	dato&x = ispan(1,OCNX,1)
	dato&y = ispan(1,OCNY,1)
	dato&x@axis = "X"
	dato&y@axis = "Y"

	return(dato)

end

;--------------------------------------------------------------------------------	
function regrid_discrete_Wrap(xi:numeric, yi:numeric, fi:numeric, xo:numeric, yo:numeric, opt:logical)
local i, j, k, n, xi, yi, fi, xo, yo, fo, mo, ntype, \
      nx, ny, howmany, ndim, siz, ltypes, maxtype, osiz, opt
begin
; Assumes that type 0 is a missing type, and all types are positive
;--------------------------------------------------------------------------------	

	if (.not.isinteger(fi)) then
		print("FATAL: fi of regrid_discrete should be of integer type")
		exit 
	end if
	
	siz = dimsizes(fi)
	ndim = dimsizes(siz)

	if (ndim.lt.2) then
		print("FATAL: regrid_discrete: fi should be atleast a 2D field")
		exit
	end if
	
	nx=siz(ndim-1)
	ny=siz(ndim-2)
	howmany = 1
	
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	nxo=dimsizes(xo)
	nyo=dimsizes(yo)
	osiz=new(ndim,integer)
	
	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do
	osiz(ndim-2) = nyo
	osiz(ndim-1) = nxo
	
	if (any(fi<0)) then
		print("FATAL: fi for regrid_discrete should be positive values")
		exit
	end if
	
	maxtype = max(fi)
	
	ntype = 0
	do i = 0, maxtype
		if(any(fi.eq.i)) then
			ntype = ntype + 1
		end if
	end do
	
	if (ntype.eq.1) then
		fo = new(osiz,integer)
		fo = maxtype
		return(fo)
	end if
	
	ltypes = new(ntype,integer)
	
	n = 0
	do i = 0, maxtype
	    if(any(fi.eq.i)) then
	      	ltypes(n) = i 
		 	n = n + 1
	    end if
	end do
	
	fii = reshape(fi,(/howmany,ny,nx/))
	fii!2 = "lon"
	fii!1 = "lat"
	fii&lon = xi
	fii&lat = yi

	buffi = new((/howmany,ny,nx/),float)
	buffo = new((/howmany,nyo,nxo,ntype/),float)
	buffo1 = new((/howmany,nyo,nxo/),integer)
	
	do i = 0, ntype-1
		buffi = where(fii.eq.ltypes(i),1.,0.)
		copy_VarCoords(fii,buffi)
		buffo(:,:,:,i) = area_conserve_remap(xi, yi, buffi, xo, yo, opt)
	end do
	
	do i = 0, nxo-1 
		do j = 0, nyo-1 
			do k = 0, howmany-1
				buffo1(k,j,i) = ltypes(maxind(buffo(k,j,i,:)))
			end do
		end do
	end do

	fo = reshape(buffo1,osiz)
	
	nnd = ndim-1
	fo!nnd = "lon"
	nnd = ndim-2
	fo!nnd = "lat"
	fo&lon = xo
	fo&lat = yo
	copy_VarCoords_2(fi,fo)
	
	return(fo)
end

;-------------------------------------------------------------------------------
function stack_and_fold_dtype(fi:numeric,xi:numeric,yi:numeric,ongrid:string,conv2int:logical)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd, conv2int, fi, dati
begin
;--------------------------------------------------------------------------------

	if (.not.isinteger(fi)) then
		if (.not.conv2int) then
			print("FATAL: input array to stack_and_fold_dtype should be of integer type")
			exit
		end if
		datf = tointeger(fi)
		copy_VarMeta(fi,datf)
	else
		datf = fi
	end if
	
	ongrd = stringtochar(ongrid)

	xxi = xi
	yyi = yi

	if (ongrd(0).eq."T") then
		delete(yyi)
		yyi = LATF
	end if

	if (ongrd(1).eq."T") then
		delete(xxi)
		xxi = LONF
	end if

	if (ongrd(0).eq."T".and.ongrd(1).eq."T") then
		dati = datf
	else
    	dati = regrid_discrete_Wrap(xxi, yyi, datf, LONF, LATF, False)
	end if

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	nx = siz(ndim-1)	
	ny = siz(ndim-2)	
	howmany = 1 
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	osiz=new(ndim,integer)

	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do

	osiz(ndim-2) = OCNY
	osiz(ndim-1) = OCNX

	datii = reshape(dati,(/howmany,ny,nx/))
	datoi = new((/howmany,OCNY,OCNX/),typeof(dati)) 

	ncopy = 4
	rfi = new((/howmany,ncopy,nx/),typeof(datii))
	latfi = latGau(ncopy, "lat", "latitude", "degrees_N")

    do j = 0, 1
        do i = 0, OCNY-1
            lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
			ip = PACK(j,i)
			do k = 0, ncopy-1
				rfi(:,k,:) = datii(:,ip,:)
			end do
            rfo = regrid_discrete_Wrap(LONF,latfi,rfi,lonc,latfi,False)
			do k = 0, 1
            	datoi(:,i,IS(j,i):IE(j,i)) = (/rfo(:,(ncopy+1)/2-1,:)/)
			end do
			delete(lonc)
			delete(rfo)
        end do
    end do

	dato = reshape(datoi,osiz)

	copy_VarCoords_2(dati,dato)

	copy_VarCoords_2(dati,dato)
	copy_VarAtts(dati,dato)

	nnd = ndim - 1
	dato!nnd = "x"
	nnd = ndim - 2
	dato!nnd = "y"
	dato&x = ispan(1,OCNX,1)
	dato&y = ispan(1,OCNY,1)
	dato&x@axis = "X"
	dato&y@axis ="Y"

	return(dato)
end

;-------------------------------------------------------------------------------
function stack_and_fold_linint2(datf:numeric,xi:numeric,yi:numeric,ongrid:string)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd
begin
;--------------------------------------------------------------------------------	

	ongrd = stringtochar(ongrid)

	xxi = xi
	yyi = yi

	if (ongrd(0).eq."T") then
		delete(yyi)
		yyi = LATF
	end if

	if (ongrd(1).eq."T") then
		delete(xxi)
		xxi = LONF
	end if

	if (ongrd(0).eq."T".and.ongrd(1).eq."T") then
		dati = datf
	else
    	dati = linint2_Wrap(xxi, yyi, datf, True, LONF, LATF, 0)
	end if

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	nx=siz(ndim-1)
	ny=siz(ndim-2)
	howmany = 1
	
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	osiz=new(ndim,integer)

	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do

	osiz(ndim-2) = OCNY
	osiz(ndim-1) = OCNX

	datii = reshape(dati,(/howmany,ny,nx/))
	datoi = new((/howmany,OCNY,OCNX/),typeof(dati)) 
	
    do j = 0, 1
        do i = 0, OCNY-1
            lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
            datoi(:,i,IS(j,i):IE(j,i)) = linint1(LONF,datii(:,PACK(j,i),:),True,lonc,0)
			delete(lonc)
        end do
    end do

	dato = reshape(datoi,osiz)

	copy_VarCoords_2(dati,dato)
	copy_VarAtts(dati,dato)
	
	nnd = ndim - 1
	dato!nnd = "x"
	nnd = ndim - 2
	dato!nnd = "y"
	dato&x = ispan(1,OCNX,1)
	dato&y = ispan(1,OCNY,1)
	dato&x@axis = "X"
	dato&y@axis = "Y"

	return(dato)

end


;--------------------------------------------------------------------------------	
procedure toOcta_and_write(fo:file, fi:file, vnm:string, xynm[*]:string, intmethd:string, ongrid:string, conv2int:logical)
local i, lonin, latin, NLATin, NLONin, dati, fo, fi, vnm, xynm, intmethd, ndim, siz, conv2int
begin
;--------------------------------------------------------------------------------	

	ndim = dimsizes(xynm) 

	if (ndim.eq.3) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:)
	else if (ndim.eq.4) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:,\$xynm(3)\$|:)
	else if (ndim.eq.2) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:)
	end if
	end if
	end if

	siz = dimsizes(dati)
	xi = fi->\$xynm(ndim-1)\$
	yi = fi->\$xynm(ndim-2)\$

	if (intmethd.eq."unstruct") then
		dato=stack_and_fold_unstrct(dati,xi,yi,ongrid)
	else if (intmethd.eq."dtype") then
		dato=stack_and_fold_dtype(dati,xi,yi,ongrid,conv2int)
	else if (intmethd.eq."linint") then
		dato=stack_and_fold_linint2(dati,xi,yi,ongrid)
	else
		print("FATAL: invalid interpmethod: "+intmethd)
		print("NOTE: available interpmethods are: unstruct, dtype, linint")
	end if
	end if
	end if

	fo->\$vnm\$ = dato

	delete(dati)
	delete(dato)
	return
end


;********************************************************************************	
;END FUNCTION SECTION



;********************************************************************************	
;--> Main SECTION


fi = addfile("$infile","r")

system("rm -rf $outfile")
fo = addfile("$outfile","c")

fvnms = getfilevarnames(fi)

do i = 0, dimsizes(fvnms)-1
	if (any(fvnms(i).eq.getfilevardims(fi,fvnms(i)))) then
		continue
	end if
	n = -1
	if (any(varlist.ne." ")) then
		do j = 0, dimsizes(varlist)-1
			if (fvnms(i).eq.varlist(j)) then
				n = j
			end if
		end do
		if (n.lt.0) then
			continue
		end if
	end if
	axnm=find_axis_nms(fvnms(i),fi)
	intpmthd=get_option_val("interpmethod","linint",n)
	ongrid=get_option_val("ongrid","FF",n)
	conv2int=get_option_val("conv2int","F",n).eq."T"
	print("stacking and folding "+fvnms(i)+" with options interpmethod = "+intpmthd+", ongrid = "+ongrid)
	toOcta_and_write(fo,fi,fvnms(i),axnm,intpmthd,ongrid,conv2int)	
	delete(axnm)
end do

end

EOF


#cat $tfile
ncl $tfile


