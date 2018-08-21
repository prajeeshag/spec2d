#!/bin/bash 

outdir="Unfold"
usage() { echo "Usage: $0 -x nlon -y nlat [-n] [-r] [-v varlist] \
[-d outdir] [-p optionlist] inputfiles" 1>&2; exit 1;}

npack=2
reduce=1

while getopts 'x:d:y:p:v:nr' flag; do
    case "${flag}" in
	p) oplist="$OPTARG" ;;
	x) NLON="$OPTARG" ;;
	y) NLAT="$OPTARG" ;;
	v) valist="$OPTARG" ;;
	d) outdir="$OPTARG" ;;
	n) npack=1 ;;
	r) reduce=0 ;;
    *)
		usage	
		;;
    esac
done

shift $(expr $OPTIND - 1)

infiles=$*

if [[ -z $infiles ]] || [[ -z $NLON ]] || [[ -z $NLAT ]]; then
	usage
fi

varlist="\" \""
if ! [[ -z $valist ]]; then
varlist="\"$valist\""
fi

filelist="\"$infiles\""

optlist="\" \""
if ! [[ -z $oplist ]]; then
optlist="\"$oplist\""
fi

tfile=$(mktemp)

mkdir -p $outdir 

echo $tfile

cat <<EOF > $tfile

begin

;********************************************************************************	
;GLOBAL SECTION
NLON = $NLON
NLAT = $NLAT

ifpack="$npack".eq.2
reduce="$reduce".eq.1

NPACK=1
if (ifpack) then
    NPACK=2
    ;reduce = True
end if

filelist = rm_single_dims(str_split_csv($filelist," ",0))

optlist = rm_single_dims(str_split_csv($optlist,":",0))

varlist = rm_single_dims(str_split_csv($varlist,",",0))

if (mod(NLAT,2).ne.0) then
	print("FATAL: NLAT should be a multiple of 2")
	exit
end if	

YUNITS=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", \
	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
XUNITS=(/"degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"/)

LONSPERLAT=new((/NLAT/),integer)

NPLON = NLON - (NLAT/2-1)*4

print("NPLON="+NPLON)

LONSPERLAT = NLON

if (reduce) then
	LONSPERLAT(0) = NPLON
	LONSPERLAT(NLAT-1) = NPLON
	do i = 1, NLAT/2-1
		LONSPERLAT(i) = LONSPERLAT(i-1) + 4
		LONSPERLAT(NLAT-1-i) = LONSPERLAT(i)
	end do
end if

NLAT@double = True
NLON@double = True
LONF = lonGlobeF(NLON, "lon", "longitude", "degrees_E")
LATF = latGau(NLAT, "lat", "latitude", "degrees_N")

OCNX = (NPACK-1)*NPLON+max(LONSPERLAT)
OCNY = NLAT/NPACK

IS=new((/NPACK,OCNY/),integer)
IE=new((/NPACK,OCNY/),integer)
ILEN=new((/NPACK,OCNY/),integer)
PACK=new((/NPACK,OCNY/),integer)

if (ifpack) then
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

else

	    IS(0,:) = 0
    PACK = -1
    do i = 0, NLAT/4-1
        PACK(0,4*i) = i
        PACK(0,4*i+1) = NLAT/2-i-1
        PACK(0,4*i+2) = NLAT-1-i
        PACK(0,4*i+3) = NLAT/2+i
    end do
    
    if (mod(NLAT,4).ne.0) then
        PACK(0,NLAT-2) = NLAT/4
        PACK(0,NLAT-1) = NLAT-NLAT/4-1
    end if

    do i = 0, OCNY-1
        IE(0,i) = IS(0,i) + LONSPERLAT(PACK(0,i)) - 1
        ILEN(0,i) = LONSPERLAT(PACK(0,i))
    end do

end if	

OCNX  = (NPACK-1)*min(LONSPERLAT)+max(LONSPERLAT)

print("OCNX="+OCNX+" OCNY="+OCNY)

OCLON = new((/OCNY,OCNX/),typeof(LONF))
OCLAT = new((/OCNY,OCNX/),typeof(LATF))
OCLON = 0.
OCLAT = 90.

do j = 0, NPACK-1
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

NOPT = dimsizes(optlist)

OPTDIC = new((/2,NOPT/),string)

do i = 0, NOPT-1
	OPTDIC(:,i) = rm_single_dims(str_split_csv(optlist(i),"=",0))
end do

delete(optlist)

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
		status_exit(1)
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

	if (any(atts.eq."units")) then
		if (any(XUNITS.eq.var@units)) then
			return("X")
		else if (any(YUNITS.eq.var@units)) then
			return("Y")
		else if (any(rm_single_dims(str_split_csv(var@units," ",0)).eq."since")) then
			return("T")
		end if
		end if
		end if
	end if

	print("FATAL: Cannot determine cartesian axis attribute for "+ax)
	status_exit(1)
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
	
	noxy=False
	if (ynm.eq."".or.xnm.eq."") then
		noxy=True
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

	axnm@noxy = noxy

	return (axnm)
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
function unstack_and_unfold_conserve(dati:numeric,xi:numeric,yi:numeric)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd
begin
;--------------------------------------------------------------------------------	

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	howmany = 1
	
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	osiz=new(ndim,integer)

	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do

	osiz(ndim-2) = NLAT
	osiz(ndim-1) = NLON

	datii = reshape(dati,(/howmany,OCNY,OCNX/))
	datoi = new((/howmany,NLAT,NLON/),typeof(dati)) 
	
	ncopy = 4
	latfi = latGau(ncopy, "lat", "latitude", "degrees_N")

    do j = 0, NPACK-1
        do i = 0, OCNY-1
            lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
			ip = PACK(j,i)
			rfi = new((/howmany,ncopy,ILEN(j,i)/),typeof(datii))
			do k = 0, ncopy-1
				rfi(:,k,:) = datii(:,i,IS(j,i):IE(j,i))
			end do
            rfo = area_conserve_remap(lonc,latfi,rfi,LONF,latfi,False)
			do k = 0, 1
            	datoi(:,ip,:) = (/rfo(:,(ncopy+1)/2-1,:)/)
			end do
			delete(lonc)
			delete(rfi)
        end do
    end do

	dato = reshape(datoi,osiz)

	copy_VarCoords_2(dati,dato)
	copy_VarAtts(dati,dato)
	
	nnd = ndim - 1
	dato!nnd = "lon"
	nnd = ndim - 2
	dato!nnd = "lat"
	dato&lon = LONF
	dato&lat = LATF
	dato&lon@axis = "X"
	dato&lat@axis = "Y"

	return(dato)
end
;end unstack_and_unfold_conserve


;-------------------------------------------------------------------------------
function unstack_and_unfold_linint2(dati:numeric,xi:numeric,yi:numeric)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd
begin
;--------------------------------------------------------------------------------	

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	howmany = 1
	
	do i = 0, ndim-3
		howmany=siz(i)*howmany
	end do
	
	osiz=new(ndim,integer)

	do i = 0, ndim-3
		osiz(i) = siz(i)
	end do

	osiz(ndim-2) = NLAT
	osiz(ndim-1) = NLON
	
	if (isinteger(dati)) then
		datii = tofloat(reshape(dati,(/howmany,OCNY,OCNX/)))
	else
		datii = reshape(dati,(/howmany,OCNY,OCNX/))
	end if

	datoi = new((/howmany,NLAT,NLON/),typeof(datii)) 

   	do j = 0, NPACK-1
   	    do i = 0, OCNY-1
   	        lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
   	        datoi(:,PACK(j,i),:) = linint1(lonc,datii(:,i,IS(j,i):IE(j,i)),True,LONF,0)
   	        ;datoi(:,PACK(j,i),0:ILEN(j,i)-1) = datii(:,i,IS(j,i):IE(j,i))
   			delete(lonc)
   	    end do
   	end do

	dato = reshape(datoi,osiz)

	copy_VarCoords_2(dati,dato)
	copy_VarAtts(dati,dato)
	
	nnd = ndim - 1
	dato!nnd = "lon"
	nnd = ndim - 2
	dato!nnd = "lat"
	dato&lon = LONF
	dato&lat = LATF
	dato&lon@axis = "X"
	dato&lat@axis = "Y"

	return(dato)
end
;end unstack_and_unfold_linint2 

;--------------------------------------------------------------------------------	
procedure fromOcta_and_write(fo:file, fi:file, vnm:string, intmethd:string, xynm:string)
local i, lonin, latin, NLATin, NLONin, dati, fo, fi, vnm, xynm, intmethd, ndim, siz, conv2int
begin
;--------------------------------------------------------------------------------	

	ndim = dimsizes(xynm) 

	if (ndim.gt.3) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:,\$xynm(3)\$|:)
	else if (ndim.eq.3) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:)
	else 
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:)
	end if
	end if
			
	siz = dimsizes(dati) 

	xi = dati&\$xynm(ndim-1)\$
	yi = dati&\$xynm(ndim-2)\$

	if (intmethd.eq."linint") then
		dato=unstack_and_unfold_linint2(dati,xi,yi)
	else if (intmethd.eq."conserve") then
		dato=unstack_and_unfold_conserve(dati,xi,yi)
	else
		print("FATAL: invalid interpmethod: "+intmethd)
		print("NOTE: available interpmethods are: 'linint', 'conserve'.")
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

do f = 0, dimsizes(filelist)-1
    filnm=filelist(f)
	print(" "+filnm)
	fi = addfile(filnm,"r")
	outfnm = "$outdir/"+filnm
	if (fileexists(outfnm)) then
		print("Error: File ("+outfnm+") already exist")
		status_exit(1)
	end if
	fo = addfile(outfnm,"c")

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
		axnm = find_axis_nms(fvnms(i),fi)
		if (axnm@noxy) then
			print("No x or y found for : "+fvnms(i))
			;fo->\$fvnms(i)\$ = fi->\$fvnms(i)\$
			continue
		end if
		intpmthd=get_option_val("interpmethod","linint",n)
		print("unstacking and unfolding "+fvnms(i)+" with options interpmethod = "+intpmthd)
		fromOcta_and_write(fo,fi,fvnms(i),intpmthd, axnm)	
		delete(axnm)
	end do
	delete(fvnms)
end do
end

EOF


#cat $tfile
ncl $tfile


