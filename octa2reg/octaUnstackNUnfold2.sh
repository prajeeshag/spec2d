#!/bin/bash 

outdir="Unfold"
usage() { echo "Usage: $0 [-x path_to_xgrd] [-v varlist] \
[-d outdir] [-p optionlist] inputfiles" 1>&2; exit 1;}

xgrd="p_xgrd.nc"

while getopts 'x:d:p:v:' flag; do
    case "${flag}" in
	x) xgrd="$OPTARG" ;;
	p) oplist="$OPTARG" ;;
	v) valist="$OPTARG" ;;
	d) outdir="$OPTARG" ;;
    *)
		usage	
		;;
    esac
done

shift $(expr $OPTIND - 1)

infiles=$*

if [[ -z $infiles ]] ; then
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
xgrd="$xgrd"

if (fileexists(xgrd)) then
	fxg = addfile(xgrd,"r")

	pxi = fxg->pxi - 1
	pxj = fxg->pxj - 1

	rxi = fxg->rxi - 1
	rxj = fxg->rxj - 1
	xf = doubletofloat(fxg->xf)

	NLON = fxg@maxlon
	NLAT = fxg@nlat
	OCNX = fxg@ocnx
	OCNY = fxg@ocny
	siz = dimsizes(pxi)
	nxgrd = siz(0)
	delete(siz)
else
	print("Error: File ("+xgrd+") does not exist")
	status_exit(1)
end if

filelist = rm_single_dims(str_split_csv($filelist," ",0))

optlist = rm_single_dims(str_split_csv($optlist,":",0))

varlist = rm_single_dims(str_split_csv($varlist,",",0))

YUNITS=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", \
	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
XUNITS=(/"degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"/)

LONF = lonGlobeF(NLON, "lon", "longitude", "degrees_E")
LATF = latGau(NLAT, "lat", "latitude", "degrees_N")

print("OCNX="+OCNX+" OCNY="+OCNY)

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

	missing = dati@missing_value		

	siz = dimsizes(dati) 

	if (siz(ndim-1) .ne. OCNX) then
        print("FATAL: x dimension size of data /= OCNX")
        status_exit(1)
    end if

    if (siz(ndim-2) .ne. OCNY) then
        print("FATAL: y dimension size of data /= OCNY")
        status_exit(1)
    end if

	osiz = siz
	osiz(ndim-1) = NLON
	osiz(ndim-2) = NLAT

	siz1 = (/1,OCNY,OCNX/)
	do n = 0, ndim-3
		siz1(0) = siz1(0) * siz(n)
	end do

	osiz1 = (/1,NLAT,NLON/)
	do n = 0, ndim-3
		osiz1(0) = osiz1(0) * osiz(n)
	end do

	rdato = new(osiz1,typeof(dati))
	rdato = 0.

	rdati = reshape(dati,siz1)
	area = new((/NLAT,NLON/),typeof(xf))

	do k = 0, siz1(0)-1
		area = 0.
		do n = 0, nxgrd-1
			if (rdati(k,pxj(n),pxi(n)) .ne. missing) then
				rdato(k,rxj(n),rxi(n)) = rdato(k,rxj(n),rxi(n)) + rdati(k,pxj(n),pxi(n))*xf(n)
				area(rxj(n),rxi(n)) =  area(rxj(n),rxi(n)) + xf(n)
			end if
		end do
		do i = 0, NLON-1
			do j = 0, NLAT-1
				if (area(j,i).gt.0.) then
					rdato(k,j,i) = rdato(k,j,i)/area(j,i)
				else
					rdato(k,j,i) = missing
				end if
			end do
		end do
	end do

	dato = reshape(rdato,osiz)
	dato@missing_value = missing
	fo->\$vnm\$ = dato

	delete(dati)
	delete(rdati)
	delete(dato)
	delete(rdato)
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


