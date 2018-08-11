#!/bin/bash 

usage() { echo "Usage: $0 -x nlon -y nlat -i inputfile.nc -o outfile.nc [-v varlist]" 1>&2; exit 1;}

while getopts 'x:y:i:o:v:' flag; do
    case "${flag}" in
    x) nlon="$OPTARG" ;;
    y) nlat="$OPTARG" ;;
    i) infile="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
	v) valist="$OPTARG" ;;
    *)
		usage	
		;;
    esac
done


if [[ -z $nlon ]] || [[ -z $outfile ]] || \
       [[ -z $nlat ]] || [[ -z $infile ]]; then
    echo Usage:
    echo $0 -x nlon -y nlat -i inputfile.nc -o outfile.nc 
    exit 1
fi

varlist="\"\""
if ! [[ -z $valist ]]; then
varlist="\"$valist\""
fi

tfile=$(mktemp)

echo $tfile

cat <<EOF > $tfile

nlon = $nlon
nlat = $nlat

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
	ie(0,i) = is(0,i) + lonsperlat(pack(0,i)) - 1
	is(1,i) = ie(0,i) + 1
	ilen(0,i) = lonsperlat(pack(0,i))
	ilen(1,i) = lonsperlat(pack(1,i))
end do

oclon = new((/ocny,ocnx/),typeof(lonf))
oclat = new((/ocny,ocnx/),typeof(latf))

do j = 0, 1
    do i = 0, ocny-1
        lonc = lonGlobeF(ilen(j,i), "lon", "longitude", "degrees_E")
        oclon(i,is(j,i):ie(j,i)) = lonc
        oclat(i,is(j,i):ie(j,i)) = latf(pack(j,i))
        delete(lonc)
    end do
end do

ocnxny = ocny*ocnx

oclon1d = reshape(oclon,(/ocnxny/))
oclat1d = reshape(oclat,(/ocnxny/))

;--------------------------------------------------------------------------------	
function get_cart(ax:string,fid:file)
local i
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
	else
		print("Cannot determine cartesian axis attribute for "+ax)
		exit
	end if
	end if
	return
end

;--------------------------------------------------------------------------------	
function find_axis_nms(varnm:string,fid:file)
local j, xnm, ynm, vdmnm, i
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
function stack_and_fold_unstrct(datf:numeric,xnm,ynm)
local siz, rsiz, dato, datii, ndim, i
begin
;--------------------------------------------------------------------------------	

	siz = dimsizes(datf)
	ndim = dimsizes(siz)

	if (ndim .eq. 2) then
		osiz=(/ocny,ocnx/)
	else if (ndim .eq. 3) then
		osiz=(/siz(0),ocny,ocnx/)
	else if (ndim .eq. 4) then
		osiz = (/siz(0),siz(1),ocny,ocnx/) 
	else
		print("FATAL: ndim cannot be < 2 and > 4")
		exit
	end if
	end if
	end if 
		
	dati = linint2_points_Wrap(datf&\$xnm\$, datf&\$ynm\$, datf, True, oclon1d, oclat1d, 0) 

	dato = reshape(dati,osiz)

	copy_VarCoords_2(dati,dato)

	return(dato)

end

;-------------------------------------------------------------------------------
function stack_and_fold_linint2(datf:numeric,xnm,ynm)
local siz, rsiz, dato, datii, ndim, i
begin
;--------------------------------------------------------------------------------	
    ;dati = linint2_Wrap(datf&\$xnm\$, datf&\$ynm\$, datf, True, lonf, latf, 0)
    dati = linint2_Wrap(datf&\$xnm\$, latf, datf, True, lonf, latf, 0)
    ;dati = area_conserve_remap_Wrap(datf&\$xnm\$, datf&\$ynm\$, datf, lonf, latf, False)

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	if (ndim .eq. 2) then
		howmany=1
		rsiz=siz(0:1)
		osiz=(/ocny,ocnx/)
	else if (ndim .eq. 3) then
		howmany=siz(0)
		rsiz=siz(1:2)
		osiz=(/siz(0),ocny,ocnx/)
	else if (ndim .eq. 4) then
		howmany=siz(0)*siz(1)
		rsiz=siz(2:3)
		osiz = (/siz(0),siz(1),ocny,ocnx/) 
	else
		print("FATAL: ndim cannot be < 2 and > 4")
		exit
	end if
	end if
	end if 
		
	datii = reshape(dati,(/howmany,rsiz(0),rsiz(1)/))
	datoi = new((/howmany,ocny,ocnx/),typeof(dati)) 
	
    do j = 0, 1
        do i = 0, ocny-1
            lonc = lonGlobeF(ilen(j,i), "lon", "longitude", "degrees_E")
            datoi(:,i,is(j,i):ie(j,i)) = linint1(lonf,datii(:,pack(j,i),:),True,lonc,0)
			delete(lonc)
        end do
    end do

	dato = reshape(datoi,osiz)

	copy_VarCoords_2(dati,dato)

	return(dato)

end


;--------------------------------------------------------------------------------	
procedure toOcta_and_write(fo:file, fi:file, vnm:string, xynm[*]:string, intmethd:string)
local axnm, i, lonin, latin, nlatin, nlonin, dati
begin
;--------------------------------------------------------------------------------	

	latin=fi->\$xynm(1)\$
	nlatin=dimsizes(latin)

	lonin=fi->\$xynm(0)\$
	nlonin=dimsizes(lonin)

	if (dimsizes(xynm).eq.3) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:)
		xnm = xynm(2)
		ynm = xynm(1)
	else if (dimsizes(xynm).eq.4) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|:,\$xynm(3)\$|:)
		xnm = xynm(3)
		ynm = xynm(2)
	else if (dimsizes(xynm).eq.2) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:)
		xnm = xynm(1)
		ynm = xynm(0)
	end if
	end if
	end if

	if (intmethd.eq."unstruct") then
		dato=stack_and_fold_unstrct(dati,xnm,ynm)
	else 
		dato=stack_and_fold_linint2(dati,xnm,ynm)
	end if

	fo->\$vnm\$ = dato

	delete(dati)
	delete(dato)
	return
end


varlist = rm_single_dims(str_split_csv($varlist,",",0))
fi = addfile("$infile","r")

system("rm -rf $outfile")
fo = addfile("$outfile","c")

fvnms = getfilevarnames(fi)

do i = 0, dimsizes(fvnms)-1
	if (any(fvnms(i).eq.getfilevardims(fi,fvnms(i)))) then
		print("skipping dim var: " + fvnms(i))
		continue
	end if
	if (any((varlist).ne."")) then
		if (.not.any(fvnms(i).eq.varlist)) then
			print("skipping var: " + fvnms(i))
			continue
		end if
	end if
	axnm=find_axis_nms(fvnms(i),fi)
	toOcta_and_write(fo,fi,fvnms(i),axnm,"linint")	
	delete(axnm)
end do

EOF


#cat $tfile
ncl $tfile
