#!/bin/bash 

usage() { echo "Usage: $0 -x NLON -y NLAT -i inputfile.nc [-n] [-r] \
[-o outfile.nc] [-m grid_spec_file] [-v vars] [-p optionlist]" 1>&2; exit 1;}

npack=2
reduce=1

outfile=""

while getopts 'x:y:i:o:v:p:m:nr' flag; do
    case "${flag}" in
    x) NLON="$OPTARG" ;;
    y) NLAT="$OPTARG" ;;
    i) infile="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
	v) valist="$OPTARG" ;;
	p) oplist="$OPTARG" ;;
	m) gridspec="$OPTARG" ;;
	n) npack=1 ;;
	r) reduce=0 ;;
    *)
		usage	
		;;
    esac
done


if [[ -z $NLON ]] || \
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

cat <<EOF > $tfile

begin

;********************************************************************************	
;GLOBAL SECTION
NLON = $NLON
NLAT = $NLAT

ifpack="$npack".eq.2
reduce="$reduce".eq.1

maskfile="$gridspec"

ifomask = True
if (maskfile.eq."") then
	ifomask=False
end if

NPACK=1
if (ifpack) then
	NPACK=2
end if

if (mod(NLAT,2).ne.0) then
	print("FATAL: NLAT should be a multiple of 2")
	status_exit(1)
end if	

YUNITS=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", \
	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
XUNITS=(/"degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE"/)

LONSPERLAT=new((/NLAT/),integer)

NPLON = NLON - (NLAT/2-1)*4

LONSPERLAT = NLON

if (reduce) then
	LONSPERLAT(0) = NPLON
	LONSPERLAT(NLAT-1) = NPLON
	
	do i = 1, NLAT/2-1
		LONSPERLAT(i) = LONSPERLAT(i-1) + 4
		LONSPERLAT(NLAT-1-i) = LONSPERLAT(i)
	end do
else
	LONSPERLAT = NLON
	NPLON = NLON
end if

NLAT@double = True
NLON@double = True
LONF = lonGlobeF(NLON, "lon", "longitude", "degrees_E")
LATF = latGau(NLAT, "lat", "latitude", "degrees_N")

OCNX = (NPACK-1)*NPLON+max(LONSPERLAT)
OCNY = NLAT/NPACK

print("OCNX="+OCNX)
print("OCNY="+OCNY)

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
	status_exit(1)
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
		print("FATAL: Cannot determine cartesian axis attribute for "+ax)
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

	if (ynm.eq."".or.xnm.eq."") then
		print("FATAL: Cannot find X or Y coordinate for variable: "+varnm)
		status_exit(1)
	end if

	axnm = (/ynm,xnm/)
	axnm@z=False
	if (znm.ne."") then
		axnm1 = axnm
		delete(axnm)
		axnm=new(dimsizes(axnm1)+1,typeof(axnm1))
		axnm(0) = znm
		axnm(1:) = axnm1
		axnm@z=True
		delete(axnm1)
	end if
	
	axnm@t=False
	if (tnm.ne."") then
		axnm1 = axnm
		delete(axnm)
		axnm=new(dimsizes(axnm1)+1,typeof(axnm1))
		axnm(0) = tnm
		axnm(1:) = axnm1
		axnm@t=True
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
		print("FATAL: packed field contains missing values!!!!")
		status_exit(1)
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
		status_exit(1) 
	end if
	
	siz = dimsizes(fi)
	ndim = dimsizes(siz)

	if (ndim.lt.2) then
		print("FATAL: regrid_discrete: fi should be atleast a 2D field")
		status_exit(1)
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
		status_exit(1)
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

;--------------------------------------------------------------------------------	
function rval(dat:numeric,i:integer,r:integer)
local dat, i, r, n
begin
;--------------------------------------------------------------------------------	

	n = i + r 
	if (n.ge.dimsizes(dat)) then
		n = n-dimsizes(dat)
	end if
	return(dat(n))
end

;--------------------------------------------------------------------------------	
function lval(dat:numeric,i:integer,r:integer)
local dat, i, r, n
begin
;--------------------------------------------------------------------------------	

	n = i - r 
	if (n.lt.0) then
		n = dimsizes(dat)+n
	end if
	return(dat(n))
end

;--------------------------------------------------------------------------------	
function search_and_fill(dat:numeric,mask:numeric)
local dat, mask, k, i, r, val, miss
begin
;--------------------------------------------------------------------------------

	siz = dimsizes(dat)
	nk = siz(0)
	miss = new(siz,logical)
	
	do k = 0, nk - 1
		miss = (dat(k,:)-mask).lt.0.
	end do
	
    do k = 0, nk - 1
		do i = 0, siz(1) - 1
			if (.not.miss(k,i)) then
				continue	
			end if
			do r = 1, siz(1)/2
				val=max((/lval(dat(k,:),i,r),rval(dat(k,:),i,r)/))
				if (val.gt.0) then
					dat(k,i) = val
					continue
				end if
			end do
		end do
	end do

	return(dat)
end

;-------------------------------------------------------------------------------
function unstack_and_unfold_conserve(dati:numeric)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd
begin
;--------------------------------------------------------------------------------	

	siz = dimsizes(dati)
	ndim = dimsizes(siz)

	if (siz(ndim-1).ne.OCNX.or.siz(ndim-2).ne.OCNY) then
		print("FATAL: unstack_and_unfold_conserve: dati size mismatch")
		status_exit(1)
	end if

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

;--------------------------------------------------------------------------------	
function winner(dati:numeric)
local dati, val, i, n, nmax, ii
begin
;--------------------------------------------------------------------------------	

	unq = get_unique_values(dati)
	if(dimsizes(unq).eq.1) then
		return(unq(0))
	end if

	nmax = 0	
	do i = 0, dimsizes(unq)-1
		n = num(dati.eq.unq(i))
		if (n.gt.nmax) then
			nmax = n
			ii = i
		end if
	end do
	
	return(unq(ii))

end 

;--------------------------------------------------------------------------------	
function assign_nearest1(dati:numeric, prblm:logical)
local dati, indz, dato, i, j, n, siz, mask, siz1, i1, j1, i2, j2, lat1, lon1
begin
;--------------------------------------------------------------------------------


	indz = ind_resolve(ind(ndtooned(prblm)),dimsizes(prblm))
	siz = dimsizes(indz)
	siz1 = dimsizes(dati)
	wnsiz = min(siz1)/2
	ilen = siz1(1)
	jlen = siz1(0)

	iar = new(siz1,integer)
	jar = new(siz1,integer)

	do i = 0, ilen-1
		do j = 0, jlen-1
			iar(j,i) = i
			jar(j,i) = j
		end do
	end do

	iar1d = ndtooned(iar)
	jar1d = ndtooned(jar)
	dat1d = ndtooned(dati)
	dato = dati

	do n = 0, siz(0)-1
		j1 = indz(n,0)
		i1 = indz(n,1)
	   	
		do w = 1, wnsiz-1	
			is = i1 - w
			ie = i1 + w
			js = j1 - w
			je = j1 + w
			if (is.lt.0) then
				is=ilen+is
			end if
			if (ie.ge.ilen) then
				ie=ie-ilen
			end if
			if (js.lt.0) then
				js = 0
			end if
			if (je.ge.jlen) then
				je = jlen-1
			end if

			if (is.gt.i1.or.ie.lt.i1) then
				ia = (iar1d.ge.is.or.iar1d.le.ie)
			else 
				ia = (iar1d.ge.is.and.iar1d.le.ie)
			end if

			ja = (jar1d.ge.js.and.jar1d.le.je)
			
			wvals = dat1d(ind(ia.and.ja))
			ind0 = ind(wvals.gt.0)
		
			if (any(ismissing(ind0))) then
				delete(wvals)
				delete(ind0)
				continue
			end if

			wval0 = wvals(ind(wvals.gt.0))

			dato(j1,i1) = winner(wval0)
		
			delete(wvals)
			delete(wval0)
			delete(ind0)
			break
		end do
	end do

	return(dato)

end 

;--------------------------------------------------------------------------------	
procedure assign_nearest2(dati:numeric, indz:numeric)
local dati, indz, dato, i, j, n, siz, mask, siz1, i1, j1, i2, j2, lat1, lon1
begin
;--------------------------------------------------------------------------------	
	siz = dimsizes(indz)
	siz1 = dimsizes(dati)
	pole = new(1,typeof(OCLAT))
	pole = 90.	

	do n = 0, siz(0)-1
		j1 = indz(n,0)
		i1 = indz(n,1)
		lat1 = OCLAT(j1,i1)
		lon1 = OCLON(j1,i1)
		min_dist = gc_latlon(-pole,0.,pole,0.,2,4)
		i2 = -1
		j2 = -1
		do i = 0, siz1(1)-1
			do j = 0, siz1(0)-1
				if (dati(j,i).lt.1) then
					continue
				end if
			    dist = gc_latlon(lat1,lon1,OCLAT(j,i),OCLON(j,i),2,4)
				if (dist.lt.min_dist) then
					min_dist = dist
					i2 = i
					j2 = j
				end if
			end do
		end do
		if (i2.eq.-1.or.j2.eq.-1) then
			print("FATAL: i2.eq.-1.or.j2.eq.-1")
			status_exit(1)
		end if
		print(" "+n+" "+dati(j2,i2)+" "+dati(j1,i1))
		dati(j1,i1) = dati(j2,i2)
	end do
end 


;-------------------------------------------------------------------------------
function stack_and_fold_dtype(fi:numeric,xi:numeric,yi:numeric,ongrid:string,conv2int:logical)
local siz, rsiz, dato, datii, ndim, i, datf, xi, yi, j, k, ongrid, ongrd, conv2int, fi, dati
begin
;--------------------------------------------------------------------------------

	if (.not.isinteger(fi)) then
		if (.not.conv2int) then
			print("FATAL: input array to stack_and_fold_dtype should be of integer type")
			status_exit(1)
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

	if (ifomask) then
		grid_spec=addfile(maskfile,"r")
		omask = grid_spec->AREA_LND
		if (any(dimsizes(omask)-(/OCNY,OCNX/).ne.0)) then
			if (all(dimsizes(omask)-(/OCNX,OCNY/).eq.0)) then
				otmp=transpose(omask)
				delete(omask)
				omask=otmp
			else
				print("FATAL: AREA_LND should be of the size [OCNY,OCNX]")
				status_exit(1)
			end if
		end if
	else
		omask = new((/OCNY,OCNX/),float)
		omask = 0. 
	end if

	omask = where(omask.gt.0.,1.,0.)

	rmask = unstack_and_unfold_conserve(omask)

	rmask = where(rmask.gt.0.,1.,0.)

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

	do k = 0, howmany-1
		prblm = (datii(k,:,:) - rmask).lt.0.
		nmprbl = num(prblm)
		print("mask problem before nearest assign1 - "+nmprbl)
		if (nmprbl.gt.0) then
			datii(k,:,:) = assign_nearest1 (datii(k,:,:),prblm)
			prblm = (datii(k,:,:) - rmask).lt.0.
			nmprbl = num(prblm)
			print("mask problem after nearest assign1 - "+nmprbl)
		end if
	end do
	delete(prblm)

	ncopy = 4
	rfi = new((/howmany,ncopy,nx/),typeof(datii))
	latfi = latGau(ncopy, "lat", "latitude", "degrees_N")

    do j = 0, NPACK-1
        do i = 0, OCNY-1
            lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
			ip = PACK(j,i)
			do k = 0, ncopy-1
				rfi(:,k,:) = datii(:,ip,:)
			end do
            rfo = regrid_discrete_Wrap(LONF,latfi,rfi,lonc,latfi,False)
            datoi(:,i,IS(j,i):IE(j,i)) = (/rfo(:,(ncopy+1)/2-1,:)/)
			delete(lonc)
			delete(rfo)
        end do
    end do

	do k = 0, howmany-1
		prblm = (datoi(k,:,:)-omask).lt.0
		nmprbl=num(prblm)
		print("mask problem before nearest assign2 - "+nmprbl)
		if (nmprbl.gt.0) then
			indz = ind_resolve(ind(ndtooned(prblm)),dimsizes(prblm))
			print(dimsizes(indz))
			assign_nearest2(datoi(k,:,:), indz)
			delete(indz)
			;status_exit(1) 
			prblm = datoi(k,:,:)-omask.lt.0
			nmprbl=num(prblm)
			print("mask problem after nearest assign2 - "+nmprbl)
		end if
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
function stack_and_fold_conserve(datf:numeric,xi:numeric,yi:numeric,ongrid:string)
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
    	dati = area_conserve_remap_Wrap(xxi, yyi, datf, LONF, LATF, False)
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
	
	ncopy = 4
	rfi = new((/howmany,ncopy,nx/),typeof(datii))
	latfi = latGau(ncopy, "lat", "latitude", "degrees_N")

    do j = 0, NPACK-1
        do i = 0, OCNY-1
            lonc = lonGlobeF(ILEN(j,i), "lon", "longitude", "degrees_E")
			ip = PACK(j,i)
			do k = 0, ncopy-1
				rfi(:,k,:) = datii(:,ip,:)
			end do
            rfo = area_conserve_remap(LONF,latfi,rfi,lonc,latfi,False)
			do k = 0, 1
            	datoi(:,i,IS(j,i):IE(j,i)) = (/rfo(:,(ncopy+1)/2-1,:)/)
			end do
			delete(lonc)
			delete(rfo)
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
;end function stack_and_fold_conserve


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
	
    do j = 0, NPACK-1
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
procedure toOcta_and_write(fo:file, fi:file, vnm:string, xynm[*]:string, intmethd:string, \
							ongrid:string, conv2int:logical)
local i, lonin, latin, NLATin, NLONin, dati, fo, fi, vnm, xynm, intmethd, ndim, siz, conv2int
begin
;--------------------------------------------------------------------------------	

	ndim = dimsizes(xynm) 
	xi = fi->\$xynm(ndim-1)\$
	yi = fi->\$xynm(ndim-2)\$

	rev=1
	if ((yi(1)-yi(0)).lt.0) then
		rev=-1
		print("Input data N->S, converting to S->N")
	end if

	if (ndim.eq.3) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|::rev,\$xynm(2)\$|:)
	else if (ndim.eq.4) then
		dati = fi->\$vnm\$(\$xynm(0)\$|:,\$xynm(1)\$|:,\$xynm(2)\$|::rev,\$xynm(3)\$|:)
	else if (ndim.eq.2) then
		dati = fi->\$vnm\$(\$xynm(0)\$|::rev,\$xynm(1)\$|:)
	end if
	end if
	end if

	siz = dimsizes(dati)
	if (intmethd.eq."unstruct") then
		dato=stack_and_fold_unstrct(dati,xi,yi,ongrid)
	else if (intmethd.eq."dtype") then
		dato=stack_and_fold_dtype(dati,xi,yi,ongrid,conv2int)
	else if (intmethd.eq."linint") then
		dato=stack_and_fold_linint2(dati,xi,yi,ongrid)
	else if (intmethd.eq."conserve") then
		dato=stack_and_fold_conserve(dati,xi,yi,ongrid)
	else
		print("FATAL: invalid interpmethod: "+intmethd)
		print("NOTE: available interpmethods are: 'unstruct', 'linint', 'conserve', and "+\
				"'dtype' (for discrete type interpolation)")
	end if
	end if
	end if
	end if

	if (any(ismissing(dato))) then
		print("FATAL: missing values in dato")
		status_exit(1)
	end if

	xnm = "x"
	ynm = "y" 

	if (xynm@z) then
		znm=xynm(ndim-3)
		dato&\$znm\$@cartesian_axis="Z"
	end if	
				
	if (ndim.gt.3) then
		znm=xynm(ndim-3)
		fo->\$vnm\$ = dato(\$xynm(0)\$|:,\$xnm\$|:,\$ynm\$|:,\$znm\$|:)
	else if (ndim.eq.3) then
		znm=xynm(ndim-3)
		if (xynm@z) then
			fo->\$vnm\$ = dato(\$xnm\$|:,\$ynm\$|:,\$znm\$|:)
		else
			fo->\$vnm\$ = dato(\$xynm(0)\$|:,\$xnm\$|:,\$ynm\$|:)
		end if
	else
		fo->\$vnm\$ = dato(\$xnm\$|:,\$ynm\$|:)
	end if
	end if

	delete(dati)
	delete(dato)
	return
end


;********************************************************************************	
;END FUNCTION SECTION



;********************************************************************************	
;--> Main SECTION

infile="$infile"
fi = addfile(infile,"r")

ofile="$outfile"
if (ofile.eq."") then
	pref="OCx"
	if (ifpack) then
		pref=pref+"Px"
	else if (reduce) then
		pref=pref+"Rx"
	end if
	end if
	ofile=pref+infile
end if

fvnms = getfilevarnames(fi)
recdimdefined=False
fileopened=False
do i = 0, dimsizes(fvnms)-1
	if (any(fvnms(i).eq.getfilevardims(fi,fvnms(i)))) then
		print("Skiping.. "+fvnms(i))
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
			print("Skiping.. "+fvnms(i))
			continue
		end if
	end if
	axnm=find_axis_nms(fvnms(i),fi)
	if (.not.fileopened) then
		system("rm -rf "+ofile)
		fo = addfile(ofile,"c")
		fileopened = True
	end if
	if(axnm@t.and..not.recdimdefined) then
		filedimdef(fo,axnm(0),-1,True)
		recdimdefined=True
	end if
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
ncl -Q $tfile


