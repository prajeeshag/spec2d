
begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Add filled boxes to an existing plot based on a range of levels
;; and their associated colors.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("add_filled_boxes")
procedure add_filled_boxes(wks,plot,data,levels,colors)
local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, index
begin
        dims = dimsizes(data)
        nrows = dims(0)
        ncols = dims(1)
    	gnid = new((/nrows,ncols/),"graphic")
        gnres = True     ; resource list for filled polygons
        do jj = 0,nrows-1
          do ii = 0,ncols-1
                boxx = (/ii,ii,ii+1,ii+1,ii/)
                boxy = (/jj,jj+1,jj+1,jj,jj/)
        		index = closest_val(data(jj,ii),levels)
                gnres@gsFillColor = colors(index,:)
				gnres@gsEdgesOn = True
				gnres@gsEdgeThicknessF = 2.
			    gnres@gsEdgeColor = "white"	
                gnid(jj,ii) = gsn_add_polygon(wks,plot,boxx,boxy,gnres)
        	end do
      	end do
      plot@gnid = gnid
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Add filled boxes to an existing plot based on a range of levels
;; and their associated colors.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("add_filled_boxes1")
procedure add_filled_boxes1(wks,plot,data,levels,colors)
local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, index
begin
        dims = dimsizes(data)
		ngrids = num(data.ne.0)
        nrows = dims(0)
        ncols = dims(1)
    	gnid = new((/ngrids/),"graphic")
        gnres = True     
	    n = -1
        do jj = 0,nrows-1
		  ncols0 = num(data(jj,:).ne.0)
          do ii = 0,ncols0-1
				fcols=tofloat(ncols)/tofloat(ncols0)
				n = n + 1
				i = ii*fcols
                ip1 = (ii+1)*fcols
                boxx = (/i,i,ip1,ip1,i/)
                boxy = (/jj,jj+1,jj+1,jj,jj/)
        		index = closest_val(data(jj,ii),levels)
                gnres@gsFillColor = colors(index,:)
				gnres@gsEdgesOn = True
				gnres@gsEdgeThicknessF = 2.
			    gnres@gsEdgeColor = "white"	
                gnid(n) = gsn_add_polygon(wks,plot,boxx,boxy,gnres)
        	end do
      	end do
      plot@gnid = gnid
end

undef("add_text_strings1")
procedure add_text_strings1(wks,plot,data,fonth)
local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, index
begin
        dims = dimsizes(data)
		ngrids = num(data.ne.0)
        nrows = dims(0)
        ncols = dims(1)
    	xpos = new((/ngrids/),"float")
    	ypos = new((/ngrids/),"float")
    	data0 = new((/ngrids/),"integer")
	    n = -1
      	txres               = True
      	txres@txFontHeightF = fonth
      	txres@txJust        = "CenterCenter"
        do jj = 0,nrows-1
		  ncols0 = num(data(jj,:).ne.0)
          do ii = 0,ncols0-1
				fcols=tofloat(ncols)/tofloat(ncols0)
				n = n + 1
				i = ii*fcols
                ip1 = (ii+1)*fcols
				ypos(n) = (i+ip1)/2
				xpos(n) = jj + 0.5
				data0(n) = data(jj,ii)
        	end do
      	end do
        plot@txid = gsn_add_text(wks,plot,sprintf("%2.2g",data0),ypos,xpos,txres)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Add text strings to every filled box showing the data value.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("add_text_strings")
procedure add_text_strings(wks,plot,data,add_offset,fonth)
local txres, xpos, ypos, dims, nrows, cols, add_offset
begin
      dims = dimsizes(data)
      nrows = dims(0)
      ncols = dims(1)
      txres               = True
      txres@txFontHeightF = fonth
      txres@txJust        = "CenterCenter"
      xpos = ndtooned(conform_dims((/nrows,ncols/),ispan(0,nrows-1,1),0) + 0.5)
      ypos = ndtooned(conform_dims((/nrows,ncols/),ispan(0,ncols-1,1),1) + 0.5)
	  data1d = ndtooned(data)
	  ind0 = ind(data1d.ne.0)
	  xpos1 = xpos(ind0)
	  ypos1 = ypos(ind0)
      data0 = data1d(ind0) 
      plot@txid = gsn_add_text(wks,plot,sprintf("%2.2g",data0+add_offset),ypos1,xpos1,txres)
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Attach a horizontal labelbar to the bottom of a plot.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("add_labelbar1")
procedure add_labelbar1(wks,plot,colors,labels)
local nboxes, vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  nboxes = dimsizes(colors)

  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False         ; Necessary to control sizes
  lbres@lbFillColors       = colors        ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.015          ; font height. Default is small
  lbres@lbLabelAlignment   = "BoxCenters"
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False

  lbres@vpWidthF           = vpw * 0.95   ; Make labelbar slightly shorter than width of plot
  lbres@vpHeightF          = vph * 0.25

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
;
; Set some annotation resources indicating how we want to attach
; the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to
; indicate where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;
  amres                  = True
  ;amres@amJust           = "BottomCenter"
  amres@amZone = 2 
  amres@amParallelPosF   = 0.5     ; Centered about X axis
  amres@amOrthogonalPosF = 0.1     ; Move labelbar down


  plot@annoid = gsn_add_annotation(plot,lbid,amres)
end



;********************************************************************************	
;GLOBAL SECTION
NPLON=4
NLAT = 12
NLON = NLAT*2+(NPLON-4)

ifpack="2".eq.2
reduce="1".eq.1

maskfile=""

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

YUNITS=(/"degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN", 	         "degrees_south", "degree_south", "degree_S", "degrees_S", "degreeS", "degreesS"/)
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

varlist = rm_single_dims(str_split_csv(" ",",",0))

optlist = rm_single_dims(str_split_csv(" ",":",0))

NOPT = dimsizes(optlist)

OPTDIC = new((/2,NOPT/),string)

do i = 0, NOPT-1
	OPTDIC(:,i) = rm_single_dims(str_split_csv(optlist(i),"=",0))
end do

TRUNC=10

regd=new((/NLAT,NLON,2/),integer)
redd=new((/NLAT,NLON,2/),integer)
ocd=new((/OCNY,OCNX,2/),integer)
fourd=new((/NLAT,TRUNC,2/),integer)
sphd=new((/TRUNC,TRUNC,2/),integer)

regd=0
redd=0
ocd=0
fourd=0
sphd=0
sphd(:,:,0) = -1


sfX=new((/NLON+1/),float)
sfY=new((/NLAT+1/),float)
sfX_four=new((/TRUNC+1/),float)

sfX_oc=new((/OCNX+1/),float)
sfY_oc=new((/OCNY+1/),float)


do j = 0, OCNY-1
	n = 0
	do i = 0, OCNX-1
		if (i.ge.ILEN(0,j)) then
			n = 1
		end if 
		jj = PACK(n,j)+1
		if (jj.le.NLAT/2) then
			jj =  jj - 1 -NLAT/2
		else
			jj = jj - NLAT/2
		end if
		ocd(j,i,0) = jj
		ocd(j,i,1) = i-ILEN(0,j)*n+1
	end do
end do

n = 1
i1 = 0
p1 = 0
do i = 0, NLON-1
	do j = 0, NLAT-1
		jj = j+1
		if (jj.le.NLAT/2) then
			jj =  jj - 1 -NLAT/2
		else
			jj = jj - NLAT/2
		end if
		regd(j,i,0) = jj
		regd(j,i,1) = i+1
		if (i.lt.LONSPERLAT(j)) then
			redd(j,i,0) = jj
			redd(j,i,1) = i+1
		end if
		if (i.lt.TRUNC) then
			fourd(j,i,0) = jj
			fourd(j,i,1) = i+1
		end if
		if ((i.lt.TRUNC).and.(j.ge.i).and.(j.lt.TRUNC)) then
			sphd(j,i,0) = i
			sphd(j,i,1) = n
			n = n + 1
		end if
	end do
end do

do i = 0, NLON
	do j = 0, NLAT
		sfY(j) = j
		sfX(i) = i
		if (i.le.TRUNC) then
			sfX_four(i) = i
		end if
	end do
end do


print(max(ocd))

do i = 0, OCNX
	do j = 0, OCNY
		sfY_oc(j) = j
		sfX_oc(i) = i
	end do
end do


wks = gsn_open_wks("eps","amfi_grid_table")


res                      = True
res@gsnMaximize          = True          ; Maximize plot in frame.
res@gsnDraw     = False
res@gsnFrame    = False
res@gsnScale    = True
res@gsnPaperOrientation = "auto"
res@tmXBMode = "Explicit"
res@tmYLMode = "Explicit"
res@tmXBMajorOutwardLengthF = 0.
res@tmXBMajorLengthF        = 0.
res@tmYLMajorOutwardLengthF = 0.
res@tmYLMajorLengthF        = 0.
res@tmXTOn                  = False
res@tmYROn                  = False


;-> G-grid regular and reduced
levels = ispan(-NLAT/2,NLAT/2,1)
rgba   = read_colormap_file("matlab_jet")
colors = span_color_rgba(rgba,dimsizes(levels))
colors(NLAT/2,0:2) = .3

res@vpWidthF = 0.8
res@vpHeightF = 0.8*tofloat(NLAT)/tofloat(NLON)
res@tmXBValues := sfX(::7)
res@tmYLValues := sfY(::8)

res@tmXMajorGrid          = True
res@tmXMajorGridLineColor = "black"
res@tmYMajorGrid          = True
res@tmYMajorGridLineColor = "black"

res@tmXMajorGridThicknessF = 3.0
res@tmYMajorGridThicknessF = 3.0

res@trXMinF                 = 0
res@trXMaxF                 = NLON
res@trYMinF                 = 0
res@trYMaxF                 = NLAT

plot1 = gsn_csm_blank_plot(wks,res)
add_filled_boxes(wks,plot1,regd(:,:,0),levels,colors)
add_text_strings(wks,plot1,regd(:,:,1),0,0.008)
add_labelbar1(wks,plot1,colors,""+levels)
draw(plot1)
frame(wks)

plot2 = gsn_csm_blank_plot(wks,res)
add_filled_boxes1(wks,plot2,redd(:,:,0),levels,colors)
add_text_strings1(wks,plot2,redd(:,:,1),0.008)
add_labelbar1(wks,plot2,colors,""+levels)
draw(plot2)
frame(wks)

plot3 = gsn_csm_blank_plot(wks,res)
add_filled_boxes(wks,plot3,redd(:,:,0),levels,colors)
add_text_strings(wks,plot3,redd(:,:,1),0,0.008)
add_labelbar1(wks,plot3,colors,""+levels)
draw(plot3)
frame(wks)

; Fourier domain
levels = ispan(-NLAT/2,NLAT/2,1)
rgba   = read_colormap_file("matlab_jet")
colors = span_color_rgba(rgba,dimsizes(levels))
colors(NLAT/2,0:2) = .3

res@vpWidthF = 0.8*tofloat(TRUNC)/tofloat(NLAT)
res@vpHeightF = 0.8
res@tmXBValues := (/0,2,4,6,8,10/)
res@tmYLValues := sfY(::8)

res@tmXMajorGrid          = True
res@tmXMajorGridLineColor = "black"
res@tmYMajorGrid          = True
res@tmYMajorGridLineColor = "black"

res@tmXMajorGridThicknessF = 3.0
res@tmYMajorGridThicknessF = 3.0

res@trXMinF                 = 0
res@trXMaxF                 = TRUNC
res@trYMinF                 = 0
res@trYMaxF                 = NLAT

plot4 = gsn_csm_blank_plot(wks,res)
add_filled_boxes(wks,plot4,fourd(:,:,0),levels,colors)
add_text_strings(wks,plot4,fourd(:,:,1),-1,0.012)
draw(plot4)
frame(wks)
plot4l = gsn_csm_blank_plot(wks,False)
add_labelbar1(wks,plot4l,colors,""+levels)
draw(plot4l)
frame(wks)

; spectral domain
levels := ispan(-1,TRUNC-1,1)
rgba   = read_colormap_file("matlab_jet")
colors := span_color_rgba(rgba,dimsizes(levels))

res@vpWidthF = 0.8
res@vpHeightF = 0.8
res@tmXBValues := (/0,2,4,6,8,10/)
res@tmYLValues := (/0,6,10/)

res@tmXMajorGrid          = True
res@tmXMajorGridLineColor = "black"
res@tmYMajorGrid          = True
res@tmYMajorGridLineColor = "black"

res@tmXMajorGridThicknessF = 3.0
res@tmYMajorGridThicknessF = 3.0

res@trXMinF                 = 0
res@trXMaxF                 = TRUNC
res@trYMinF                 = 0
res@trYMaxF                 = TRUNC

plot5 = gsn_csm_blank_plot(wks,res)
add_filled_boxes(wks,plot5,sphd(:,:,0),levels,colors)
add_text_strings(wks,plot5,sphd(:,:,1),-1,0.015)
draw(plot5)
frame(wks)

plot5l = gsn_csm_blank_plot(wks,False)
add_labelbar1(wks,plot5l,colors,""+levels)
draw(plot5l)
frame(wks)


;P-grid
levels := ispan(-NLAT/2,NLAT/2,1)
rgba   := read_colormap_file("matlab_jet")
colors := span_color_rgba(rgba,dimsizes(levels))
colors(NLAT/2,0:2) = .3

res@vpWidthF = 0.8
res@vpHeightF = 0.8*tofloat(OCNY)/tofloat(OCNX)
res@tmXBValues := sfX_oc(::6)
res@tmYLValues := sfY_oc(::4)

res@tmXMajorGrid          = True
res@tmXMajorGridLineColor = "black"
res@tmYMajorGrid          = True
res@tmYMajorGridLineColor = "black"

res@tmXMajorGridThicknessF = 2.0
res@tmYMajorGridThicknessF = 2.0

res@trXMinF                 = 0
res@trXMaxF                 = OCNX
res@trYMinF                 = 0
res@trYMaxF                 = OCNY

plot6 = gsn_csm_blank_plot(wks,res)
add_filled_boxes(wks,plot6,ocd(:,:,0),levels,colors)
add_text_strings(wks,plot6,ocd(:,:,1),0,0.008)
add_labelbar1(wks,plot6,colors,""+levels)
draw(plot6)
frame(wks)


end
