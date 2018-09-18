begin

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
  lbres@vpHeightF          = vph * 0.08

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
  amres@amZone           = 2
  ;amres@amJust           = "BottomCenter"
  amres@amParallelPosF   = 0.5     ; Centered about X axis
  amres@amOrthogonalPosF = 0.0     ; Move labelbar down

  plot@annoid = gsn_add_annotation(plot,lbid,amres)
end

nlat=64
nplon=20
maxlon = nplon+(nlat/2-1)*4
ocnx = maxlon+nplon
ocny = nlat/2
npoly = ocnx*ocny
PI=get_pi("double")
R2D=get_r2d("double")

lat = latGau(nlat,"lat","lat","degrees_N")
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")

print(lat);

latb = new(nlat+1,double)

latb(0) = -0.5*PI
sumwts = 0.
do j = 0, nlat-2
    sumwts = sumwts + gwt(j)
    latb(j+1) = asin(sumwts-1.)
end do
latb(nlat) = 0.5*PI
latb=latb*R2D

do j = 0, nlat-1
end do

print(latb)

x = new((/4,npoly/),double)
y = new((/4,npoly/),double)
xc = new((/npoly/),double)
yc = new((/npoly/),double)
area = new((/npoly/),double)

n = 0
do j = 1, nlat
  jj = j
  if (jj.gt.nlat/2) then
    jj = nlat - jj + 1
  end if
  nlon = nplon + (jj-1)*4
  lonb = new(nlon+1,typeof(latb))
  dlon = 360./nlon
  lonb(0) = -dlon/2
  do i = 1, nlon
    lonb(i) = lonb(i-1)+dlon
  end do
  do i = 1, nlon
  	 area(n) = sqrt(2*PI*6371*6371*gwt(j-1)/nlon)
     x(0,n) = lonb(i-1) 
     y(0,n) = latb(j-1)
     x(1,n) = lonb(i)
     y(1,n) = latb(j-1)
     x(2,n) = lonb(i)
     y(2,n) = latb(j)
     x(3,n) = lonb(i-1)
     y(3,n) = latb(j)
     xc(n) = (lonb(i-1)+lonb(i))*0.5
     yc(n) = lat(j-1) 
     n = n + 1
  end do
  delete(lonb)
end do

minval=floor(min(area))
maxval=ceil(max(area))
levels=fspan(minval,maxval,nlat/4)
;cmap_name="BlueRed"
;cmap_name="MPL_YlOrRd"
cmap_name="MPL_GnBu"

;************************************************
; create the plot
;************************************************
 ;wks = gsn_open_wks ("png","polyg")     ; send graphics to PNG file
 ;wks = gsn_open_wks ("eps","red_sphere")     ; send graphics to PNG file
 wks = gsn_open_wks ("x11","red_sphere")     ; send graphics to PNG file
 ;proj = "Hammer" 
 proj = "Satellite" 
 ;proj = "Orthographic" 
 gsn_define_colormap(wks,"default")     ; set color map

 res=True
 res@mpProjection = proj
 res@mpCenterLonF               = 80.0       ; choose center lon
 res@mpCenterLatF               = 35.         ; choose center lat
 ;res@mpSatelliteDistF           = 3.0         ; choose satellite view 
; map resources
 res@mpDataBaseVersion   = "MediumRes"  ; change database
 res@mpFillOn            = False         ; fill continents
 res@mpPerimOn           = False        ; no box around map

 res@mpGridAndLimbOn        = False                ; no lat/lon lines
   
; increase viewport size (make map bigger)
; res@vpWidthF  = 0.93
; res@vpHeightF = 0.93
; res@vpXF      = 0.05
; res@vpYF      = 0.99   
      
 res@gsnDraw  = False       ; don't draw the plots now
 res@gsnFrame = False       ; or advance the frame

 plot = gsn_map(wks,proj,res) ; create the plot   

;polygon resources
 res_poly             = True
 res_poly@gsEdgesOn   = True       ; draw border around polygons
 res_poly@gsEdgeColor = "black"    
 ;res_poly@gsFillColor = (/0,0,0,0/)        
 colors = get_color_rgba(cmap_name,levels,area)

 restxt=True
 restxt@txFontHeightF = 0.02
 do n = 0, npoly-1
  	 res_poly@gsFillColor = colors(n,:)
     gsn_polygon (wks,plot,x(:,n),y(:,n),res_poly)
    if (xc(n).eq.0.) then
     	;gsn_text(wks,plot,sprintf("%.2f",yc(n)),xc(n),yc(n),restxt)
     end if
 end do

 colors:= get_color_rgba(cmap_name,levels,levels)
 
 add_labelbar1(wks,plot,colors,sprintf("%.0f",levels))  
 
 draw(plot)
 frame(wks)    
end
