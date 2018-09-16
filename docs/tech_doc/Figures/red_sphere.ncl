begin

nlat=94
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

print(latb)

x = new((/4,npoly/),double)
y = new((/4,npoly/),double)
xc = new((/npoly/),double)
yc = new((/npoly/),double)

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


;************************************************
; create the plot
;************************************************
 ;wks = gsn_open_wks ("png","polyg")     ; send graphics to PNG file
 wks = gsn_open_wks ("x11","polyg")     ; send graphics to PNG file
 ;proj = "Hammer" 
 proj = "Satellite" 
 ;proj = "Orthographic" 
 gsn_define_colormap(wks,"default")     ; set color map

 res=True
 res@mpProjection = proj
 res@mpCenterLonF               = 0.0       ; choose center lon
 res@mpCenterLatF               = 30.         ; choose center lat
 ;res@mpSatelliteDistF           = 3.0         ; choose satellite view 
; map resources
 res@mpDataBaseVersion   = "MediumRes"  ; change database
 res@mpFillOn            = True         ; fill continents
 res@mpPerimOn           = False        ; no box around map

 res@mpGridAndLimbOn        = False                ; no lat/lon lines
   
; increase viewport size (make map bigger)
 res@vpWidthF  = 0.93
 res@vpHeightF = 0.93
 res@vpXF      = 0.05
 res@vpYF      = 0.99   
      
 res@gsnDraw  = False       ; don't draw the plots now
 res@gsnFrame = False       ; or advance the frame

 plot = gsn_map(wks,proj,res) ; create the plot   

; polygon resources
  res_poly             = True
  res_poly@gsEdgesOn   = True       ; draw border around polygons
  res_poly@gsEdgeColor = "black"    
  res_poly@gsFillColor = (/0,0,0,0/)        ; VA clim div 1

 restxt=True
 restxt@txFontHeightF = 0.02
  do n = 0, npoly-1
     gsn_polygon (wks,plot,x(:,n),y(:,n),res_poly)
     if (xc(n).eq.0.) then
     	;gsn_text(wks,plot,sprintf("%.2f",yc(n)),xc(n),yc(n),restxt)
     end if
  end do
  frame(wks)    
end
