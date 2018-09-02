
nlon = 360
nlat = 180

emisref=(/ 0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99 /)

dat1 = asciiread("emis.txt",(/nlat,nlon/),"integer")
dat = dat1(::-1,:)

emis = new((/nlat,nlon/),"float")

do i = 0, nlon - 1
	do j = 0, nlat - 1 
		idx = max((/2,dat(j,i)/))
		if ( idx .ge. 7 ) then
			idx = 2
		end if
		idx = idx - 1
		emis(j,i) = emisref(idx)
	end do
end do
		
glon = lonGlobeF(nlon,"lon","longitude","degrees_east")
glat = latGlobeF(nlat,"lat","latitude","degrees_north")

emis!0 = "lat"
emis!1 = "lon"

emis&lon = glon
emis&lat = glat

system("rm -f emis_1x1.nc")

fo=addfile("emis_1x1.nc","c")

fo->emis=emis

