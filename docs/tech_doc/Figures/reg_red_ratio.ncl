

n=(1024-64+2)/2
nlats = new(n,integer)
npts_red = new(n,integer)
npts_reg = new(n,integer)

n = 0
do nlat = 64, 1024, 2
	maxlon = 20+(nlat/2-1)*4
	ocnx = maxlon+20
	ocny = nlat/2
	nlats(n) = nlat
	npts_red(n)=ocnx*ocny
	npts_reg(n)=maxlon*nlat
	n = n + 1
end do

pred=(npts_reg-npts_red)*100./npts_reg

pred@long_name="% reduction in number of grid points"
nlats@long_name="number of latitude circles"

wks = gsn_open_wks ("eps","pred")
res = True
plot  = gsn_csm_xy (wks, nlats, pred, res)

