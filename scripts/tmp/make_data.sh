
cdo -ymonmean -selyear,1980/2010 HadISST_ice.nc out2.nc
mv out2.nc out1.nc

cdo -remapbil,work/INPUT/tg3.nc -selvar,sic out1.nc out2.nc
mv out2.nc out1.nc

cdo -setmisstoc,0. out1.nc out2.nc
mv out2.nc out1.nc

#ncpdq -a -lat out1.nc out2.nc
#mv out2.nc out1.nc

ncpdq -a lon,lat out1.nc out2.nc
mv out2.nc out1.nc

ncatted -a modulo,time,c,c," " out1.nc

ncatted -a calendar,time,m,c,"julian" out1.nc

mv out1.nc sic_monclim.nc


cdo -ymonmean -selyear,1980/2010 HadISST_sst.nc out2.nc
mv out2.nc out1.nc

cdo -remapbil,work/INPUT/tg3.nc -selvar,sst out1.nc out2.nc
mv out2.nc out1.nc

cdo -fillmiss out1.nc out2.nc
mv out2.nc out1.nc

cdo -addc,273.15 out1.nc out2.nc
mv out2.nc out1.nc

#ncpdq -a -lat out1.nc out2.nc
#mv out2.nc out1.nc

ncpdq -a lon,lat out1.nc out2.nc
mv out2.nc out1.nc

ncatted -a modulo,time,c,c," " out1.nc

ncatted -a calendar,time,m,c,"julian" out1.nc

mv out1.nc sst_monclim.nc


