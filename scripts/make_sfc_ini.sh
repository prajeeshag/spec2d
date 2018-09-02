dat=gfs/work/atm_3hr_inst_1870_01_01_00.nc

ncks -v snwdph,sheleg,snc,tskin,smsoil,stsoil,slsoil,trans,canopy $dat out1.nc

cdo -ltc,9e20 out1.nc mask.nc

cdo -ifthen mask.nc out1.nc out2.nc
mv out2.nc out1.nc

cdo -fillmiss out1.nc out2.nc
mv out2.nc out1.nc

ncpdq -a lon,lat,lev4 out1.nc out2.nc
mv out2.nc out1.nc

ncpdq -a -lat out1.nc out2.nc
mv out2.nc out1.nc

ncrename -d time,Time -v snc,sncovr -v tskin,tslnd -v stsoil,stc -v slsoil,slc -v smsoil,smc out1.nc out2.nc
mv out2.nc out1.nc

#ncwa -O -a time -d time,0,0 out1.nc out2.nc
#mv out2.nc out1.nc

mv out1.nc sfc_res.nc


