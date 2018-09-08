
OCSNF=/moes/home/prajeesh/spec2d/octa2reg/octaStackNFold.sh
ocnx=204
ocny=94

EMIS=/moes/home/prajeesh/spec2d/data/global_emis.1x1.nc 
MTN=/moes/home/prajeesh/spec2d/data/global_mtn.nc
OZONE=/moes/home/prajeesh/spec2d/data/o3clim.nc
TG3=/moes/home/prajeesh/spec2d/data/global_tg3clim.nc
VEGFRAC=/moes/home/prajeesh/spec2d/data/global_vegfrac.1x1.nc
ALBEDO=/moes/home/prajeesh/spec2d/data/global_albedo4.1x1.nc
ZORL=/moes/home/prajeesh/spec2d/data/global_zorclim.1x1.nc
SLMSK=/moes/home/prajeesh/spec2d/data/global_slmask.nc
SLOPETYPE=/moes/home/prajeesh/spec2d/data/global_slope.1x1.nc
SOILTYPE=/moes/home/prajeesh/spec2d/data/global_soiltype.1x1.nc
VEGTYPE=/moes/home/prajeesh/spec2d/data/global_vegtype.1x1.nc

SST=/iitm2/cccr-res/prajeesh/observations/HadISST/sst_monthly.nc
ICE=/iitm2/cccr-res/prajeesh/observations/HadISST/ice_monthly.nc



OUTDIR=${ocnx}x${ocny}

mkdir -p $OUTDIR

#$OCSNF -x $ocnx -y $ocny -i $EMIS     -o $OUTDIR/emis_ref.nc -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $MTN      -o $OUTDIR/mtn.nc      -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $OZONE    -o $OUTDIR/ozone.nc    -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $TG3      -o $OUTDIR/tg3.nc      -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $VEGFRAC  -o $OUTDIR/vegfrac.nc  -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $ALBEDO   -o $OUTDIR/albedo.nc   -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $ZORL     -o $OUTDIR/zorl.nc     -p interpmethod=linint

#$OCSNF -x $ocnx -y $ocny -i $ICE -o $OUTDIR/ice_monthly.nc   -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i $SST -o $OUTDIR/sst_monthly.nc   -p interpmethod=linint


#$OCSNF -x $ocnx -y $ocny -i $SLMSK -o _tem.nc -v AREA_LND -p interpmethod=conserve
#cdo -gec,0.5 _tem.nc $OUTDIR/grid_spec.nc
#rm -f _tem.nc
#cdo -chname,AREA_LND,mask $OUTDIR/grid_spec.nc $OUTDIR/mask.nc
#cdo -chname,AREA_LND,AREA_LND_CELL -addc,1. -mulc,0. $OUTDIR/grid_spec.nc AREA_LND_CELL.nc
#ncks -A AREA_LND_CELL.nc $OUTDIR/grid_spec.nc
#rm -f AREA_LND_CELL.nc
#
#$OCSNF -x $ocnx -y $ocny -i $SLOPETYPE -m $OUTDIR/mask.nc -o $OUTDIR/slopetype.nc -p interpmethod=dtype:conv2int=T
#$OCSNF -x $ocnx -y $ocny -i $SOILTYPE -m $OUTDIR/mask.nc -o $OUTDIR/soiltype.nc -p interpmethod=dtype:conv2int=T
#$OCSNF -x $ocnx -y $ocny -i $VEGTYPE -m $OUTDIR/mask.nc -o $OUTDIR/vegtype.nc -p interpmethod=dtype:conv2int=T


#$OCSNF -x $ocnx -y $ocny -i gt_aero_time.nc -o octa/gt_aero_time.nc -p interpmethod=linint
#$OCSNF -x $ocnx -y $ocny -i ../work1/INPUT/amfi_res.nc  -o octa/amfi_res.nc -p interpmethod=conserve \
#	-v canopy,sheleg,slc,smc,sncovr,snwdph,stc,trans,tslnd,psp,psp1,qp,qp1,tp,tp1

