#!/bin/bash

# rootdir -  Absolute path to the root directory
rootdir=/moes/home/prajeesh/spec2d

# NLAT - number of latitudes (should be a even number)
NLAT=94

# MAXLON - maximum number of longitudes in octahedral reduced grid, usually MAXLON=20+(NLAT/2-1)*4
MAXLON=$((20+(NLAT/2-1)*4))

# ocean_grid - MOM ocean_grid, this required to create grid_spec.nc, for a Atmosphere only run
# any arbitrary ocean_grid can be used. The land area frac, as a result the land grid points are
# defined on the basis of non-wet points of ocean_grid. 
# ocean_grid is prepared by fms preproccessing/generate_grids.

ocean_grid=$rootdir/data/ocean_grid_360x200x50.nc


p_amfigrid=$rootdir/exec/amfi_grid/amfi_grid
p_xgrid=$rootdir/exec/xgrid/xgrid
p_p2rxgrid=$rootdir/exec/p2r_xgrid/p2r_xgrid
stackNfold=$rootdir/src/preprocessing/StackNFold/StackNFold.sh

datadir=$rootdir/data

SLOPETYPE=$datadir/global_slope.1x1.nc
SOILTYPE=$datadir/global_soiltype.1x1.nc
VEGTYPE=$datadir/global_vegtype.1x1.nc
EMIS=$datadir/global_emis.1x1.nc
MTN=$datadir/global_mtn.nc
TG3=$datadir/global_tg3clim.nc
VEGFRAC=$datadir/global_vegfrac.1x1.nc
ALBEDO=$datadir/global_albedo4.1x1.nc
ZORL=$datadir/global_zorclim.1x1.nc
ZORLV=zorl
OZONE=$datadir/o3clim.nc
TOPO=$datadir/ETOPO1_topography.nc
TEMP_PRES=$datadir/temp_pres.nc

SST=$datadir/sst_monclim.nc
ICE=$datadir/ice_monclim.nc

cat <<EOF >input.nml
&fms_io_nml
  threading_read='multi'
  threading_write='single'
  fileset_write='single'
  max_files_r = 200
  max_files_w = 200
/

&fms_nml
  clock_grain='loop' ! 'component' ! 'routine' !
  domains_stack_size = 8000000
  stack_size =0
/
EOF


if [ ! -f amfi_grid.nc ]; then
	echo "----------Preparing AMFI grid ${MAXLON}x${NLAT}------------------------"
	$p_amfigrid <<< $NLAT
	if [[ $? -ne 0 ]] ; then
		rm -f amfi_grid.nc
	    exit 1
	fi
fi

if [ ! -f grid_spec.nc ]; then
	echo "----------Preparing grid_spec------------------------"
	$p_xgrid -o $ocean_grid -a amfi_grid.nc
	if [[ $? -ne 0 ]] ; then
		rm -f grid_spec.nc
	    exit 1
	fi
fi

echo "----------Preparing INPUT files------------------------"

if [ ! -f emis_ref.nc ]; then
	echo "$EMIS"
	$stackNfold -x $MAXLON -y $NLAT -i $EMIS     -o emis_ref.nc -p interpmethod=linint
	if [[ $? -ne 0 ]] ; then
		rm -f emis_ref.nc
	    exit 1
	fi
fi

if [ ! -f mtn.nc ]; then
	echo "$MTN"
	$stackNfold -x $MAXLON -y $NLAT -i $MTN -o mtn.nc -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f mtn.nc
	    exit 1
	fi
fi

if [ ! -f ozone.nc ]; then
	echo "$OZONE"
	$stackNfold -x $MAXLON -y $NLAT -i $OZONE -o ozone.nc -p interpmethod=linint
	if [[ $? -ne 0 ]] ; then
		rm -f ozone.nc
	    exit 1
	fi
fi

if [ ! -f tg3.nc ]; then
	echo "$TG3"
	$stackNfold -x $MAXLON -y $NLAT -i $TG3 -o tg3.nc -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f tg3.nc
	    exit 1
	fi
fi

if [ ! -f vegfrac.nc ]; then
	echo "$VEGFRAC"
	$stackNfold -x $MAXLON -y $NLAT -i $VEGFRAC  -o vegfrac.nc -v vegfrac  -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f vegfrac.nc
	    exit 1
	fi
fi

if [ ! -f albedo.nc ]; then
	echo "$ALBEDO"
	$stackNfold -x $MAXLON -y $NLAT -i $ALBEDO   -o albedo.nc -v alvwf,alnwf,facwf,facsf,alnsf,alvsf  -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f albedo.nc
	    exit 1
	fi
fi

if [ ! -f zorl.nc ]; then
	echo "$ZORL"
	$stackNfold -x $MAXLON -y $NLAT -i $ZORL -o zorl.nc  -v $ZORLV   -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f zorl.nc
	    exit 1
	fi
fi

if [ ! -f topography.nc ]; then
	echo "$TOPO"
	$stackNfold -x $MAXLON -y $NLAT -i $TOPO     -o topography.nc  -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f topography.nc
	    exit 1
	fi
fi

if [ ! -f temp_pres.nc ]; then
	echo "$TEMP_PRES"
	$stackNfold -x $MAXLON -y $NLAT -i $TEMP_PRES     -o temp_pres.nc     -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f temp_pres.nc
	    exit 1
	fi
fi

if [ ! -f sea_ice_forcing.nc ]; then
	echo "$ICE"
	$stackNfold -x $MAXLON -y $NLAT -i $ICE -o sea_ice_forcing.nc  -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f sea_ice_forcing.nc
	    exit 1
	fi
fi

if [ ! -f sst_forcing.nc ]; then
	echo "$SST"
	$stackNfold -x $MAXLON -y $NLAT -i $SST -o sst_forcing.nc  -p interpmethod=conserve
	if [[ $? -ne 0 ]] ; then
		rm -f sst_forcing.nc
	    exit 1
	fi
fi

if [ ! -f vegtype.nc ]; then
	echo "$VEGTYPE"
	$stackNfold -x $MAXLON -y $NLAT -i $VEGTYPE -m grid_spec.nc -o vegtype.nc -p interpmethod=dtype:conv2int=T
	if [[ $? -ne 0 ]] ; then
		rm -f vegtype.nc
	    exit 1
	fi
fi

if [ ! -f soiltype.nc ]; then
	echo "$SOILTYPE"
	$stackNfold -x $MAXLON -y $NLAT -i $SOILTYPE -m grid_spec.nc -o soiltype.nc -p interpmethod=dtype:conv2int=T
	if [[ $? -ne 0 ]] ; then
		rm -f soiltype.nc
	    exit 1
	fi
fi

if [ ! -f slopetype.nc ]; then
	echo "$SLOPETYPE"
	$stackNfold -x $MAXLON -y $NLAT -i $SLOPETYPE -m grid_spec.nc -o slopetype.nc -p interpmethod=dtype:conv2int=T
	if [[ $? -ne 0 ]] ; then
		rm -f slopetype.nc
	    exit 1
	fi
fi

cat <<EOF >atm.res
1        (Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)
1950     1     1     0     0     0        Model start time
1950     1     1     0     0     0        Model Current time
EOF


cat <<EOF >data_table

"ATM",  "facsf",    "facsf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "facwf",    "facwf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "alvsf",    "alvsf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "alnsf",    "alnsf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "alvwf",    "alvwf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "alnwf",    "alnwf",    "./INPUT/albedo.nc",       .true.,  1.0
"ATM",  "ozone",    "o3",       "./INPUT/ozone.nc",        .true.,  1.0
"ATM",  "zorl",     "zorl",     "./INPUT/zorl.nc",         .true.,  1.0
"ATM",  "vegfrac",  "vegfrac",  "./INPUT/vegfrac.nc",      .true.,  1.0
"ATM",  "sst",      "sst",      "./INPUT/sst_forcing.nc",  .true.,  1.0
"ATM",  "fice",     "sic",      "./INPUT/sea_ice_forcing.nc",  .true.,  1.0

EOF

cat <<EOF >field_table

 "TRACER", "atmos_mod", "sphum"
           "longname",     "specific humidity"
           "units",        "kg/kg"
           "profile_type", "fixed", "surface_value=3.e-6" /

 "TRACER", "atmos_mod", "clw"
           "longname",     "cloud liquid water"
           "units",        "kg/kg"
           "profile_type", "fixed", "surface_value=3.e-6" /

EOF

cat <<EOF >diag_table
AMFI
1900 01 01 00 00 00

atm_out%4yr%2mo%2dy%2hr%2mi, 1,  days, 1, days, time, 1, days,

am_phys,   rlds,      rlds,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,   rlus,      rlus,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,   sfcemis,   sfcemis,   atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,   tskin,     tskin,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  dsphumphy,   dsphumphy,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  dtphy,   dtphy,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  dtrd,     dtrd,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  duphy,   duphy,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  dvphy,   dvphy,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  fpr,     fpr,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  fprmp,     fprmp,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  hpbl,     hpbl,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  kcnv,     kcnv,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  lhflx,     lhflx,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  lpr,     lpr,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  lprcu,     lprcu,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  lprmp,     lprmp,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  pr,      pr,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  rsds,      rsds,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  rsus,      rsus,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  shflx,     shflx,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  ta,      ta,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  taux,     taux,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_phys,  tauy,     tauy,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2

am_rad, albedo,    albedo,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, clouds01,  clouds01,  atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, coszen,    coszen,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, ozone,     ozone,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rlut,      rlut,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rlutc,     rlutc,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rsdsc,     rsdsc,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rsdt,      rsdt,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rsusc,     rsusc,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rsut,      rsut,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, rsutc,     rsutc,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, ruvds,     ruvds,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_rad, ruvdsc,    ruvdsc,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2

am_sfc,  cd,        cd,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  cdq,       cdq,       atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  evap,      evap,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  hflx,      hflx,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  qss,       qss,       atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  rb,        rb,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  stress,    stress,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  uustar,    uustar,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  wind,      wind,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  canopy,     canopy,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  ec,         ec,         atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  edir,       edir,       atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  esnow,      esnow,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  et,         et,         atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  eta,        eta,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  etp,        etp,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  runoff1,    runoff1,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  sheat,      sheat,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  sheleg,     sheleg,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  slc,        slc,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  slopetype,  slopetype,  atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  smc,        smc,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  smcdry,     smcdry,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  smcmax,     smcmax,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  smcref,     smcref,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  smcwlt,     smcwlt,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  sncovr,     sncovr,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  snomlt,     snomlt,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  snwdph,     snwdph,     atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  soilm,      soilm,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  soiltype,   soiltype,   atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  soilw,      soilw,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  ssoil,      ssoil,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  stc,        stc,        atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  trans,      trans,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  tslnd,      tslnd,      atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  vegfrac,    vegfrac,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2
am_sfc,  vegtype,    vegtype,    atm_out%4yr%2mo%2dy%2hr%2mi,  all,  .true.,  none,  2

EOF

cat <<EOF >input.nml

&amfi_nml
	months=12
	days=0
	dt_atmos=600
/

&phys_nml
	dt_rad=3600
/

&atmos_nml
	trunc=62
	num_lat = 94
  	nlev = 64
	packed = .true.
	reduced = .true.
	layout = 8,32
/

&spectral_dynamics_nml
/

&cu_conv_nml
  conv_scheme='SASNEW'
/

&spherical_nml
  debug = .false.
/

&grid_fourier_nml
  plan_level = 3
  debug = .false.
/


&radiation_nml
ozone_fnm="INPUT/ozone.nc"
/

&set_aerosols_nml
use_this_module=.false.
/

&albedo_nml
/ 

&diag_manager_nml
  max_axes = 100,
  max_num_axis_sets = 100,
  max_input_fields = 699
  max_output_fields = 699
  mix_snapshot_average_fields=.false.
  issue_oor_warnings = .false.
  !do_diag_field_log = .true.
/

&fms_io_nml
  threading_read='multi'
  threading_write='single'
  fileset_write='single'
  max_files_r = 200
  max_files_w = 200
/

&fms_nml
  clock_grain='loop' ! 'component' ! 'routine' !
  domains_stack_size = 8000000
  stack_size =0
/

EOF
