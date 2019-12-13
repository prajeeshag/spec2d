#!/bin/bash --login
set -e

# rootdir -  Absolute path to the root directory
rootdir=_ROOTDIR_

# NLAT - number of latitudes (should be a even number)
NLAT=94

TRUNC=62

usage() { echo "Usage: $0 [-l nlat] [-t tuncation]" 1>&2; exit 1;}


while getopts 'l:t:' flag; do
    case "${flag}" in
    l) NLAT="$OPTARG" ;;
    t) TRUNC="$OPTARG" ;;
    *)
        echo "error"
        usage
        ;;
    esac
done


# MAXLON - maximum number of longitudes in octahedral reduced grid, usually MAXLON=20+(NLAT/2-1)*4
MAXLON=$((20+(NLAT/2-1)*4))

# AQUAPLANET - Whether an aquaplanet experiment?
AQUAPLANET=False

# ocean_grid - MOM ocean_grid, this required to create grid_spec.nc, for a Atmosphere only run
# any arbitrary ocean_grid can be used. The land area frac, as a result the land grid points are
# defined on the basis of non-wet points of ocean_grid. 
# ocean_grid is prepared by fms preproccessing/generate_grids.

ocean_grid=$rootdir/data/ocean_grid_360x200x50.nc


#------------------------------------No need to edit beyond this---------------------------------

machine=$(cat $rootdir/bin/._machine_)
source $rootdir/bin/env.$machine

scriptdir=$rootdir/scripts

p_amfigrid=$rootdir/exec/amfi_grid/amfi_grid
p_xgrid=$rootdir/exec/xgrid/xgrid
p_p2rxgrid=$rootdir/exec/p2r_xgrid/p2r_xgrid
p_listlayout=$rootdir/exec/listlayout/listlayout

stackNfold=$rootdir/src/preprocessing/StackNFold/StackNFold.sh

datadir=$rootdir/data


inDataOpt=()
inDataOpt+=("$datadir/global_slope.1x1.nc    slopetype.nc 		T dtype    FF T")
inDataOpt+=("$datadir/global_soiltype.1x1.nc soiltype.nc  		T dtype    FF T")
inDataOpt+=("$datadir/global_vegtype.1x1.nc  vegtype.nc   		T dtype    FF T")
inDataOpt+=("$datadir/global_emis.1x1.nc     emis_ref.nc  		F linint   FF F")
inDataOpt+=("$datadir/global_mtn.nc          mtn.nc       		F conserve FF F")
inDataOpt+=("$datadir/global_tg3clim.nc      tg3.nc       		F conserve FF F")
inDataOpt+=("$datadir/o3clim.nc              ozone.nc     		F linint   FF F")
inDataOpt+=("$datadir/global_vegfrac.1x1.nc  vegfrac.nc   		F conserve FF F vegfrac")
inDataOpt+=("$datadir/global_albedo4.1x1.nc  albedo.nc    		F conserve FF F alvwf,alnwf,facwf,facsf,alnsf,alvsf")
inDataOpt+=("$datadir/global_zorclim.1x1.nc  zorl.nc      		F conserve FF F zorl")
inDataOpt+=("$datadir/ETOPO1_topography.nc   topography.nc		F conserve FF F")
inDataOpt+=("$datadir/temp_pres.nc   		 temp_pres.nc		F conserve FF F")
inDataOpt+=("$datadir/sst_monclim.nc   		 sst_forcing.nc		F conserve FF F")
inDataOpt+=("$datadir/ice_monclim.nc   		 sea_ice_forcing.nc	F conserve FF F")


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

if [ ! -f p_xgrd.nc ]; then
	echo "----------Preparing p_xgrd ${MAXLON}x${NLAT}------------------------"
	$p_p2rxgrid <<< $NLAT
	if [[ $? -ne 0 ]] ; then
		rm -f p_xgrd.nc
	    exit 1
	fi
fi

if [ ! "$AQUAPLANET" == "True" ]; then

	if [ ! -f grid_spec.nc ]; then
		echo "----------Preparing grid_spec------------------------"
		$p_xgrid -o $ocean_grid -a amfi_grid.nc
		if [[ $? -ne 0 ]] ; then
			rm -f grid_spec.nc
			rm -f *X*.nc
		    exit 1
		fi
	fi
	
	
	echo "----------Preparing INPUT files------------------------"
	
	for ((i = 0; i < ${#inDataOpt[@]}; i++)); do
	    #inDataOpt+=("$datadir/global_slope.1x1.nc slopetype.nc T dtype FF T")
		opt="${inDataOpt[$i]}"
		opta=($opt)
		ifnm=${opta[0]}
		ofnm=${opta[1]}
		if [ "${opta[2]}" == "T" ]; then
			maskopt="-m grid_spec.nc"
		else
			maskopt=""
		fi
		intp=${opta[3]}
		ongrid=${opta[4]}
		toint=${opta[5]}
		if [ ! -z "${opta[6]}" ]; then
			varsopt="-v ${opta[6]}"
		else
			varsopt=""
		fi
	
		if [ ! -f $ofnm ]; then
			if [ -f $datadir/fix/t${TRUNC}_${MAXLON}_${NLAT}_$ofnm ]; then
				echo "linking $ofnm from $datadir/fix/t${TRUNC}_${MAXLON}_${NLAT}_$ofnm"
			    ln -sf $datadir/fix/t${TRUNC}_${MAXLON}_${NLAT}_$ofnm $ofnm
				if [[ $? -ne 0 ]] ; then
					rm -f $ofnm 
				    exit 1
				fi
			else
				echo "$opt"
				$stackNfold -x $MAXLON -y $NLAT -i $ifnm $maskopt -o $ofnm $varsopt -p interpmethod=$intp:conv2int=$toint:ongrid=$ongrid
				if [[ $? -ne 0 ]] ; then
					rm -f $ofnm
				    exit 1
				fi
			fi
		fi
	done

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

EXENAME=spec2d
if [ "$AQUAPLANET" == "True" ]; then
	EXENAME=spec2d_AQUAPLANET
fi

EXE=$rootdir/exec/$EXENAME/${EXENAME}.exe

#Making Diag Table
$scriptdir/dflo_to_dtable.sh $scriptdir/diag_field_log.out

cp $scriptdir/input.nml .
TRUNC=$(((NLAT-1)*2/3))
if [ "$(($TRUNC%2))" -ne "0" ]; then
	TRUNC=$((TRUNC-1))
fi

echo " "
echo "valid spherical truncations for NLAT = "$NLAT" is TRUNC = "$TRUNC 
echo " "

sed -i "s|_TRUNC_|$TRUNC|g" input.nml
sed -i "s|_NLAT_|$NLAT|g" input.nml  

mkdir -p INPUT
mkdir -p RESTART

mv *.nc INPUT/
mv atm.res INPUT/

cp $scriptdir/run_amfi_${JOBSCDLR}.sh .

sed -i "s|_EXE_|$EXE|g" run_amfi_${JOBSCDLR}.sh

$p_listlayout <<< "$NLAT $TRUNC" > valid_pe_layouts_${NLAT}_$TRUNC

echo " "
echo " "
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "You need to edit some informations in the run script run_amfi_lsf.sh"
echo "As well as, should set namelist values for 'trunc', 'num_lat', 'layout'"
echo "in input.nml accordingly."
echo "All possible valid pe layouts are given in the file valid_pe_layouts_${NLAT}_T$TRUNC "
echo "--------------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------------"
