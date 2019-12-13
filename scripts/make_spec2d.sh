#!/bin/bash --login
set -e

AQUAPLANET=False


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#----------------------------No editing needed beyond this-----------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------



# Path to the root directory.
rootdir=_ROOTDIR_

machine=$(cat $rootdir/bin/._machine_)
execdir="$rootdir/exec"
srcdir="$rootdir/src"
mkmf="$rootdir/bin/mkmf"

mkmftemp="$rootdir/bin/mkmf.template"
envir="$rootdir/bin/env.$machine"

source $envir

numproc=16

debug=False

usage() {echo "Usage: $0 [-g] [-j numproc]" 1>&2; exit 1;}


while getopts 'gj:' flag; do
    case "${flag}" in
	g) debug=True ;;
	j) numproc="$OPTARG" ;;
    *)
		echo "error"
		usage
        ;;
    esac
done

if [[ "$debug" == True ]]; then
	mkmftemplate=${mkmftemp}_debug.$machine
else
	mkmftemplate=${mkmftemp}.$machine
fi

shift $(expr $OPTIND - 1)

cat $mkmftemplate

echo "#--------------------------MAKE MPPNCCOMBINEP2R-----------------------------------"
cppDef=""
exe=mppnccp2r
paths="$srcdir/postprocessing/mppnccombinep2r/mppnccombinep2r.c"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p $exe -t $mkmftemplate $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


# MAKE FFTW

if [ -z "$FFTW_DIR" ]; then
	echo "#-------------------------------MAKE FFTW--------------------------------------"
	cppDef=""
	export LD=$FC
	if [ ! -f $execdir/fftw/lib/libfftw3.a ]; then
		cd $srcdir/shared/fftw-3.3.8
		./configure --prefix=$execdir/fftw --enable-mpi --enable-openmp --enable-threads --disable-doc
		make clean
		make -j $numproc
		make install
		make clean
	fi
	export FFTW_DIR=$execdir/fftw
	echo "#------------------------------------------------------------------------------"
fi


echo "#-----------------------------MAKE FMS LIBRARY NO-MPI---------------------------------"
paths="$srcdir/shared/mpp $srcdir/shared/include \
       $srcdir/shared/mpp/include \
       $srcdir/shared/fms $srcdir/shared/platform \
       $srcdir/shared/memutils $srcdir/shared/constants \
       $srcdir/shared/horiz_interp $srcdir/shared/mosaic \
	   $srcdir/shared/time_manager $srcdir/shared/data_override \
       $srcdir/shared/time_interp $srcdir/shared/axis_utils \
       $srcdir/shared/astronomy $srcdir/shared/diag_manager \
       $srcdir/shared/sat_vapor_pres $srcdir/shared/mersenne_twister \
	   $srcdir/shared/tracer_manager $srcdir/shared/field_manager \
	   $srcdir/shared/strman"

cppDef="-Duse_netCDF"
mkdir -p $execdir/lib_fms_nompi
cd $execdir/lib_fms_nompi
$mkmf -c "$cppDef" -f -p lib_fms_nompi.a -t $mkmftemplate $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


echo "#-----------------------------MAKE FMS LIBRARY MPI---------------------------------"
paths="$srcdir/shared/mpp $srcdir/shared/include \
       $srcdir/shared/mpp/include \
       $srcdir/shared/fms $srcdir/shared/platform \
       $srcdir/shared/memutils $srcdir/shared/constants \
       $srcdir/shared/horiz_interp $srcdir/shared/mosaic \
	   $srcdir/shared/time_manager $srcdir/shared/data_override \
       $srcdir/shared/time_interp $srcdir/shared/axis_utils \
       $srcdir/shared/astronomy $srcdir/shared/diag_manager \
       $srcdir/shared/sat_vapor_pres $srcdir/shared/mersenne_twister \
	   $srcdir/shared/tracer_manager $srcdir/shared/field_manager \
	   $srcdir/shared/strman"

cppDef="-Duse_netCDF -Duse_libMPI -DOVERLOAD_C8"
mkdir -p $execdir/lib_fms
cd $execdir/lib_fms
$mkmf -c "$cppDef" -f -p lib_fms.a -t $mkmftemplate $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"

#echo "#--------------------------MAKE topo_regularization-----------------------------------"
#cppDef="-Duse_netCDF -Duse_libMPI" # -Dtest_gauss_legendre"
#exe=topo_regularization
#paths="$srcdir/preprocessing/topo_regularization/ $srcdir/shared/fft"
#export LD=$FC
#mkdir -p $execdir/$exe
#cd $execdir/$exe
#OPTS="-I$execdir/lib_fms"
#LIBS="$execdir/lib_fms/lib_fms.a"
#
#$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
#make -j $numproc
#echo "#------------------------------------------------------------------------------"



echo "#--------------------------MAKE AMFI_GRID-----------------------------------"
cppDef="-Duse_netCDF "
exe=amfi_grid
paths="$srcdir/preprocessing/make_grids/amfi $srcdir/amfi/ocpack \
		$srcdir/amfi/transforms/gauss_and_legendre.F90"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms_nompi"
LIBS="$execdir/lib_fms_nompi/lib_fms_nompi.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


echo "#--------------------------MAKE XGRID-----------------------------------"
cppDef="-Duse_netCDF"
exe=xgrid
paths="$srcdir/preprocessing/make_grids/xgrid"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


#echo "#--------------------------MAKE XREGRID-----------------------------------"
#cppDef="-Duse_netCDF"
#exe=xregrid
#paths="$srcdir/preprocessing/make_grids/regrid"
#export LD=$FC
#mkdir -p $execdir/$exe
#cd $execdir/$exe
#OPTS="-I$execdir/lib_fms"
#LIBS="$execdir/lib_fms/lib_fms.a"
#$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
#make -j $numproc
#echo "#------------------------------------------------------------------------------"


echo "#--------------------------MAKE P2R_XGRID-----------------------------------"
cppDef="-Duse_netCDF -Dlib_xgrid"
exe=p2r_xgrid
paths="$srcdir/preprocessing/make_grids/xgrid \
       $srcdir/preprocessing/make_grids/p2r_xgrid \
		$srcdir/amfi/ocpack $srcdir/amfi/transforms/gauss_and_legendre.F90"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms_nompi"
LIBS="$execdir/lib_fms_nompi/lib_fms_nompi.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


echo "#--------------------------MAKE spectral_topo-----------------------------------"
cppDef="-Duse_netCDF -Duse_libMPI" # -Dtest_gauss_legendre"
exe=spectral_topo
paths="$srcdir/preprocessing/spectral_topo/ \
		$srcdir/amfi/ocpack $srcdir/amfi/transforms/"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"

if [ ! -z "$FFTW_DIR" ]; then
	OPTS="$OPTS -I$FFTW_DIR/include"
	LIBS="$LIBS $FFTW_DIR/lib/libfftw3_mpi.a $FFTW_DIR/lib/libfftw3.a"
fi

$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"


echo "#--------------------------listlayout-----------------------------------"
cppDef=""
exe=listlayout
paths="$srcdir/preprocessing/listlayout \
	   $srcdir/amfi/ocpack" 
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate $paths
make -j $numproc
echo "#------------------------------------------------------------------------------"

echo "#-------------------------MAKE RUN_NCCOMBINEP2R--------------------------------------"
cppDef="-Dlib_mppnccp2r -Duse_libMPI"
exe=run_mppnccp2r
paths="$srcdir/postprocessing/mppnccombinep2r"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms"

LIBS="$execdir/lib_fms/lib_fms.a"

$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths
make -j $numproc
echo "#--------------------------------------------------------------------------------"


echo "#----------------------------MAKE AMFI-----------------------------------------"
amfi="$srcdir/amfi"
#paths="$amfi/model $amfi/driver $amfi/radiation $amfi/spec_dyn"
paths=$(find $amfi -type d)

cppDef="-Duse_netCDF -Duse_libMPI"

if [[ "$MPI3" == True ]]; then
	cppDef=$cppDef" -DMPI3"
fi

exe=spec2d
if [[ "$AQUAPLANET" == True ]]; then
	cppDef=$cppDef" -DAQUAPLANET"
	exe=$exe"_AQUAPLANET"
fi

mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"

if [ ! -z "$FFTW_DIR" ]; then
	OPTS="$OPTS -I$FFTW_DIR/include"
	LIBS="$LIBS $FFTW_DIR/lib/libfftw3_mpi.a $FFTW_DIR/lib/libfftw3.a"
fi

$mkmf -c "$cppDef" -f -p ${exe}.exe -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths

make -j $numproc
echo "#---------------------COMPILATION COMPLETED------------------------------------"


