#!/bin/bash
set -e

export FC=mpiifort
export CC=mpiicc
export MPICC=mpiicc
export F77=mpiifort


rootdir=$(pwd)
execdir="$rootdir/exec"
srcdir="$rootdir/src"
mkmf="$rootdir/bin/mkmf"

mkmftemplate="$rootdir/bin/mkmf.template"

numproc=16

while getopts 'gj:' flag; do
    case "${flag}" in
    g) mkmftemplate="$rootdir/bin/mkmf.template.debug" ;;
	j) numproc="$OPTARG" ;;
    *)
		echo "error"
		exit 1
        ;;
    esac
done

shift $(expr $OPTIND - 1)

echo $mkmftemplate

#-------------------------MAKE MPPNCCOMBINE--------------------------------------
cppDef="-Duse_netCDF -Duse_libMPI"  
exe=mppncc
paths="$srcdir/postprocessing/mppnccombine"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p $exe -t $mkmftemplate $paths
make -j $numproc
#--------------------------------------------------------------------------------	


# MAKE FFTW
#--------------------------------------------------------------------------------	
cppDef=""  
export LD=$FC
if [ ! -f $execdir/fftw/lib/libfftw3.a ]; then
	cd $srcdir/shared/fftw-3.3.8
	./configure --prefix=$execdir/fftw --enable-mpi --enable-openmp --enable-threads
	make clean
	make -j $numproc 
	make install
fi
#--------------------------------------------------------------------------------	


# MAKE FMS LIBRARY
#--------------------------------------------------------------------------------	
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

cppDef="-Duse_netCDF -Duse_libMPI"  
mkdir -p $execdir/lib_fms
cd $execdir/lib_fms
$mkmf -c "$cppDef" -f -p lib_fms.a -t $mkmftemplate $paths
make -j $numproc 
#--------------------------------------------------------------------------------	


#-------------------------MAKE AMFI_GRID--------------------------------------
cppDef="-Duse_netCDF -Duse_libMPI"  
exe=amfi_grid
paths="$srcdir/preprocessing/make_grids/amfi $srcdir/amfi/ocpack \
		$srcdir/amfi/transforms/gauss_and_legendre.F90"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc 
#--------------------------------------------------------------------------------	


#-------------------------MAKE XGRID--------------------------------------
cppDef="-Duse_netCDF"  
exe=xgrid
paths="$srcdir/preprocessing/make_grids/xgrid"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate $paths
make -j $numproc 
#--------------------------------------------------------------------------------	


#-------------------------MAKE XREGRID--------------------------------------
cppDef="-Duse_netCDF"  
exe=xregrid
paths="$srcdir/preprocessing/make_grids/regrid"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc 
#--------------------------------------------------------------------------------	


##-------------------------make regrid_p2r--------------------------------------
#cppDef="-Duse_netCDF"  
#exe=regrid_p2r
#paths="$srcdir/postprocessing/regrid_p2r"
#export LD=$FC
#mkdir -p $execdir/$exe
#cd $execdir/$exe
#OPTS="-I$execdir/lib_fms"
#LIBS="$execdir/lib_fms/lib_fms.a"
#$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
#make 
##--------------------------------------------------------------------------------	


#-------------------------MAKE P2R_XGRID--------------------------------------
cppDef="-Duse_netCDF -Dlib_xgrid"  
exe=p2r_xgrid
paths="$srcdir/preprocessing/make_grids/xgrid $srcdir/preprocessing/make_grids/p2r_xgrid \
		$srcdir/amfi/ocpack"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc
#--------------------------------------------------------------------------------	


# MAKE AMFI
#--------------------------------------------------------------------------------	
amfi="$srcdir/amfi"
#paths="$amfi/model $amfi/driver $amfi/radiation $amfi/spec_dyn"
paths=$(find $amfi -type d)

#cppDef="-Duse_netCDF -Duse_libMPI -DMPI3"  
cppDef="-Duse_netCDF -Duse_libMPI"  
exe=spec2d
mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms -I$execdir/fftw/include"

LIBS="$execdir/lib_fms/lib_fms.a $execdir/fftw/lib/libfftw3_mpi.a $execdir/fftw/lib/libfftw3.a"

$mkmf -c "$cppDef" -f -p ${exe}.exe -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths

make -j $numproc
#--------------------------------------------------------------------------------	

cd $rootdir/work/ocy

rm -f fftw.*

mpirun -np 8 -prepend-rank $rootdir/exec/spec2d/spec2d.exe

rm -f atm_out.nc

./mppncc -r atm_out.nc

