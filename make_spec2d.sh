#!/bin/bash
set -e

export FC=mpiifort
export CC=mpiicc
export MPICC=mpiicc
export F77=mpiifort

thisdir=$(pwd)

mkmftemplate="$thisdir/bin/mkmf.template"
numproc=16

while getopts 'gj:' flag; do
    case "${flag}" in
    g) mkmftemplate="$thisdir/bin/mkmf.template.debug" ;;
	j) numproc="$OPTARG" ;;
    *)
		echo "error"
		exit 1
        ;;
    esac
done

shift $(expr $OPTIND - 1)

echo $mkmftemplate

EXE="spec2d"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

amfi="$thisdir/amfi"

#------------AMFI SRC------------------------------------------
#paths="$amfi/model $amfi/driver $amfi/radiation $amfi/spec_dyn"
paths=$(find $amfi -type d)
#--------------------------------------------------------------------------------	

#------------------------libFMS SRC----------------------------------------------	
libfmspaths="$thisdir/shared/mpp $thisdir/shared/include \
       $thisdir/shared/mpp/include \
       $thisdir/shared/fms $thisdir/shared/platform \
       $thisdir/shared/memutils $thisdir/shared/constants \
       $thisdir/shared/horiz_interp $thisdir/shared/mosaic \
	   $thisdir/shared/time_manager $thisdir/shared/data_override \
       $thisdir/shared/time_interp $thisdir/shared/axis_utils \
       $thisdir/shared/astronomy $thisdir/shared/diag_manager \
       $thisdir/shared/sat_vapor_pres $thisdir/shared/mersenne_twister \
	   $thisdir/shared/tracer_manager $thisdir/shared/field_manager \
	   $thisdir/shared/strman"
#--------------------------------------------------------------------------------	



#-------------------------make mppnccombine--------------------------------------
cppDef="-Duse_netCDF -Duse_libMPI"  
exe=mppncc
tpath="$thisdir/mppnccombine"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p $exe -t $mkmftemplate $tpath
make 
#--------------------------------------------------------------------------------	

# make FFTW
#--------------------------------------------------------------------------------	
cppDef="-Duse_netCDF -Duse_libMPI"  
export LD=$FC
if [ ! -f $execdir/fftw/lib/libfftw3.a ]; then
	cd $thisdir/shared/fftw-3.3.8
	./configure --prefix=$execdir/fftw --enable-mpi --enable-openmp --enable-threads
	make clean
	make -j 16
	make install
fi
#--------------------------------------------------------------------------------	


#make FMS Library
#--------------------------------------------------------------------------------	
cppDef="-Duse_netCDF -Duse_libMPI"  
mkdir -p $execdir/lib_fms
cd $execdir/lib_fms
$mkmf -c "$cppDef" -f -p lib_fms.a -t $mkmftemplate $libfmspaths
make -j 16
#--------------------------------------------------------------------------------	


#-------------------------make amfi_grid--------------------------------------
cppDef="-Duse_netCDF -Duse_libMPI"  
exe=amfi_grid
tpath="$thisdir/make_grids/amfi $thisdir/amfi/ocpack \
		$thisdir/amfi/transforms/gauss_and_legendre.F90"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $tpath
make 
#--------------------------------------------------------------------------------	

#-------------------------make xgrid--------------------------------------
cppDef="-Duse_netCDF"  
exe=xgrid
tpath="$thisdir/make_grids/xgrid"
export LD=$CC
mkdir -p $execdir/$exe
cd $execdir/$exe
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate $tpath
make 
#--------------------------------------------------------------------------------	

#-------------------------make xregrid--------------------------------------
cppDef="-Duse_netCDF"  
exe=xregrid
tpath="$thisdir/make_grids/regrid"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $tpath
make 
#--------------------------------------------------------------------------------	

#-------------------------make xregrid--------------------------------------
cppDef="-Duse_netCDF"  
exe=xregrid
tpath="$thisdir/make_grids/regrid"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $tpath
make 
#--------------------------------------------------------------------------------	

#-------------------------make p2r_xgrid--------------------------------------
cppDef="-Duse_netCDF -Dlib_xgrid"  
exe=p2r_xgrid
tpath="$thisdir/make_grids/xgrid $thisdir/make_grids/p2r_xgrid $thisdir/amfi/ocpack"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $tpath
make 
#--------------------------------------------------------------------------------	


#make AMFI
#--------------------------------------------------------------------------------	
cppDef="-Duse_netCDF -Duse_libMPI -DMPI3"  
mkdir -p $execdir/$EXE
cd $execdir/$EXE

OPTS="-I$execdir/lib_fms -I$execdir/fftw/include"

LIBS="$execdir/lib_fms/lib_fms.a $execdir/fftw/lib/libfftw3_mpi.a $execdir/fftw/lib/libfftw3.a"

$mkmf -c "$cppDef" -f -p ${EXE}.exe -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths

make -j $numproc
#--------------------------------------------------------------------------------	

#cd $thisdir/work_ocy

#rm -f fftw.*

#mpirun -np 8 -prepend-rank $thisdir/exec/spec2d/spec2d.exe

#rm -f atm_out.nc

#./mppncc -r atm_out.nc

#cd $thisdir

#$execdir/amfi_grid/amfi_grid  <<< 94
