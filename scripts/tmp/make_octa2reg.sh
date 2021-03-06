#!/bin/bash
set -e

export FC=mpiifort
export CC=mpiicc
export MPICC=mpiicc
export F77=mpiifort

cppDef="-Duse_netCDF -Duse_libMPI " #-Dtest_grid_to_fourier " 

thisdir=$(pwd)

EXE="octa2reg"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

mkmftemplate="$thisdir/bin/mkmf.template"
#mkmftemplate="$thisdir/bin/mkmf.template.debug"

#------------AMFI SRC------------------------------------------
paths="$thisdir/octa2reg $thisdir/amfi/spec_dyn/transforms/gauss_and_legendre.F90"
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
	   $thisdir/shared/tracer_manager $thisdir/shared/field_manager"
#--------------------------------------------------------------------------------	

#-------------------------mppnccombine SRC---------------------------------------	
mppnccpath="$thisdir/mppnccombine"
#-------------------------------------------------------------------------------



#-------------------------make mppnccombine--------------------------------------
export LD=mpiicc
mkdir -p $execdir/mppncc
cd $execdir/mppncc
$mkmf -c "$cppDef" -f -p mppncc -t $mkmftemplate $mppnccpath
make 
#--------------------------------------------------------------------------------	



export LD=mpiifort
# make FFTW
#--------------------------------------------------------------------------------	
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
mkdir -p $execdir/lib_fms
cd $execdir/lib_fms
$mkmf -c "$cppDef" -f -p lib_fms.a -t $mkmftemplate $libfmspaths
make -j 16
#--------------------------------------------------------------------------------	



#make AMFI
#--------------------------------------------------------------------------------	
mkdir -p $execdir/$EXE
cd $execdir/$EXE

OPTS="-I$execdir/lib_fms -I$execdir/fftw/include"

LIBS="$execdir/lib_fms/lib_fms.a $execdir/fftw/lib/libfftw3_mpi.a $execdir/fftw/lib/libfftw3.a"

$mkmf -c "$cppDef" -f -p ${EXE}.exe -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths

make $@
#--------------------------------------------------------------------------------	

cd $thisdir

$thisdir/exec/$EXE/${EXE}.exe

