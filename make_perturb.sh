#!/bin/bash
set -e

export FC=mpiifort
export CC=mpiicc
export MPICC=mpiicc
export F77=mpiifort

cppDef="-Duse_netCDF -Duse_libMPI -DOVERLOAD_C8 -Dgloopa" #-Dtest_grid_to_fourier #-Dcheck_mpi" # -Dtest_interp"

thisdir=$(pwd)

EXE="perturb_ini"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

#mkmftemplate="$thisdir/bin/mkmf.template.debug"
mkmftemplate="$thisdir/bin/mkmf.template"

paths="$thisdir/perturb_ini $thisdir/amfi/spec_dyn"

libsigiopaths=$thisdir/libsigio

#------------------------libFMS SRC----------------------------------------------	
libfmspaths="$thisdir/shared/mpp $thisdir/shared/include \
       $thisdir/shared/mpp/include \
       $thisdir/shared/fms $thisdir/shared/platform \
       $thisdir/shared/memutils $thisdir/shared/constants \
       $thisdir/shared/horiz_interp $thisdir/shared/mosaic \
	   $thisdir/shared/time_manager $thisdir/shared/data_override \
       $thisdir/shared/time_interp $thisdir/shared/axis_utils \
       $thisdir/shared/astronomy $thisdir/shared/diag_manager \
       $thisdir/shared/sat_vapor_pres"
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


#make 
#--------------------------------------------------------------------------------	
mkdir -p $execdir/lib_sigio
cd $execdir/lib_sigio

$mkmf -c "$cppDef" -f -p lib_sigio.a -t ${mkmftemplate}_wor8 $libsigiopaths

make $@
#--------------------------------------------------------------------------------	

#make 
#--------------------------------------------------------------------------------	
mkdir -p $execdir/$EXE
cd $execdir/$EXE

OPTS="-I$execdir/lib_fms -I$execdir/fftw/include -I$execdir/lib_sigio"

LIBS="$execdir/lib_fms/lib_fms.a $execdir/fftw/lib/libfftw3_mpi.a $execdir/fftw/lib/libfftw3.a $execdir/lib_sigio/lib_sigio.a -L/iitm1/cccr/sandeep/preind-aerosol/library -lw3emc_v2.0.5_d -lw3nco_v2.0.5_d -lbacio_d"

$mkmf -c "$cppDef" -f -p ${EXE}.exe -t ${mkmftemplate} -o "$OPTS" -l "$LIBS"  $paths

make $@
#--------------------------------------------------------------------------------	

cd $thisdir

exec/perturb_ini/perturb_ini.exe sig_ini
