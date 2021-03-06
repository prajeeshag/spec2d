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

cppDef="-Duse_netCDF -Duse_libMPI"  

EXE="cfs_aero"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

amfi="$thisdir/cfs_aero"

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

cppDef="-Duse_netCDF -Duse_libMPI -DMPI3"  
#make AMFI
#--------------------------------------------------------------------------------	
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
