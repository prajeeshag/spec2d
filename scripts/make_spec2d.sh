#!/bin/bash
set -e

# Path to the root directory.
rootdir=/moes/home/prajeesh/spec2d

# Fortran compiler
export FC=mpiifort
export F77=mpiifort

# C compiler
export CC=mpiicc
export MPICC=mpiicc

# if mpi library version is 3 or above
MPI3=True

#netcdf library path
NETCDF=/gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.2.1

#Fortran compiler options
FFLAGS="-r8 -O2 -fp-model precise -convert big_endian -align array32byte -I$NETCDF/include"
#Fortran compiler debug options
DFFLAGS="-g -traceback -fpe0 -fp-stack-check -check all -check noarg_temp_created"

# C compiler options
CFLAGS="-O2 -I$NETCDF/include"
# C compiler debug options
DCFLAGS="-g -traceback"

#Linker options
LDFLAGS="-L$NETCDF/lib -lnetcdf -lnetcdff -mkl -lrt -lstdc++ -lm"





#--------------------------------------------------------------------------------	
#--------------------------------------------------------------------------------	
#----------------------------No editing needed beyond this-----------------------
#--------------------------------------------------------------------------------	
#--------------------------------------------------------------------------------	

execdir="$rootdir/exec"
srcdir="$rootdir/src"
mkmf="$rootdir/bin/mkmf"

mkmftemplate="$rootdir/bin/mkmf.template"

numproc=16

debug=False
while getopts 'gj:' flag; do
    case "${flag}" in
	g) debug=True ;;
	j) numproc="$OPTARG" ;;
    *)
		echo "error"
		exit 1
        ;;
    esac
done

if [[ "$debug" == True ]]; then
cat<<EOF > $mkmftemplate
FFLAGS = $DFFLAGS $FFLAGS
CFLAGS = $DCFLAGS $CFLAGS
LDFLAGS = $LDFLAGS
EOF
else
cat<<EOF > $mkmftemplate
FFLAGS = $FFLAGS
CFLAGS = $CFLAGS
LDFLAGS = $LDFLAGS
EOF
fi

shift $(expr $OPTIND - 1)

cat $mkmftemplate

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

#-------------------------MAKE MPPNCCOMBINEP2R--------------------------------------
cppDef=""  
exe=mppnccp2r
paths="$srcdir/postprocessing/mppnccombinep2r/mppnccombinep2r.c"
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
	make clean
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


#-------------------------MAKE P2R_XGRID--------------------------------------
cppDef="-Duse_netCDF -Dlib_xgrid"  
exe=p2r_xgrid
paths="$srcdir/preprocessing/make_grids/xgrid \
       $srcdir/preprocessing/make_grids/p2r_xgrid \
		$srcdir/amfi/ocpack $srcdir/amfi/transforms/gauss_and_legendre.F90"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe
OPTS="-I$execdir/lib_fms"
LIBS="$execdir/lib_fms/lib_fms.a"
$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS" $paths
make -j $numproc
#--------------------------------------------------------------------------------	

#-------------------------MAKE RUN_NCCOMBINEP2R--------------------------------------
echo "#-------------------------MAKE RUN_NCCOMBINEP2R--------------------------------------"
cppDef="-Dlib_mppnccp2r"  
exe=run_mppnccp2r
paths="$srcdir/postprocessing/mppnccombinep2r"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms" 

LIBS="$execdir/lib_fms/lib_fms.a"

$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths
make -j $numproc
#--------------------------------------------------------------------------------	


# MAKE AMFI
#--------------------------------------------------------------------------------	
amfi="$srcdir/amfi"
#paths="$amfi/model $amfi/driver $amfi/radiation $amfi/spec_dyn"
paths=$(find $amfi -type d)

if [[ "$MPI3" == True ]]; then
	cppDef="-Duse_netCDF -Duse_libMPI -DMPI3"  
else
	cppDef="-Duse_netCDF -Duse_libMPI"
fi
exe=spec2d
mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms -I$execdir/fftw/include"

LIBS="$execdir/lib_fms/lib_fms.a $execdir/fftw/lib/libfftw3_mpi.a $execdir/fftw/lib/libfftw3.a"

$mkmf -c "$cppDef" -f -p ${exe}.exe -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths

make -j $numproc
#--------------------------------------------------------------------------------	



