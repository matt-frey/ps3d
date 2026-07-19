#!/bin/bash

export NJOBS=4

if [[ "$PREFIX" == "" ]]; then
    echo "Missing 'PREFIX' environment variable. Please set this variable."
    exit
fi

if [[ "$DOWNLOADS_DIR" == "" ]]; then
    echo "Missing 'DOWNLOADS_DIR' environment variable. Please set this variable."
    exit
fi

if [[ "$SRC_DIR" == "" ]]; then
    echo "Missing 'SRC_DIR' environment variable. Please set this variable."
    exit
fi

mkdir -p $DOWNLOADS_DIR
mkdir -p $SRC_DIR

# -----------------------------------------------------------------------------
gcc_compiler=$(which gcc)

if [[ ! "$gcc_compiler" ]]; then
    echo "Missing GNU gcc compiler"
    exit
fi

echo "Found gcc: $gcc_compiler"

gfortran_compiler=$(which gfortran)

if [[ ! "$gfortran_compiler" ]]; then
    echo "Missing GNU gfortran compiler"
    exit
fi

echo "Found gfortran: $gfortran_compiler"

gpp_compiler=$(which g++)

if [[ ! "$gpp_compiler" ]]; then
    echo "Missing GNU g++ compiler"
    exit
fi

echo "Found g++ $gpp_compiler"

export FC=$gfortran_compiler
export CC=$gcc_compiler
export CXX=$gpp_compiler

mpif90_compiler=$(which mpif90)

if [[ ! "$mpif90_compiler" ]]; then
    echo "MPI not found. Installing MPI now."
    bash openmpi.sh

    MPI_PREFIX=$PREFIX/openmpi


    if [ ! -f "$MPI_PREFIX/bin/mpif90" ]; then
        echo "MPI not properly installed"
        exit
    fi
else
    echo "Found MPI compiler: $mpif90_compiler"
    MPI_PREFIX=${mpif90_compiler%/*/*}
fi

export MPI_DIR=$MPI_PREFIX
export MPIF90=$MPI_DIR/bin/mpif90
export MPICC=$MPI_DIR/bin/mpicc
export MPICXX=$MPI_DIR/bin/mpicxx

export LD_LIBRARY_PATH=${MPI_DIR}/lib:$LD_LIBRARY_PATH
export LDFLAGS="$LDFLAGS -L${MPI_DIR}/lib"
export CPPFLAGS="$CPPFLAGS -I${MPI_DIR}/include"

# -----------------------------------------------------------------------------
h5cpp_compiler=$(which h5pcc) # h5c++)

is_h5_parallel=''
if [[ "$h5cpp_compiler" ]]; then
    is_h5_parallel=$($h5cpp_compiler -showconfig | grep "Parallel HDF5: yes")
fi

if [[ ! "$is_h5_parallel" ]]; then
    HDF5_PREFIX=$PREFIX/hdf5

    bash hdf5.sh

    if [ ! -f "$HDF5_PREFIX/bin/h5dump" ]; then
        echo "HDF5 not properly installed"
        exit
    fi
else
    export HDF5_PREFIX=${h5cpp_compiler%/*/*}
    echo "Found HDF5: $HDF5_PREFIX"
fi

# we need to append these paths
export LD_LIBRARY_PATH=${HDF5_PREFIX}/lib:$LD_LIBRARY_PATH
export LDFLAGS="$LDFLAGS -L${HDF5_PREFIX}/lib"
export CPPFLAGS="$CPPFLAGS -I${HDF5_PREFIX}/include"

# -----------------------------------------------------------------------------
nc_config=$(which nc-config)

is_nc_parallel='no'

if [[ "$nc_config" ]]; then
    is_nc_parallel=$($nc_config --has-parallel)
fi

if [[ "$is_nc_parallel" == "no" ]]; then

    NETCDF_PREFIX=$PREFIX/netcdf

    bash netcdf.sh

    if [ ! -f "$NETCDF_PREFIX/bin/nc-config" ]; then
        echo "netCDF-C not properly installed"
        exit
    fi
else
    NETCDF_PREFIX=${nc_config%/*/*}
    echo "Found netCDF-C $NETCDF_PREFIX"
fi

# -----------------------------------------------------------------------------

if [ ! -f "$NETCDF_PREFIX/lib/libnetcdff.a" ]; then
    echo "Could not find netCDF-Fortran."
fi

export NETCDF_C_DIR=$NETCDF_PREFIX
export NETCDF_FORTRAN_DIR=$NETCDF_PREFIX


#nf_config=$(which nf-config)

#if [[ "$nf_config" ]]; then
#    NETCDFF_PREFIX=${nf_config%/*/*}
#    echo "Found netCDF-Fortran $NETCDFF_PREFIX"
#    export NETCDF_FORTRAN_DIR=$NETCDFF_PREFIX
#else
#    echo "Did not find netcdf-Fortran"
#fi

# -----------------------------------------------------------------------------


#export EPIC_PREFIX=$PREFIX/epic

bash ps3d.sh
