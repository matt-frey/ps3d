#!/bin/bash

# recipe for:
NETCDF_C=netcdf-c
NETCDF_C_VERSION=4.10.0 #4.8.1
NETCDF_F=netcdf-fortran
NETCDF_F_VERSION=4.6.3 #4.5.4


NETCDF_PREFIX=$PREFIX/netcdf
mkdir -p "$NETCDF_PREFIX"

# download NetCDF-C
if [ ! -f "${DOWNLOADS_DIR}/$NETCDF_C-$NETCDF_C_VERSION.tar.gz" ]; then
    curl -L \
        --output "${DOWNLOADS_DIR}/$NETCDF_C-$NETCDF_C_VERSION.tar.gz" \
        "https://downloads.unidata.ucar.edu/$NETCDF_C/$NETCDF_C_VERSION/$NETCDF_C-$NETCDF_C_VERSION.tar.gz"
fi

# unpack
mkdir -p "${SRC_DIR}/$NETCDF_C" && cd "$_"
tar xvf "${DOWNLOADS_DIR}/$NETCDF_C-$NETCDF_C_VERSION.tar.gz"
#

# configure
mkdir -p "${SRC_DIR}/$NETCDF_C/build" && cd "$_"
FC=$MPIF90 CC=$MPICC CXX=$MPICXX ${SRC_DIR}/$NETCDF_C/$NETCDF_C-$NETCDF_C_VERSION/configure      \
    --enable-netcdf-4                                           \
    --with-pic                                                  \
    --prefix=${NETCDF_PREFIX}

# compile & install
make -j ${NJOBS}
make install




# we need to append these paths
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NETCDF_PREFIX}/lib
export LDFLAGS="$LDFLAGS -L${NETCDF_PREFIX}/lib"
export CPPFLAGS="$CPPFLAGS -I${NETCDF_PREFIX}/include"

# download NetCDF-Fortran
if [ ! -f "${DOWNLOADS_DIR}/$NETCDF_F-$NETCDF_F_VERSION.tar.gz" ]; then
curl -L \
            --output "${DOWNLOADS_DIR}/$NETCDF_F-$NETCDF_F_VERSION.tar.gz" \
            "https://github.com/Unidata/$NETCDF_F/archive/refs/tags/v$NETCDF_F_VERSION.tar.gz"
fi

# unpack
mkdir -p "${SRC_DIR}/$NETCDF_F" && cd "$_"
tar xvf "${DOWNLOADS_DIR}/$NETCDF_F-$NETCDF_F_VERSION.tar.gz"


# configure
mkdir -p "${SRC_DIR}/$NETCDF_F/build" && cd "$_"
FC=$MPIF90 CC=$MPICC CXX=$MPICXX ${SRC_DIR}/$NETCDF_F/$NETCDF_F-$NETCDF_F_VERSION/configure      \
    --with-pic                                                  \
    --prefix=${NETCDF_PREFIX}

# compile & install
make -j ${NJOBS}
make install
