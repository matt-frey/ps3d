#!/bin/bash

# recipe for:
H5=hdf5
H5_MAJOR=1
H5_MINOR=12
H5_VER=${H5_MAJOR}.${H5_MINOR}.1

HDF5_PREFIX=$PREFIX/hdf5

mkdir -p $HDF5_PREFIX

#download HDF5
if [ ! -f "${DOWNLOADS_DIR}/$H5-$H5_VER.tar.gz" ]; then
    curl -L \
        --output "${DOWNLOADS_DIR}/$H5-$H5_VER.tar.gz" \
        "https://support.hdfgroup.org/ftp/HDF5/prev-releases/$H5-$H5_MAJOR.$H5_MINOR/$H5-$H5_VER/src/$H5-$H5_VER.tar.gz"
fi

# unpack
mkdir -p "${SRC_DIR}/$H5" && cd "$_"
tar xvf "${DOWNLOADS_DIR}/$H5-$H5_VER.tar.gz"

# guess system type
system_type="$(${SRC_DIR}/$H5/$H5-$H5_VER/bin/config.guess)"

# configure
mkdir -p "${SRC_DIR}/$H5/build" && cd "$_"
FC=$MPIF90 CC=$MPICC CXX=$MPICXX ${SRC_DIR}/$H5/$H5-$H5_VER/configure        \
    --build="$system_type"                  \
    --host="$system_type"                   \
    --target="$system_type"                 \
    --enable-fortran                        \
    --enable-shared                         \
    --enable-parallel                       \
    --with-pic                              \
    --prefix=${HDF5_PREFIX}

# compile & install
make -j ${NJOBS}
make install
