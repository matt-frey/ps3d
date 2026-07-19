#!/bin/bash

# recipe for:
P=openmpi
V_MAJOR=4
V_MINOR=1
V=${V_MAJOR}.${V_MINOR}.5

# download
curl -L \
    --output "${DOWNLOADS_DIR}/$P-$V.tar.gz" \
    "https://download.open-mpi.org/release/open-mpi/v4.1/$P-$V.tar.gz"

# unpack
mkdir -p "${SRC_DIR}/$P" && cd "$_"
tar xvf "${DOWNLOADS_DIR}/$P-$V.tar.gz"

MPI_PREFIX=$PREFIX/openmpi

# configure
mkdir -p "${SRC_DIR}/$P/build" && cd "$_"
${SRC_DIR}/$P/$P-$V/configure       \
    --enable-mpi-cxx                \
    --enable-fortran                \
    --enable-shared                 \
    --enable-static                 \
    --without-verbs                 \
    --with-pic                      \
    --prefix=${MPI_PREFIX}

# compile & install
make -j ${NJOBS}
make install
