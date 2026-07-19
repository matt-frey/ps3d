# Set up environment for Open MPI with gnu compilers
module load openmpi/5.0.7
module load ucx/1.16.0
module load hdf5/1.14.6

export UCX_WARN_UNUSED_ENV_VARS=n

export DOWNLOADS_DIR=/home/mf248/work/downloads
export SRC_DIR=/home/mf248/work/sources
export PREFIX=/home/mf248/work/gcc/11.4.1
