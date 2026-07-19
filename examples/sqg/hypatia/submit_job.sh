#!/bin/bash
#SBATCH --job-name=ps3d
#SBATCH --partition=large-long
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --output=%x_%j.log

# Load required modules
module load openmpi/5.0.7
module load ucx/1.16.0
module load hdf5/1.14.6

export UCX_WARN_UNUSED_ENV_VARS=n

echo "PREFIX: ${PREFIX}"

export LD_LIBRARY_PATH=${PREFIX}/netcdf/lib:$LD_LIBRARY_PATH

PS3D_BIN=${PREFIX}/ps3d/bin

srun --mpi=pmix_v3 $PS3D_BIN/ps3d --config ps3d_sqg.config
