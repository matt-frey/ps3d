#!/bin/bash

# 29 June 2022
# https://stackoverflow.com/questions/6482377/check-existence-of-input-argument-in-a-bash-shell-script
if [ $# -eq 0 ]; then
      echo "Missing number of grid cells: bash $0 [ncells]"
      exit
fi


echo "Submit $1x$1x$1"
export OMP_NUM_THREADS=6; nohup ps3d --config "beltrami_$1.config" > "nohup_$1.out" &
