#!/bin/bash

echo "Running prediss scan ..."
for i in $(seq 10 10 100); do
    export OMP_NUM_THREADS=1; nohup ps3d --config "beltrami_32_pred$i.config" > "nohup_32_pred$i.out" &
done
