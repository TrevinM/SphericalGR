#!/bin/bash
# SBATCH --mail-type=BEGIN,END,FAIL

export OMP_NUM_THREADS=16
if [ ! -d output ]; then mkdir output; fi
rm SphericalGR
cp ../../build/SphericalGR .

./SphericalGR Input
