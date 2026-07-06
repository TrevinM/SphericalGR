#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL

export OMP_NUM_THREADS=16

./SphericalGR Input
