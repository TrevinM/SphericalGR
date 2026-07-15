#!/bin/bash
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=16
if [ ! -d output ]; then mkdir output; fi

../../SphericalGR Input