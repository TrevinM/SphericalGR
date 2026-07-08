#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL

export OMP_NUM_THREADS=16

/mnt/research/tbaumgar/Students/pharris/github/build/SphericalGR /mnt/research/tbaumgar/Students/pharris/github/build/Input
