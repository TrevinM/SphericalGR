#!/bin/bash
#SBATCH --mail-type=FAIL


# ========================= #
#           INPUTS          #
# ========================= #

#Name of critical collapse suite to be run (with sbatch) or created (with bash)
export suite=$1

export max_runs=$2

export verbose=$3

export reconverge=$4 #True or false (0/1)

# ======================== #
#         Constants        #
# ======================== #

#Code directory
export dir=/mnt/research/tbaumgar/Students/tmacomber/SphericalGR

#Directory to make suites in
export suite_path=$dir/test/CRIT_$suite

#File to store bug reports
export report_file=$suite_path/output

# ======================== #
#           Setup          #
# ======================== #

if [ ! -d $suite_path ]
then
    echo "Need to first create the test suite by using make_suite file"
    exit 1
fi

warning_prev=0
export warning_prev

if [[ $SLURM_JOB_ID != "" ]]
then
    sbatch -J "${suite}_manager" -o $report_file $dir/manager/manager.bash
fi 

# ========================= #
#          Clean up         #
# ========================= #

#Remove empty launch Slurms
if [ -e "slurm-${SLURM_JOB_ID}.out" ]
then
    slurm_empty=1
    while read -r line
    do
        if [[ ! -z $line && ! $line =~ (^Submitted batch job )([0-9]+) ]]; then slurm_empty=0 ; fi
    done < "slurm-${SLURM_JOB_ID}.out"
    if [[ $slurm_empty -eq 1 ]]; then rm "slurm-${SLURM_JOB_ID}.out"; fi
fi
