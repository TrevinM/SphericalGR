#!/bin/bash
# SBATCH --mail-type=END,FAIL

export OMP_NUM_THREADS=16


# ===================== #
#         INPUTS        #
# ===================== #

#Location of code directory
export dir=/mnt/research/tbaumgar/Students/tmacomber/SphericalGR

#Name of critical collapse suite to be run
suite=QxDO
suite_path=$dir/test/$suite

#Default upper and lower bounds
lower_bound=0
upper_bound=100

#Resolution of float math
decimal_precision=10


# ====================== #
#       INITIALIZE       #
# ====================== #

#Create test suite
if [ ! -d $suite_path ]
then
    mkdir $suite_path
    wait
    echo $lower_bound $upper_bound >> $suite_path/constraints
    cp $dir/build/SphericalGR $suite_path/
fi

echo $suite

#Count previous data
num_data=0
if [ -e $suite_path/data ]
then
    while read -r line
    do
        let num_data++
    done < $suite_path/data
fi


# ======================= #
#         EXECUTE         #
# ======================= #

for (( n=0 ; n<15 ; n++ ))
do
    #Find iteration and constraints
    i=-1
    constraint=()

    while read -r line
    do
        let i++
        constraint=($line)

    done < $suite_path/constraints

    echo ""
    echo "========= Iteration $i ========="
    echo `date`
    echo ""


    #Set iteration variables
    export it_path=$suite_path/iteration_$i

    lower_bound=${constraint[0]}
    upper_bound=${constraint[1]}

    delta="$(echo "scale=$decimal_precision;($upper_bound-$lower_bound)/10" | bc)"
    delta_fast="$(echo "scale=$decimal_precision;($upper_bound-$lower_bound)/100" | bc)"


    #Check for errors
    if [ "$delta" == "0" ]
    then
        echo "ERROR: Delta = 0. Resolution may have exceeded decimal precision of $decimal_precision..."
        exit 1
    fi

    if [ ! -d $it_path ]
    then
        mkdir $it_path
    else 
        echo "ERROR: Iteration $i already created but didn't yield new constraints..."
        exit 1
    fi
    
    echo "Bounds: $lower_bound and $upper_bound"
    echo "Previous data: $num_data"

    #Estimate eta critical
    if [[ $num_data -lt 5 ]]
    then
        fast_convergence=0
    else
        cd $suite_path
        eta_c="$(python3 $dir/manager/power_fit.py $lower_bound $upper_bound)"
        eta_c="$(echo "scale=$decimal_precision;$eta_c/1" | bc)"
        echo "Eta critical: $eta_c"
        fast_convergence=1
    fi

    #Calculate etas
    etas=()
    if [[ $fast_convergence -eq 0 || "$delta_fast" == "0" ]]
    then
        #Simply divide constraints into 10 cells
        echo "Using standard convergence"

        for j in 1 2 3 4 5 6 7 8 9
        do
            etas=("${etas[@]}" $(echo "scale=$decimal_precision;$lower_bound+$j*$delta" | bc))
        done
    else
        #Use previous black hole data to try and find critical value faster
        overflow=1

        #Add values around estimated critical
        for j in -10 -5 -2 -1 0 1 2 5 10
        do
            eta="$(echo "scale=$decimal_precision;$eta_c+$j*$delta_fast" | bc)"
            if [[ 1 -eq "$(echo "$eta < $upper_bound" | bc)" && 1 -eq "$(echo "$eta > $lower_bound" | bc)" ]]
            then
                etas=("${etas[@]}" $eta)
            else
                let overflow++
            fi
        done

        #Divide overflow tests evenly across constrints
        echo "Using fast convergence with $overflow overflow"

        delta="$(echo "scale=$decimal_precision;($upper_bound-$lower_bound)/($overflow + 1)" | bc)"
        for j in $(seq 1 $overflow)
        do
            eta="$(echo "scale=$decimal_precision;$lower_bound+$j*$delta" | bc)"
            etas=("${etas[@]}" $eta)
        done
    fi

    #Create all jobs
    job_dirs=()
    jobs=()
    for eta in "${etas[@]}"
    do
        if [ ! -d $it_path/$eta ]
        then
            mkdir $it_path/$eta
            wait
            cp $dir/manager/moosehead.sh $it_path/$eta/$suite\_$i\_$eta.sh
            cp $dir/test/example_DualMax/* $it_path/$eta/
            wait

            #Edit Inputs
            export ETA=$eta
            bash $dir/manager/$suite.sh

            #Add job
            job_dirs=("${job_dirs[@]}" $it_path/$eta)
            jobs=("${jobs[@]}" $suite\_$i\_$eta.sh)
        else
            echo "The test $suite\_$i\_$eta had already been created"
        fi
    done

    #Start all jobs
    SECONDS=0

    let "len_jobs = ${#jobs[@]} - 1"
    for k in $( seq 0 $len_jobs )
    do
        cd ${job_dirs[k]}
        sbatch -W -N 1 -n 4 ${jobs[k]} &
    done
    wait

    echo "Finished in $(echo "scale=2;$SECONDS/3600" | bc) hours"


    # ======================= #
    #     NEW CONSTRAINTS     #
    # ======================= #

    upper_new=$upper_bound
    lower_new=$lower_bound

    if [ $fast_convergence -eq 1 ]
    then
        echo "New constraints found at $(echo "scale=1;($lower_new - $eta_c)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_c)/$delta_fast" | bc) fast delta"
    fi

    convergence=$(echo "($upper_bound-$lower_bound)/($upper_new-$lower_new)" | bc)
    echo "Converged by ${convergence}x"

    if [[ "$upper_new" != "$upper_bound" || "$lower_new" != "$lower_bound" ]]
    then
        echo $lower_new $upper_new >> $suite_path/constraints
    else
        echo "ERROR: Did not converge..."
        exit 1
    fi
done

echo ""
echo "Maximum number of iterations for 1 job has been reached..."