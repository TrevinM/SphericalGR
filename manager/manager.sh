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
decimal_precision=12


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
num_data_sub=0
num_data_sup=0
if [ -e $suite_path/data_sup ]
then
    while read -r line
    do
        let num_data_sup++
    done < $suite_path/data_sup
fi


# ======================= #
#         EXECUTE         #
# ======================= #

for (( n=0 ; n<14 ; n++ ))
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
        echo "ERROR: Iteration $i already created but didn't add new constraints..."
        exit 1
    fi
    
    echo "Bounds: $lower_bound and $upper_bound"
    echo "Previous data (sub): $num_data_sub"
    echo "Previous data (sup): $num_data_sup"
    conv_type=0

    #Estimate eta critical
    if [[ $num_data_sub -ge 5 ]]
    then
        cd $suite_path
        py_out="$(python3 $dir/manager/sub_fit.py $lower_bound $upper_bound)"
        py_out=( $py_out )
        eta_sub=${py_out[0]}
        eta_c=$eta_sub
        echo "Eta (sub): $eta_sub"
        let conv_type++ 
    fi
    if [[ $num_data_sup -ge 5 ]]
    then
        cd $suite_path
        py_out="$(python3 $dir/manager/sup_fit.py $lower_bound $upper_bound)"
        py_out=( $py_out )
        eta_sup=${py_out[0]}
        eta_c=$eta_sup
        echo "Eta (sup): $eta_sup"
        let conv_type++
        let conv_type++
    fi
    if [[ $conv_type -eq 3 ]]
    then
        eta_avg="$(echo "scale=$decimal_precision;($eta_sup+$eta_sub)/2" | bc)"
        if [[ 1 -eq "$(echo "($eta_sup-$eta_sub)/$delta_fast < 2*$delta_fast" | bc)" && 1 -eq "$(echo "($eta_sub-$eta_sup)/$delta_fast < 2*$delta_fast" | bc)" ]]
        then
            #Average convergence. eta_sup and eta_sub close together
            eta_c=$eta_avg
        else
            #Dual convergence.
            conv_type=4
        fi
    fi

    #Calculate etas
    etas=()
    if [[ $conv_type -eq 0 || "$delta_fast" == "0" ]]
    then
        #Simply divide constraints into 10 cells
        echo "Delta: $delta"
        echo "Using standard convergence"

        for j in 1 2 3 4 5 6 7 8 9
        do
            etas=("${etas[@]}" $(echo "scale=$decimal_precision;$lower_bound+$j*$delta" | bc))
        done
    elif [[ $conv_type -eq 1 || $conv_type -eq 2 || $conv_type -eq 3 ]]
    then
        #Use previous data to try and find critical value faster
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

        echo "Delta2: $delta_fast"
        if [[ $conv_type -eq 1 ]]
        then
            echo "Using fast (sub) convergence with $overflow overflow"
        elif [[ $conv_type -eq 2 ]]
        then
            echo "Using fast (sup) convergence with $overflow overflow"
        else
            echo "Using fast (avg) convergence with $overflow overflow"
        fi

        #Divide overflow tests evenly across constrints
        delta="$(echo "scale=$decimal_precision;($upper_bound-$lower_bound)/($overflow + 1)" | bc)"
        for j in $(seq 1 $overflow)
        do
            eta="$(echo "scale=$decimal_precision;$lower_bound+$j*$delta" | bc)"
            etas=("${etas[@]}" $eta)
        done
    else
        #Use previous data to try and find critical value faster
        overflow=0

        #Add values around eta sub
        for j in -4 -1 0 1 4
        do
            eta="$(echo "scale=$decimal_precision;$eta_sub+$j*$delta_fast" | bc)"
            if [[ 1 -eq "$(echo "$eta < $upper_bound" | bc)" && 1 -eq "$(echo "$eta > $lower_bound" | bc)" ]]
            then
                etas=("${etas[@]}" $eta)
            else
                let overflow++
            fi
        done    

        #Add values around eta sup
        for j in -4 -1 0 1 4
        do
            eta="$(echo "scale=$decimal_precision;$eta_sup+$j*$delta_fast" | bc)"
            if [[ 1 -eq "$(echo "$eta < $upper_bound" | bc)" && 1 -eq "$(echo "$eta > $lower_bound" | bc)" ]]
            then
                etas=("${etas[@]}" $eta)
            else
                let overflow++
            fi
        done          
        
        echo "Using dual convergence with $overflow overflow"
        
        #Divide overflow tests evenly across constrints
        if [[ overflow -gt 0 ]]
        then
            delta="$(echo "scale=$decimal_precision;($upper_bound-$lower_bound)/($overflow + 1)" | bc)"
            for j in $(seq 1 $overflow)
            do
                eta="$(echo "scale=$decimal_precision;$lower_bound+$j*$delta" | bc)"
                etas=("${etas[@]}" $eta)
            done
        fi
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

    if [[ "$upper_new" != "$upper_bound" || "$lower_new" != "$lower_bound" ]]
    then
        echo $lower_new $upper_new >> $suite_path/constraints
    else
        echo "ERROR: Did not converge..."
        exit 1
    fi

    if [ $conv_type -eq 0 ]
    then
        echo "New constraints: $(echo "scale=1;($lower_new - $lower_bound)/$delta" | bc) and $(echo "scale=1;($upper_new - $lower_bound)/$delta" | bc)"
    elif [ $conv_type -eq 1 ]
    then
        echo "New constraints: $(echo "scale=1;($lower_new - $eta_sub)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sub)/$delta_fast" | bc) (sub)"
    elif [ $conv_type -eq 2 ]
    then
        echo "New constraints: $(echo "scale=1;($lower_new - $eta_sup)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sup)/$delta_fast" | bc) (sup)"
    else
        echo "New constraints: $(echo "scale=1;($lower_new - $eta_sub)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sub)/$delta_fast" | bc) (sub)"
        echo "                 $(echo "scale=1;($lower_new - $eta_sup)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sup)/$delta_fast" | bc) (sup)"
        echo "                 $(echo "scale=1;($lower_new - $eta_avg)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_avg)/$delta_fast" | bc) (avg)"

    fi

    convergence=$(echo "($upper_bound-$lower_bound)/($upper_new-$lower_new)" | bc)
    echo "Converged by ${convergence}x"

done

echo ""
echo "Maximum number of runs for 1 job has been reached..."