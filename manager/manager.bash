#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --open-mode=append

# This code is to be run via the 'launch' script

export OMP_NUM_THREADS=16

if [[ $SLURM_JOB_ID == "" ]]
then
    echo "Cannot run locally. Use SBATCH..."
    exit 1
fi


# ========================== #
#         INITIALIZE         #
# ========================== #


bug () { [[ verbose -ge $2 ]] && [[ verbose -gt 0 ]] && echo "`date`" "    " $1 ; }
export -f bug
log () { echo $1 >> $suite_path/log; }
export -f log
finish () {
    #Remove empty manager Slurms
    cd $dir/manager
    if [ -e $report_file ]
    then
        slurm_empty=1
        while read -r line
        do
            if [[ ! -z $line ]]; then slurm_empty=0 ; fi
        done < "slurm-${SLURM_JOB_ID}.out"
        if [[ $slurm_empty -eq 1 ]]; then rm "slurm-${SLURM_JOB_ID}.out"; fi

    fi
    bug "Quitting Manager..."
    exit $1
}

#Read Settings
settings=( $(cat $suite_path/settings) )
family="${settings[1]}"
max_precision="${settings[3]}"
num_cores="${settings[5]}"
sub_test="${settings[7]}"
sup_test="${settings[9]}"
fast_conv="${settings[11]}"

bug "Suite:     $suite"    
bug "Family:    $family"
bug "Precision: $max_precision"
bug "Cores:     $num_cores"
bug "SubTest:   $sub_test"
bug "SupTest:   $sup_test"
bug "FastConv:  $fast_conv"
bug "Slurm ID:  $SLURM_JOB_ID"


num_data_rho=0
num_data_sup=0

#Count previous data
if [ -e $suite_path/data_rho ]
then
    while read -r line
    do
        let num_data_rho++
    done < $suite_path/data_rho
fi
if [ -e $suite_path/data_sup ]
then
    while read -r line
    do
        let num_data_sup++
    done < $suite_path/data_sup
fi


if [ ! -e $suite_path/SphericalGR ]
then
    cp $dir/build/SphericalGR $suite_path/
    bug "New SphericalGR"
fi


# =========================== #
#           EXECUTE           #
# =========================== #

#Find iteration and constraints
i=0
constraint=()

while read -r line
do
    if [[ ! -z $line ]]
    then
        let i++
        constraint=($line)
    fi

done < $suite_path/constraints

lower_bound=${constraint[0]}
upper_bound=${constraint[1]}

#Increase precision until bounds fully represented
for (( curr_precision=0 ; curr_precision<$max_precision ; curr_precision++ ))
do
    rounded_l="$(echo "scale=$curr_precision;$lower_bound/1" | bc)"
    rounded_u="$(echo "scale=$curr_precision;$upper_bound/1" | bc)"
    if [[ 1 -eq "$(echo "$rounded_l==$lower_bound" | bc)" && 1 -eq "$(echo "$rounded_u==$upper_bound" | bc)" ]]
    then
        #Add precision for delta
        let curr_precision++
        break
    fi
done

if [[ $reconverge -eq 0 ]]
then
    # Regular Iteration

    export it_path=$suite_path/iteration_$i
    cores=$num_cores

    if [ ! -d $it_path ]
    then
        mkdir $it_path
    else 
        echo "ERROR: Iteration $i already created but didn't add new constraints..."
        exit 1
    fi

    bug "========= Iteration $i ========="

    log ""
    log "========= Iteration $i ========="
    log "`date`"
    log ""


    #Set iteration variables
    delta="$(echo "scale=$curr_precision;($upper_bound-$lower_bound)/10" | bc)"
    delta_fast="$(echo "scale=$curr_precision+1;($upper_bound-$lower_bound)/100" | bc)"


    #Check for precision
    if [ "$delta" == "0" ]
    then
        log "ERROR: Delta = 0. Resolution may have exceeded decimal precision of $curr_precision..."
        finish 1
    fi
    
    log "Bounds: $lower_bound and $upper_bound"
    log "Previous data (rho): $num_data_rho"
    # log "Previous data (sup): $num_data_sup"
    conv_type=0

    #Estimate eta critical
    if [[ $num_data_rho -ge 5 ]]
    then
        cd $suite_path
        py_out="$(python3 $dir/manager/sub_fit.py $lower_bound $upper_bound)"
        py_out=( $py_out )
        eta_rho=${py_out[0]}
        eta_c=$eta_rho
        log "Eta (rho): $eta_rho"
        let conv_type++ 
    fi
    # if [[ $num_data_sup -ge 5 ]]
    # then
    #     cd $suite_path
    #     py_out="$(python3 $dir/manager/sup_fit.py $lower_bound $upper_bound)"
    #     py_out=( $py_out )
    #     eta_sup=${py_out[0]}
    #     eta_c=$eta_sup
    #     log "Eta (sup): $eta_sup"
    #     let conv_type++
    #     let conv_type++
    # fi
    if [[ $conv_type -eq 3 ]]
    then
        eta_avg="$(echo "scale=$max_precision;($eta_sup+$eta_sub)/2" | bc)"
        log "Eta (avg): $eta_avg"

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
    if [[ $conv_type -eq 0 || "$delta_fast" == "0" || $fast_conv -eq 0 ]]
    then
        #Simply divide constraints into 10 cells
        # log "Delta: $delta"
        log "Using standard convergence"

        for j in 1 2 3 4 5 6 7 8 9
        do
            etas=("${etas[@]}" $(echo "scale=$curr_precision;$lower_bound+$j*$delta" | bc))
        done
    elif [[ $conv_type -eq 1 || $conv_type -eq 2 || $conv_type -eq 3 ]]
    then
        #Use previous data to try and find critical value faster
        overflow=1

        #Add values around estimated critical
        for j in -10 -5 -2 -1 0 1 2 5 10
        do
            eta="$(echo "scale=$curr_precision+1;$eta_c+$j*$delta_fast" | bc)"
            if [[ 1 -eq "$(echo "$eta < $upper_bound" | bc)" && 1 -eq "$(echo "$eta > $lower_bound" | bc)" ]]
            then
                etas=("${etas[@]}" $eta)
            else
                let overflow++
            fi
        done

        # log "Delta2: $delta_fast"
        if [[ $conv_type -eq 1 ]]
        then
            log "Using fast (rho) convergence with $overflow overflow"
        elif [[ $conv_type -eq 2 ]]
        then
            log "Using fast (sup) convergence with $overflow overflow"
        else
            log "Using fast (avg) convergence with $overflow overflow"
        fi

        #Divide overflow tests evenly across constrints
        delta="$(echo "scale=$curr_precision+1;($upper_bound-$lower_bound)/($overflow + 1)" | bc)"
        for j in $(seq 1 $overflow)
        do
            eta="$(echo "scale=$curr_precision+1;$lower_bound+$j*$delta" | bc)"
            etas=("${etas[@]}" $eta)
        done
    else
        #Use previous data to try and find critical value faster
        overflow=0

        #Add values around eta sub
        for j in -4 -1 0 1 4
        do
            eta="$(echo "scale=$curr_precision+1;$eta_sub+$j*$delta_fast" | bc)"
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
            eta="$(echo "scale=$curr_precision+1;$eta_sup+$j*$delta_fast" | bc)"
            if [[ 1 -eq "$(echo "$eta < $upper_bound" | bc)" && 1 -eq "$(echo "$eta > $lower_bound" | bc)" ]]
            then
                etas=("${etas[@]}" $eta)
            else
                let overflow++
            fi
        done          
        
        log "Using dual convergence with $overflow overflow"
        
        #Divide overflow tests evenly across constrints
        if [[ overflow -gt 0 ]]
        then
            delta="$(echo "scale=$curr_precision+1;($upper_bound-$lower_bound)/($overflow + 1)" | bc)"
            for j in $(seq 1 $overflow)
            do
                eta="$(echo "scale=$curr_precision+1;$lower_bound+$j*$delta" | bc)"
                etas=("${etas[@]}" $eta)
            done
        fi
    fi
else

    # =========================== #
    #        RECONVERGENCE        #
    # =========================== #

    #Number the reconvergence
    for (( r=1 ; i<10 ; r++ ))
    do
        if [ ! -d $suite_path/reconverge_${r} ]
        then
            export it_path=$suite_path/reconverge_$r
            mkdir $it_path
            wait
            break
        fi
    done

    #Move over all old iteration to old directories
    for i in $( seq 0 15 )
    do
        if [[ -d ${suite_path}/iteration_$i ]]
        then
            if [[ ! -d ${suite_path}/old_$r ]]
            then
                mkdir ${suite_path}/old_$r
                wait
            fi
            mv ${suite_path}/iteration_$i ${suite_path}/old_$r/iteration_$i
        fi
    done

    bug "========= Reconverge $r ========="

    log ""
    log "========= Reconverge $r ========="
    log "`date`"
    log ""

    #Really only matters for naming of jobs
    i="r${r}"

    #Run Lower bounds at very precisions and +1 to the last digit
    etas=(10)
    delta=10
    let curr_precision--
    for s in $( seq 0 $curr_precision )
    do
        delta="$(echo "scale=$s;$delta/10" | bc)"
        eta_l="$(echo "scale=$s;$lower_bound/1" | bc)"
        eta_u="$(echo "scale=$s;$eta_l+$delta" | bc)"
        [[ $s -eq 0 ]] && log "Bounds: $eta_l and $eta_u"
        etas=("${etas[@]}" $eta_l)
        etas=("${etas[@]}" $eta_u)
    done
    upper_bound=100
    lower_bound=0

    #Define cores to not exceed regular core usage (and be even)
    let "cores=($num_cores*9)/(${#etas[@]}*2)"
    let "cores=$cores*2"
    conv_type=-1

fi


# ========================== #
#          RUN JOBS          #
# ========================== #

bug "Etas: ${etas[*]}"

#Create all jobs
job_dirs=()
jobs=()
for eta in "${etas[@]}"
do
    if [ ! -d $it_path/$eta ]
    then
        mkdir $it_path/$eta
        wait
        cp $dir/manager/moosehead.sh $it_path/$eta/${suite}_${i}_${eta}.sh
        cp $suite_path/example/* $it_path/$eta/
        wait

        #Copy and edit input files
        export eta=$eta
        bash $dir/manager/eta_$family.sh

        #Add job
        job_dirs=("${job_dirs[@]}" $it_path/$eta)
        jobs=("${jobs[@]}" ${suite}_${i}_${eta}.sh)
    else
        log "WARNING: The test $suite\_$i\_$eta had already been created"
    fi
done

#Start all jobs
SECONDS=0

let "len_jobs = ${#jobs[@]} - 1"
for k in $( seq 0 $len_jobs )
do
    cd ${job_dirs[k]}
    sbatch -W -N 1 -n $cores ${jobs[k]} >> $suite_path/log &
done
wait

log "Finished in $(echo "scale=2;$SECONDS/3600" | bc) hours on $cores cores"
bug "Finished in $(echo "scale=2;$SECONDS/3600" | bc) hours on $cores cores"


# ========================== #
#       CRITICAL TESTS       #
# ========================== #


upper_new=$upper_bound
lower_new=$lower_bound

missing_etas=()
weird_etas=()
uncon_etas=()
dual_etas=()

for k in $( seq 0 $len_jobs )
do
    cd ${job_dirs[k]}
    eta=${etas[k]}

    # Check for a monitor file
    if [ ! -e **/*.mon ]
    then
        # Couldn't Open a monitor
        bug "$eta missing monitor in ${job_dirs[k]}..."
        missing_etas=("${missing_etas[@]}" $eta)
        continue
    fi

    # # Check for bad end behaviors
    # if [ -e $dir/manager/test_$family ]
    # then
    #     test="$( bash $dir/manager/aux_$family )"
    #     if [[ $test -eq 1 ]]
    #     then
    #         #Weird End Behavior
    #         bug "$eta weird end behavior..."
    #         weird_etas=("${weird_etas[@]}" $eta)s
    #         continue
    #     fi
    # else
    #     bug "No aux test found for family $family"
    # fi

    # Sub or Super critical
    #  1 : Test passed
    #  0 : Test Failed
    # -1 : Error

    sub_crit="$( bash $dir/manager/$sub_test )"
    sup_crit="$( bash $dir/manager/$sup_test )"

    if [[ $sub_crit -eq -1 || $sup_crit -eq -1 ]]
    then
        bug "$eta failed tests..."
        missing_etas=("${missing_etas[@]}" $eta)
        continue

    elif [[ $sub_crit -eq 1 && $sup_crit -eq 1 ]]
    then
        #Both sub and sup crit
        bug "$eta dual behavior..."
        dual_etas=("${dual_etas[@]}" $eta)
        continue

    elif [[ $sub_crit -eq 1 ]]
    then
        #Sub critical
        bug "$eta sub critical"
        [[ 1 -eq "$(echo "$eta > $lower_new" | bc)" ]] && lower_new=$eta && bug "changed lower bound"

        #Log Data from .dualmaxwell_mon
        if [ -e **/DualMaxwell*.dualmaxwell_mon ]
        then
            while read -r line
            do
                if [[ $line != \#* ]]
                then
                    dual_split=( $line )
                fi
            done < **/DualMaxwell*.dualmaxwell_mon
                   
            echo "$(echo "scale=$max_precision;$eta/1" |bc)" ${dual_split[5]} >> $suite_path/data_rho
            let num_data_rho++
        else
            log "WARNING: Couldn't open dualmaxwell $eta monitor..." 
            bug "WARNING: Couldn't open dualmaxwell $eta monitor..." 
        fi

    elif [[ $sup_crit -eq 1 ]]
    then
        #Sup critical
        bug "$eta SUPER critical"
        [[ 1 -eq "$(echo "$eta < $upper_new" | bc)" ]] && upper_new=$eta && bug "changed upper bound"

    else
        #Neither
        bug "$eta neither"
        uncon_etas=("${uncon_etas[@]}" $eta)
        continue
    fi

    #Log general data
    while read -r line
    do
        if [[ $line != \#* ]]
        then
            mon_split=( $line )
        fi
    done < **/*.mon

    echo "$(echo "scale=$max_precision;$eta/1" |bc)" ${mon_split[0]} ${mon_split[1]} >> $suite_path/data_time

done


if [[ $missing_etas != "" ]]
then
    log "ERROR: Etas ${missing_etas[*]} did not have monitor files..."
fi
if [[ $weird_etas != "" ]]
then
    log "ERROR: Etas ${weird_etas[*]} weird end behavior..."
fi
if [[ $dual_etas != "" ]]
then
    log "ERROR: Etas ${dual_etas[*]} dual end behavior..."
fi
if [[ $uncon_etas != "" ]]
then
    log "ERROR: Etas ${uncon_etas[*]} unconclusive..."
fi
if [[ $missing_etas != "" || $weird_etas != "" || $dual_etas != "" || $uncon_etas != "" ]]
then
    finish 1
fi


# =========================== #
#      VERIFY NEW BOUNDS      #
# =========================== #


#Check for Convergence
if [[ "$upper_new" = "$upper_bound" && "$lower_new" = "$lower_bound" ]]
then
    log "ERROR: Did not converge..."
    finish 1
fi

#Add warning if one of the boundaries didn't move
if [[ "$upper_new" = "$upper_bound" ]]
then
    log "WARNING: All tests were subcritical"
    let warning_prev++
elif [[ "$lower_new" = "$lower_bound" ]]
then
    log "WARNING: All tests were supercritical"
    let warning_prev++
else
    warning_prev=0
fi

#Print out specifics on new constraints
if [[ $conv_type -eq 0 || "$delta_fast" == "0" || $fast_conv -eq 0 ]]
then
    log "New constraints: $(echo "scale=1;($lower_new - $lower_bound)/$delta" | bc) and $(echo "scale=1;($upper_new - $lower_bound)/$delta" | bc)"
fi
if [ $conv_type -eq 1 ]
then
    log "New constraints: $(echo "scale=1;($lower_new - $eta_sub)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sub)/$delta_fast" | bc) (rho)"
elif [ $conv_type -eq 2 ]
then
    log "New constraints: $(echo "scale=1;($lower_new - $eta_sup)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sup)/$delta_fast" | bc) (sup)"
elif [ $conv_type -ge 3 ]
then
    log "New constraints: $(echo "scale=1;($lower_new - $eta_sub)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sub)/$delta_fast" | bc) (rho)"
    log "                 $(echo "scale=1;($lower_new - $eta_sup)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_sup)/$delta_fast" | bc) (sup)"
    log "                 $(echo "scale=1;($lower_new - $eta_avg)/$delta_fast" | bc) and $(echo "scale=1;($upper_new - $eta_avg)/$delta_fast" | bc) (avg)"

fi

if [[ 1 -eq "$(echo "$upper_new <= $lower_new" | bc)" ]]
then
    log "ERROR: Bounds flipped..."
    finish 1
fi

#Print out convergence factor
convergence=$(echo "($upper_bound-$lower_bound)/($upper_new-$lower_new)" | bc)
log "Converged by ${convergence}x"

if [ $reconverge -eq 0 ]
then
    echo $lower_new $upper_new >> $suite_path/constraints
else
    echo $lower_new $upper_new > $suite_path/constraints
    reconverge=0
fi

if [[ $warning_prev -ge 3 ]]
then
    log "ERROR: 3 Warning Strikes..."
    finish 1
fi

if [[ $max_runs -gt 1 ]]
then
    let max_runs--
    sbatch -J "${suite}_manager" -o $report_file $dir/manager/manager.bash
else
    log "Maximum number of runs ($max_runs) for 1 job has been reached..."
fi

finish 0
