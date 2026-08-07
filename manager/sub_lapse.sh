#Check for subcritical result using lapse above 0.85
#Returns 2 if lapse decreases by over 0.2 in one step

bug "Running sub_lapse.sh"

if [ -e output/DualMaxwell*.mon ]
then
    curr_lapse=0
    #Get last minimum lapse
    while read -r line
    do
        split=( $line )
        if [[ ${split[0]} != "#" ]]
        then
            new_lapse=${split[7]}
            #Check for weird end behavior
            if [[ 1 -eq "$(echo "$curr_lapse - $new_lapse > 0.2" | bc)" ]]
            then
                exit 2
            fi
            curr_lapse=$new_lapse
        fi
    done < output/DualMaxwell*.mon
    
    #Check Lapse
    if [[ 1 -eq "$(echo "$curr_lapse > 0.85" | bc)" ]]
    then
        exit 1
    fi

    exit 0
else
    log "ERROR: Couldn't open monitor..."
    exit -1
fi