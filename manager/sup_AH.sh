#Check if any horizons were found
bug "Running sup_AH.sh"

if [ -e output/DualMaxwell*.hor_mon ]
then
    #Check for non header lines
    while read -r line
    do
        if [[ ! $line = \#* ]]
        then
            exit 1
            break
        fi
    done < output/DualMaxwell*.hor_mon

    exit 0
else
    log "ERROR: Couldn't open horizon monitor..."
    exit -1
fi