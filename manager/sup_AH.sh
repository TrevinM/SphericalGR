#Check if any horizons were found

if [ -e **/*.hor_mon ]
then
    #Check for non header lines
    while read -r line
    do
        if [[ ! $line = \#* ]]
        then
            exit 1
            break
        fi
    done < **/*.hor_mon

    exit 0
else
    exit -1
fi