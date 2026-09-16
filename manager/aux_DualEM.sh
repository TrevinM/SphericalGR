#Returns 1 if the simulation doesn't meet tests

#Check for abrupt drop in lapse
i=0
curr_lapse=0
curr_delta=0
curr_t=0
while read -r line
do
    if [[ $line != \#* ]]
    then

        split=( $line )

        #Check that the new line isn't just a regridding
        new_t=${split[0]}
        if [[ 0 -eq "$(echo "$new_t == $curr_t" | bc)" ]]
        then
            new_lapse=${split[7]}

            # Check for scientific notation
            if [[ $new_lapse = *e-* ]]
            then
                new_lapse=0
            fi

            # Check for weird end behavior
            let "i++"
            new_delta="$(echo "$new_lapse - $curr_lapse" | bc)"

            if [[ i -gt 2 && 1 -eq "$(echo "$new_delta - $curr_delta < -0.03" | bc)" ]]
            then
                echo 1
                exit 0
            fi

            curr_delta=$new_delta
            curr_lapse=$new_lapse
            curr_t=$new_t
        fi
    fi
done        

echo 0
