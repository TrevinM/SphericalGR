#Check for subcritical result using lapse above 0.75

while read -r line
do
    if [[ $line != \#* ]]
    then
        mon_split=( $line )
    fi
done < **/*.mon
eta=${mon_split[7]}

#Scientific Notation Check
if [[ $eta = -*e* || $eta = *e-* ]]
then
    eta=0
elif [[ $eta = e*+* ]]
then
    eta=1000
fi

if [[ 1 -eq "$(echo "$eta > 0.75" | bc)" ]]
then
    echo 1
else
    echo 0
fi
