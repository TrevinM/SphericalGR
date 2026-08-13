#If any of the following variables are written into the example input file they will get replaced.
export a1=$eta
export a2=$eta
export a3=$eta

export as1=$eta
export as2=$eta
export as3=$eta

#Modify Files
echo "$(envsubst '$a1, $a2, $a3, $as1, $as2, $as3, $eta' < $it_path/$eta/DualEM_Input)" > $it_path/$eta/DualEM_Input
