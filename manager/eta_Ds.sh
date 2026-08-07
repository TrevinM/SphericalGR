export a1=0.0
export a2=0.0
export a3=0.0

export as1=$eta
export as2=0.0
export as3=0.0

#Modify Files
cat $dir/manager/edit_DualEM | envsubst '$a1, $a2, $a3, $as1, $as2, $as3' > $it_path/$eta/DualEM_Input
