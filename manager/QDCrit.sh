#Copy an example test that fits needs
cp $dir/test/example_DualMax/* $it_path/$ETA/
wait

#Using the factor between eta critical for dipole and quadropole found by Mendoza
export D=$ETA
export Q="$(echo "scale=16;3.87031907035*$ETA" | bc)"

#Modify files
cat $dir/manager/QD_DualEM | envsubst '$ETA, $Q' > $it_path/$ETA/DualEM_Input
