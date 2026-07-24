#Copy an example test that fits needs
cp $dir/test/example_DualMax/* $it_path/$ETA/
wait

#Using optimization of J/E to determine relation
export D=$ETA
export Q="$(echo "scale=16;2.5354627641855497*$ETA" | bc)"

#Modify files
cat $dir/manager/QD_DualEM | envsubst '$D, $Q' > $it_path/$ETA/DualEM_Input
