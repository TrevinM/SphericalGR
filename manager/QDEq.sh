#Copy an example test that fits needs
cp $dir/test/example_DualMax/* $it_path/$ETA/
wait

export D=$ETA
export Q=$ETA

#Modify files
cat $dir/manager/QD_DualEM | envsubst '$D, $Q' > $it_path/$ETA/DualEM_Input
