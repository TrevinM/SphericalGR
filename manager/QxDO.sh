#Copy an example test that fits needs
cp $dir/test/example_DualMax/* $it_path/$ETA/
wait

#Modify files
cat $dir/manager/QxDO_DualEM | envsubst '$ETA' > $it_path/$ETA/DualEM_Input
