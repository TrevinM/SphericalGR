
#Modify Files
echo "$(envsubst '$eta' < $it_path/$eta/DualEM_Input)" > $it_path/$eta/DualEM_Input
