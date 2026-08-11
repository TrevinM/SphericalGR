nr_list=( 128 256 384 512 768 1024 )
nt_list=( 8 16 24 32 48 64 )

#If any of the following variables are written into the example input file they will get replaced.
export nr=${nr_list[$eta]}
export nt=${nt_list[$eta]}
export sp=0
export rf=0

#Modify Files
echo "$(envsubst '$nr, $nt, $sp, $rf, $eta' < $it_path/$eta/Grid_Input)" > $it_path/$eta/GridInput
