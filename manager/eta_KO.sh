#If any of the following variables are written into the example input file they will get replaced.
export KO=$eta

echo "$(envsubst '$KO, $eta' < $it_path/$eta/Input)" > $it_path/$eta/Input
