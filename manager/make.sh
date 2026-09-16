# Bash this script to automatically make critical convergence test suites.

# ========================= #
#           INPUTS          #
# ========================= #


read -p "Suite: " -i suite

#Code directory
dir=/mnt/research/tbaumgar/Students/tmacomber/SphericalGR

#Directory to make suites in
suite_path=$dir/test/CRIT_$suite

if [[ ! -d $suite_path ]]
then
    #Create test suite
    mkdir $suite_path
    wait
    mkdir $suite_path/example
    wait

    #Get Settings
    read -p "Eta Family: " -i family
    read -e -p "Example: " -i "example_DualMax" example
    read -e -p "Lower Bound: " -i "0" lower_bound
    read -e -p "Upper Bound: " -i "10" upper_bound
    read -e -p "Max Decimals: " -i "12" max_precision
    read -p "Number of Cores: " num_cores
    read -e -p "Subcritical Test: " -i "sub_lapse.sh" sub_test
    read -e -p "Supercritical Test: " -i "sup_AH.sh" sup_test
    read -e -p "Fast Convergence (0/1): " -i "1" fast_conv

    #Copy Example
    cp $dir/test/$example/* $suite_path/example/

    #Write Settings
    echo $lower_bound $upper_bound      >> $suite_path/constraints
    echo Family     $family             >> $suite_path/settings
    echo Precision  $max_precision      >> $suite_path/settings
    echo Cores      $num_cores          >> $suite_path/settings
    echo SubTest    $sub_test           >> $suite_path/settings
    echo SupTest    $sup_test           >> $suite_path/settings
    echo FastConv   $fast_conv          >> $suite_path/settings

    echo "Test suite $suite created!"
    echo ""
    echo "!!! Don't forget to add variables to the suite Inputs !!!"
else
    echo "Suite ${suite} already created."
    echo "SBATCH launch to run the suite"
fi
