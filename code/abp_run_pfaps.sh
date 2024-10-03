#!/bin/bash
{
# start=0.025
# ending=0.200
# space=0.020
start=0.1
ending=1.4
space=0.1
# start=0.100
# ending=0.200
# space=0.020
name_all="pfaps_test44_harmonic_measure_largeeps"
mode="pfap"
num_segments=390
# input_density="pfaps_homo/pfaps_phase_diagram1_histo_phase_diagram"
# input_density="pfaps_homo/pfaps_phase_diagram1_histo_phase_diagram"
init=homo
homo_fit="gauss"
vars=$(seq $start $space $ending)
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="0.1 0.3 0.5 0.7 0.9 1.0 1.1 1.2 1.3"
# vars="0.02 0.04 0.05 0.06 0.07 0.08 0.10 0.12 0.15 0.20 0.25 0.30"
num=$(echo "$vars" | wc -w)
slab_fit="average"

{ 
time {
# tail -n +2 "$input_density" | xargs -P9 -n 3 ./abp_pfaps.sh
# printf "$vars" | xargs -P$num -d ' ' -I{} ./abp_pfaps.sh $name_all {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps.sh $name_all {}
# python3 analyze_veff.py $name_all $mode "$vars" $num_segments
python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaA"
python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaIK"
# python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit $slab_fit
} \
} 2> "$name_all"_output

exit
}