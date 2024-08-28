#!/bin/bash
{
# start=0.025
# ending=0.200
# space=0.020
start=0.1
ending=1.2
space=0.1
# start=0.100
# ending=0.200
# space=0.020
name_all="pfaps_test48_harmonic_wall"
mode="pfap"
num_segments=100
input_density="pfaps_homo/pfaps_phase_diagram1_histo_phase_diagram"
init=homo
homo_fit="gauss"
num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
vars=$(seq $start $space $ending)
# vars="0.02 0.05 0.08 0.10 0.15 0.20 0.25 0.30 0.40 0.50"
vars="0.2 0.5 0.8 1.0 1.1 1.2 1.3 1.4"
num=$(echo "$vars" | wc -w)
slab_fit="average"
{ 
time {
# tail -n +2 "$input_density" | xargs -P9 -n 3 ./abp_pfaps_walls.sh $name_all
# printf "$vars" | xargs -P$num -d ' ' -I{} ./abp_pfaps_walls.sh $name_all {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps_walls.sh $name_all {}
# python3 analyze_veff.py $name_all $mode "$vars" $num_segments
python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 1
# python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit $slab_fit
} \
} 2> "$name_all"_output

exit
}