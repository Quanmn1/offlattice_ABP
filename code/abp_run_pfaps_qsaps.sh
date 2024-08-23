#!/bin/bash

# start=5
# ending=60
# space=5
# start=0.10
# ending=0.40
# space=0.02
start=0.06
ending=0.14
space=0.02
name_all="pfaps_qsaps_test16_harmonic_largescaleeps"
mode="pfqs"
num_segments=5
init="homo"
input_density="pfaps_qsaps_test15_harmonic_scaleeps_histo_phase_diagram"
homo_fit="gauss"
slab_fit="average"
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="0.10 0.12 0.13 0.14 0.16 0.18 0.20 0.22 0.24"
# vars=$(seq $start $space $ending)
num=$(echo $vars | wc -w)
{ 
time {
# tail -n +2 "$input_density" | awk '{print $1, $2, $3}' | xargs -P$num -n 3 ./abp_pfaps.sh $name_all
# printf "$vars" | xargs -P$num -d ' ' -I {} ./abp_pfaps_qsaps.sh {} $name_all
seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps_qsaps.sh {} $name_all
python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit $slab_fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode $start $ending $space
# python3 analyze_sigma.py $name_all $mode $start $ending $space
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} \
2> "$name_all"_output