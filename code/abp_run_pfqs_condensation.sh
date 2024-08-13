#!/bin/bash

# start=5
# ending=60
# space=5
# start=0.0
# ending=1.0
# space=0.2
# start=0.14
# ending=0.30
# space=0.04
name_all="pfqs_condensation_test3_scaleeps"
mode="pfqs"
num_segments=5
init="homo"
# input_density="pfaps_phase_diagram1_histo_phase_diagram"
homo_fit="gauss"
slab_fit="average"
# vars=$(seq $start $space $ending)
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
vars="0.08 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32"
num=$(echo $vars | wc -w)

{ 
time {
# tail -n +2 "$input_density" | xargs -P11 -n 3 ./abp_pfaps.sh
# printf "$vars" | xargs -P$num -d ' ' -I{} ./abp_pfqs_condensation.sh {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfqs_condensation.sh {}
python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode $start $ending $space
# python3 analyze_sigma.py $name_all $mode $start $ending $space
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} \
2> "$name_all"_output