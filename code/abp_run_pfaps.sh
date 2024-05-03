#!/bin/bash

start=0.018
ending=0.038
space=0.002
name_all="pfaps_test18_diagram"
mode="pfap"
num_segments=10
input_density="pfaps_phase_diagram1_histo_phase_diagram"
fit="gauss"
{ 
time {
tail -n +2 "$input_density" | xargs -P11 -n 3 ./abp_pfaps.sh
# seq $start $space $ending | xargs -P11 -I{} ./abp_pfaps.sh {}

python3 analyze_phases.py $name_all $mode $start $ending $space $num_segments $fit
} \
} \
2> "$name_all"_output