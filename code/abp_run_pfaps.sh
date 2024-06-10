#!/bin/bash

start=0.024
ending=0.032
space=0.004
name_all="pfaps_test23_diagram"
mode="pfap"
num_segments=10
input_density="pfaps_homo/pfaps_phase_diagram1_histo_phase_diagram"
init=homo
fit="max"
{ 
time {
# tail -n +2 "$input_density" | xargs -P3 -n 3 ./abp_pfaps.sh
seq $start $space $ending | xargs -P3 -I{} ./abp_pfaps.sh {}

python3 analyze_phases.py $name_all $mode $start $ending $space $num_segments $init $fit
} \
} \
2> "$name_all"_output