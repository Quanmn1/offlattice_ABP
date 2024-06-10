#!/bin/bash

start=0.04
ending=0.15
space=0.01
name_all="pfaps_qsaps_test4"
mode="pfqs"
num_segments=6
init="homo"
# input_density="pfaps_phase_diagram1_histo_phase_diagram"
fit="max"
{ 
time {
# tail -n +2 "$input_density" | xargs -P11 -n 3 ./abp_pfaps.sh
# seq $start $space $ending | xargs -P12 -I{} ./abp_pfaps_qsaps.sh {}
python3 analyze_phases.py $name_all $mode $start $ending $space $num_segments $init $fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode $start $ending $space
} \
} \
2> "$name_all"_output