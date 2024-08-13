#!/bin/bash

start=0.7
ending=0.9
space=0.1
name_all="qsaps_condensation_test0"
mode="qsap"
num_segments=5
init="homo"
fit="max"
{
time {

# seq $start $space $ending | xargs -P6 -I{} ./abp_qsaps_zero.sh {} 

python3 analyze_phases.py $name_all $mode $start $ending $space $num_segments $init $fit

} \
} \
2> "$name_all"_output