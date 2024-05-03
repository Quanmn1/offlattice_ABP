#!/bin/bash

start=20
ending=100
space=10
name_all="qsaps_test10"
mode="qsap"
num_segments=5

{
time {

seq $start $space $ending | xargs -P9 -I{} ./abp_qsaps_exp.sh {} 

python3 analyze_phases.py $name_all $mode $start $ending $space $num_segments

} \
} \
2> "$name_all"_output