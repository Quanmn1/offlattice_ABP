#!/bin/bash
{
start=50
ending=55
space=5
# start=0.15
# ending=0.23
# space=0.01
# start=10
# ending=160
# space=10
name_all="pfaps_qsaps_test26_r16_noncrit"
mode="pfqs"
num_segments=100
init="homo"
input_density="pfaps_qsaps_blob5_init_densities"
pad=0
homo_fit="gauss"
slab_fit="tanh"
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="9"
# vars="0.040 0.050 0.055 0.060 0.065 0.070 0.075 0.080"
# vars="0.065 0.070 0.075 0.080 0.085 0.090 0.095 0.100 0.110 0.120 0.130 0.140 0.150 0.160 0.170 0.180 0.190 0.200 0.210 0.220 0.230 0.240"
# vars="0.160 0.170 0.180 0.190 0.200 0.210 0.220 0.230 0.240"
# vars="0.060 0.065"
# vars="0.040 0.050 0.055 0.060 0.065"
vars=$(seq $start $space $ending)
num=$(echo $vars | wc -w)
{ 
time {
# tail -n +1: show every line in the file
# tail -n +1 "$input_density" | awk '{print $1, $2, $3, NR}' | xargs -P$num -n 4 ./abp_pfaps_qsaps.sh $name_all
# printf "$vars" | xargs -P$num -d ' ' -I {} ./abp_pfaps_qsaps.sh $name_all {}
seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps_qsaps.sh $name_all {}
# python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $pad $homo_fit $slab_fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
python3 analyze_veff.py $name_all $mode "$vars" $num_segments
python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaA"
python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaIK"
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} 2> "$name_all"_output

exit
}