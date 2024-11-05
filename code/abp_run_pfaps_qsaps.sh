#!/bin/bash
{
# start=5
# ending=60
# space=5
# start=0.10
# ending=0.40
# space=0.02
# start=10
# ending=180
# space=10
name_all="pfaps_qsaps_test17_metastable_0.05"
mode="pfqs"
num_segments=10
init="homo"
# input_density="pfaps_qsaps_test16_harmonic_largescaleeps_histo_phase_diagram"
homo_fit="gauss"
slab_fit="average"
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.14 0.16 0.18 0.20 0.24 0.28 0.32"
# vars="5 10 15 20 25 30 40 50 100 200 300 400 500 600"
vars="0.05"
# vars=$(seq $start $space $ending)
num=$(echo $vars | wc -w)
{ 
time {
# tail -n +2 "$input_density" | awk '{print $1, $2, $3}' | xargs -P$num -n 3 ./abp_pfaps_qsaps.sh $name_all
printf "$vars" | xargs -P$num -d ' ' -I {} ./abp_pfaps_qsaps.sh $name_all {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps_qsaps.sh $name_all {}
python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit $slab_fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode "$vars" $num_segments
# python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaA"
# python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaIK"
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} 2> "$name_all"_output

exit
}