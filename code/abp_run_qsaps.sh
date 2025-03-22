#!/bin/bash
{
start=0.7
ending=1.3
space=0.1
name_all="qsaps_lambda1_phase_diagram"
mode="qsap"
num_segments=5
init="homo"
# input_density="pfaps_phase_diagram1_histo_phase_diagram"
homo_fit="gauss"
slab_fit="average"
vars=$(seq $start $space $ending)
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="0.18 0.20 0.22 0.24 0.26 0.28 0.32"
num=$(echo $vars | wc -w)

{ 
time {
# tail -n +2 "$input_density" | xargs -P11 -n 3 ./abp_qsaps_zero.sh
# printf "$vars" | xargs -P$num -d ' ' -I{} ./abp_qsaps_exp.sh $name_all {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_qsaps_exp.sh $name_all {}
python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $homo_fit $slab_fit
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode $start $ending $space
# python3 analyze_sigma.py $name_all $mode $start $ending $space
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} 2> "$name_all"_output

exit
}