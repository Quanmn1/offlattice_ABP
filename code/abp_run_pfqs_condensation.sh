#!/bin/bash
{
# start=5
# ending=60
# space=5
# start=0.0
# ending=1.0
# space=0.2
# start=0.04
# ending=0.14
# space=0.02
name_all="pfaps_qsaps_condensed_consistent4"
mode="pfqs"
num_segments=5
init="homo"
# input_density="pfaps_qsaps_absorbing_hard"
homo_fit="gauss"
slab_fit="tanh"
# vars=$(seq $start $space $ending)
# num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars="0.32 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58"
# vars="0.19 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50"
# vars="0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.19 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50"
# vars="0.19 0.20"
# vars="0.0000 0.0200 0.0400 0.0600 0.0800 0.1000 0.1200 0.1425 0.1475 0.1500"
vars="0.2000"
# vars="0.000 0.020 0.040 0.060 0.080 0.100 0.120 0.140 0.240 0.260 0.280 0.300 0.320 0.340 0.360 0.380 0.400 0.420 0.440 0.460 0.480 0.500"
# vars="0.1600 0.1800 0.2000 0.2200 0.2400 0.2600 0.2800"
# vars="0.2800 0.3000 0.3200 0.3400 0.3600 0.3800 0.4000 0.4200 0.4400 0.4600 0.4800 0.5000"
# vars="0.1600 0.1800 0.2000 0.2200 0.2400 0.2600 0.2800 0.3000 0.3200 0.3400 0.3600 0.3800 0.4000 0.4200 0.4400 0.4600 0.4800 0.5000"
# vars="0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.19 0.20"
# vars="0.150"
pad=4
num=$(echo $vars | wc -w)

{ 
time {
# tail -n +1 "$input_density" | awk '{print $1, $2, $3}' | xargs -P$num -n 3 ./abp_pfqs_condensation.sh $name_all
printf "$vars" | xargs -P$num -d ' ' -I{} ./abp_pfqs_condensation.sh $name_all {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfqs_condensation.sh {}
# python3 analyze_histogram_coarsened.py $name_all "$vars" 2 4
# python3 analyze_phases.py $name_all $mode "$vars" $num_segments $init $pad $homo_fit $slab_fit
# python3 analyze_histogram_randomizedgrid.py $name_all "0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.1425 0.145" $pad no
# python3 analyze_histogram_randomizedgrid.py $name_all "0.1475 0.1500" $pad yes
# seq $start $space $ending | xargs -P6 -I{} ./get_last.sh {}
# python3 analyze_veff.py $name_all $mode "$vars" $num_segments
# python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaA"
# python3 analyze_sigma.py $name_all $mode "$vars" $num_segments 0 "sigmaIK"
# python3 plot_guides.py $name_all pfaps_test24_largePe
} \
} 2> "$name_all"_output

exit
}