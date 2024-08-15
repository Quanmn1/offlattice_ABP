start=0.025
ending=0.200
space=0.020
name_all="pfaps_test46_harmonic_yongfeng"
num=$(echo "scale=0; $ending/$space - $start/$space + 1" | bc)
# vars=$(seq $start $space $ending)
# vars="0.1 0.3 0.5 0.7 0.9 1.0 1.1 1.2 1.3"
vars="0.01 0.02 0.05 0.08 0.10 0.15"
num=$(echo "$vars" | wc -w)

# tail -n +2 "$input_density" | xargs -P9 -n 3 ./abp_pfaps.sh
printf "$vars" | xargs -P$num -d ' ' -I{} ./get_last.sh {}
# seq $start $space $ending | xargs -P$num -I{} ./abp_pfaps.sh {}

