#!/bin/bash

name_all="pfaps_test11_diagram"
dt=0.05
# N=36000
liquid_fraction=0.5
Lx=400
Ly=200
v=0.5
epsilon=0.125
Dr=$1
# v_min=5
# v_max=$1
# rho_m=10
# rho_large=$(echo "scale=1; (($v_max/$v_min-2)/10*1.5+1.5)*$rho_m"  | bc)
# rho_small=$(echo "scale=1; (1-($v_max/$v_min-2)/10)*$rho_m"  | bc)
rho_large=$(echo "scale=2; 1.21 - 2.5 * $Dr"  | bc)
rho_small=$(echo "scale=2; 15.5 * $Dr - 0.2"  | bc)
echo "rho_large = $rho_large"
echo "rho_small = $rho_small"
epsilon=1
final_time=5000
density_box_size=5
rmax=1
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=20
data_store=20
update_histo=1
histo_store=20
start_time=1000

name="$name_all"_"$Dr"

# name="pfaps_test_16"
file="$name"_data
dir="$name"_video

# Number of times
M=$(awk 'NF==1 {m++} END{print m}' $file)
pad=$(echo ${#M} | awk '{print $1+1}')

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $4>rho {rho=$4} END{print rho}' $file)

# i is the marker of the file

# t is a time of the next block to read. eps is a small increment used to "compare" the times.

# iread is turned off when a new file is made and turned to 1 when data are stored in this file

# NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4,$6 >> file}}
# -> if it is a data line (NF>2), then if $1 is larger than the current time ($1>t-eps) then you should record

# NF==1 {if(iread==1) {i+=1;iread=0;file=sprintf("data%0'"$pad"'d",i)}}
# if you are between two blocks and you have just read the last block
# (iread==1) then you should increment i and open a new file, and
# increment t up to the next time you want to record

last=$(awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file}}
NF>2 {iread=1;print $1,$2,$3,$4  >> file}
END {printf("%d", i)}
' "$file")

# Delete every file except the last one, for continuing simulation
last_file=$(printf "%s/data%0${pad}d" "$dir" "$last")
new_file="$dir"/last_state
mv "$last_file" "$new_file"
rm "$dir"/data*
