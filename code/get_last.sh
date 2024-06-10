#!/bin/bash

name_all="pfaps_qsaps_test4"
dt=0.01
N=17500
Lx=30
Ly=30
rho_m=25
v=0.5
lambda=1
phi=10
rmax_qsap=1
epsilon=0.125
rmax_pfap=$1
Dr=0.1
final_time=1000
density_box_size=2
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=50
data_store=$timestep
update_histo=5
histo_store=$timestep
start_time=0
resume="no"

name="$name_all"_"$rmax_pfap"

# name="pfaps_test_16"
file="$name"_data
dir="$name"_video

if [ ! -d "$dir" ]; then
    mkdir "$dir"
fi

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
NF>2 {iread=1;print $1,$2,$3,$4,$5  >> file}
END {printf("%d", i)}
' "$file")

# Delete every file except the last one, for continuing simulation
last_file=$(printf "%s/data%0${pad}d" "$dir" "$last")
new_file="$dir"/last_state
mv "$last_file" "$new_file"
rm "$dir"/data*