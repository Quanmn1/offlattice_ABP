#!/bin/bash

# syntax: ./good_videos.sh test_name variable

{
name_all=$1
Lx=20
Ly=20
v=5
Dr=1
rho_m=25
rmax_pfap=$2
lambda=1
phi=10
rmax_qsap=1
epsilon=$(echo "scale=1; 100 * $rmax_pfap"  | bc)
final_time=1000
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=10
terminal_x=1500
terminal_y=1500

name="$name_all"_"$rmax_pfap"

file="$name"_data
dir="$name"_video

if [ ! -d "$dir" ]; then
    mkdir "$dir"
else
    rm "$dir"/data*
    rm "$dir"/histogram*
fi

# Number of times
M=$(awk 'NF==1 {m++} END{print m}' $file)
pad=$(echo ${#M} | awk '{print $1+1}')

# Calcul of max density
# rho=$(awk 'BEGIN{rho=0} $4>rho {rho=$4} END{print rho}' $file)
rho=90
echo "rho max is $rho"

last=$(awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0  >> file_out}
END {printf("%d", i)}
' "$file")

python3 good_videos.py $name $Lx $Ly $rmax_pfap $rho $terminal_x $terminal_y $timestep $last

# ffmpeg -loglevel fatal -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$name"-python.mp4 

last_file=$(printf "%s/data%0${pad}d.png" "$dir" "$last")
new_file="$dir"/"$name"_last_state.png
mv "$last_file" "$new_file"
rm "$dir"/data*

exit
}