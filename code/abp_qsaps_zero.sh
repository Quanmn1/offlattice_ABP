#!/bin/bash

name_exe=abp_qsaps_zero_linear
gcc ABP.c -o $name_exe -lm -O3 -Wall

name_all="qsaps_condensation_test2_stop"
dt=0.02
N=20000
# liquid_fraction=0.5
Lx=40
Ly=40
lambda=1
# rho0=$(echo "scale=1; $N / $Lx / $Ly"  | bc)
# rho_large=$(echo "scale=1; 52 + 50 * $lambda"  | bc)
# rho_small=$(echo "scale=1; 46 - 30 * $lambda"  | bc)
# echo "rho_large = $rho_large"
# echo "rho_small = $rho_small"
rho_m=25
v=5
phi=1
rmax=1
Dr=1
final_time=200
density_box_size=4
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=1
data_store=$timestep
update_histo=1
histo_store=$timestep
start_time=0
resume="no"

name="$name_all"
{
    time ./$name_exe $dt $N $Lx $Ly $rho_m $v $lambda $phi $rmax $Dr $final_time \
    $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 1234
# ./$name_exe $dt $rho_small $rho_large $liquid_fraction $Lx $Ly $rho_m $v $lambda $phi $rmax $Dr $final_time \
#     $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 10
} \
2>> "$name"_param

file="$name"_data
dir="$name"_video

if [ ! -d "$dir" ]; then
    mkdir "$dir"
else
    rm "$dir"/data*
fi

# Number of times
M=$(awk 'NF==1 {m++} END{print m}' $file)
pad=$(echo ${#M} | awk '{print $1+1}')

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $4>rho {rho=$4} END{print rho}' $file)

echo "rho max is $rho"

# in HOMOGENEOUS simulations, prepare a file for measurements of v(rho)
# file_out="$name"_v
# rm $file_out
# awk 'NF>2 {print $4,$5  >> "'"$file_out"'"}' "$file"

# make files for each time snapshots for plotting
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
NF>2 {iread=1;print $0  >> file}
END {printf("%d", i)}
' "$file")

for i in "$dir"/data*;
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    echo "file $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$rho]

    # set limits to x and y axes
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size 1600,1000
    set output "$i.png"
    unset key
    size=$(echo "scale=1; $rmax / 10" | bc)
    set style fill solid
    set term png font ",25"

    # x:y:size:color
    # skip the first line which contains time
    # pt 7 gives you a filled circle and ps 10 is the size, lt -1 solid line
    # us 1:2 w p pt 7 ps .2 lt -1, "$i" 
    pl "$i" skip 1 us 1:2:(size):4 w circles palette
    
    set output
EOF
done

# http://guide.debianizzati.org/index.php/FFmpeg 
# https://en.wikibooks.org/wiki/FFMPEG_An_Intermediate_Guide/image_sequence
# http://www.ffmpeg.org/ffmpeg.html
# Syntax ffmpeg { comandi  inputfile} {comandi outputfile}
# -r 10 is the framerate  -i is for the input  file
# data%0'$pad'd.png is to ordinate images with leading zeros
# before % and after d is taaken literaly while $pad is thenumber or file that will be taken
# -y overwrite output file without asking
# libx264 delivers better quality

# ffmpeg -r 10 -i "$dir"/data%0"$pad"d.png -b:a 16M -vcodec libx264 "$name"-gnuplot-density.mp4 

ffmpeg -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$name"-gnuplot-density.mp4 

# ffmpeg makes nicer movies

# Delete every file except the last one, for continuing simulation
last_file=$(printf "%s/data%0${pad}d" "$dir" "$last")
new_file="$dir"/"$name"_last_state
mv "$last_file" "$new_file"
rm "$dir"/data*