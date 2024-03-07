#!/bin/bash

name_exe=abp_qsaps_density2
gcc ABP.c -o $name_exe -lm

## $1 is the filename (*-data, with the last column corresponding to the local density)
## $2 is Lx
## $3 is Ly
## $4 is the particle size
## $5 is time increment

name_all="qsaps_phase_diagram2"
dt=0.05
N=125000
Lx=50
Ly=50
rho_m=50
# v_min=5
# v_max=$1
v=5
lambda=$1
phi=20
Dr=1
rmax=1
final_time=5000
density_box_size=5
ratio=1
timestep=20
data_store=20
update_histo=1
histo_store=40

name="$name_all"_"$lambda"
./$name_exe $dt $N $Lx $Ly $rho_m $v $lambda $phi $Dr $rmax $final_time \
    $density_box_size 0 $update_histo 0 $histo_store 0 $data_store $name 10

file="$name"_data
dir="$name"_video

if [ ! -d "$dir" ]; then
    mkdir "$dir"
else
    rm "$dir"/data*
fi

# Number of times
M=$(awk 'NF==1 {m++} END{print m}' $file)
echo "M=$M"
pad=$(echo ${#M} | awk '{print $1+1}')
echo "pad $pad"

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $4>rho {rho=$4} END{print rho}' $file)

echo "rho max is $rho"

# i is the marker of the file

# t is a time of the next block to read. eps is a small increment used to "compare" the times.

# iread is turned off when a new file is made and turned to 1 when data are stored in this file

# NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4,$6 >> file}}
# -> if it is a data line (NF>2), then if $1 is larger than the current time ($1>t-eps) then you should record

# NF==1 {if(iread==1) {i+=1;iread=0;file=sprintf("data%0'"$pad"'d",i)}}
# if you are between two blocks and you have just read the last block
# (iread==1) then you should increment i and open a new file, and
# increment t up to the next time you want to record

awk 'BEGIN{iread=0;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file}}
NF>2 {iread=1;print $1,$2,$3,$4  >> file}
' $file

for i in "$dir"/data*;
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    echo "file $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$rho]

    # set limits to x and y axes
    # set xr[$(echo "scale=1; -$Lx/2"  | bc):$(echo "scale=1; $Lx/2"  | bc)] 
    # set yr[$(echo "scale=1; -$Ly/2"  | bc):$(echo "scale=1; $Ly/2"  | bc)]
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size 1600,1600
    set output "$i.png"
    unset key
    size=$rmax
    set style fill solid
    set term png font ",25"

    # x:y:size:color
    # skip the first line which contains time
    # pt 7 gives you a filled circle and ps 10 is the size, lt -1 solid line
    # us 1:2 w p pt 7 ps .2 lt -1, "$i" 
    pl "$i" skip 1 us 1:2:(size):4 w circles lc palette 
    
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

rm "$dir"/data*

# Plot histogram

file="$name"_histogram

awk 'BEGIN{iread=0;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file=sprintf("'"$dir"/'histogram%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file=sprintf("'"$dir"/'histogram%0'"$pad"'d",i)};print $1 >> file}
NF==2 {iread=1;print $1,$2  >> file}
' $file

for i in "$dir"/histogram*;
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    echo "file $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'

    # set limits to x and y axes
    # set xr[$(echo "scale=1; -$Lx/2"  | bc):$(echo "scale=1; $Lx/2"  | bc)] 
    # set yr[$(echo "scale=1; -$Ly/2"  | bc):$(echo "scale=1; $Ly/2"  | bc)]
    set size ratio $ratio
    set terminal png size 1600,1600
    set output "$i.png"
    unset key
    set style fill solid
    # set style data histograms
    set boxwidth 0.5 relative
    set term png font ",25"
    set xtics 0.2

    # x:y:size:color
    # skip the first line which contains time
    # pt 7 gives you a filled circle and ps 10 is the size, lt -1 solid line
    # pl "$i" skip 1 us 2:xticlabels(sprintf("%.2f",column(1)))
    pl "$i" skip 1 us 1:2 with boxes
    
    set term png font ",25" 
    set output
EOF
done

find "$dir" -type f ! -name "*.*" -exec rm -f {} +
