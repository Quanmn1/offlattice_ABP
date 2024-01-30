#!/bin/bash

## $1 is the filename (*-data, with the last column corresponding to the local density)
## $2 is Lx
## $3 is Ly
## $4 is the particle size
## $5 is time increment

name="pfaps_test6"
dt=0.01
N=200
Lx=20
Ly=20
v=10
epsilon=0.5
Dr=1
final_time=100
density_box_size=2
rmax=1
ratio=1
size=0.1
xlim=

# ./abp_pfaps_density $dt $N $Lx $Ly $v $epsilon $Dr 1 $final_time $density_box_size 0 1 0 1 0 1 $name 10

file="$name"_data
dir="$name"_video

if [ ! -d "$dir" ]; then
    # If not, create it
    mkdir "$dir"
fi
rm "$dir"/data*

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

# sus about $1>t-eps. check again

awk 'BEGIN{iread=0;i=0;t=0;t_increment='"$dt"';eps=0.000001;file=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file}}
NF>2 {if($1>t-eps) {iread=1;print $1,$2,$3,$4  >> file}}
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
    set xr[$(echo "scale=1; -$2/2"  | bc):$(echo "scale=1; -$2/2"  | bc)] 
    set yr[$(echo "scale=1; -$3/2"  | bc):$(echo "scale=1; -$3/2"  | bc)]
    set size ratio $ratio
    set terminal png size 1600,1600
    set output "$i.png"
    unset key
    size=$size
    set style fill solid
    set term png font ",25"

    # x:y:size:color
    # skip the first line which contains time
    # pt 7 gives you a filled circle and ps 10 is the size, lt -1 solid line
    pl "$i" skip 1 us 1:2 w p pt 7 ps .2 lt -1, "$i" us 1:2:(size):4 w circles lc palette 
    
    set term png font ",25" 
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

ffmpeg -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$name"-gnuplot-density.mp4 

# ffmpeg makes nicer movies

rm "$dir"/data*

