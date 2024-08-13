#!/bin/bash

name_exe="abp_pfaps_harmonic_wall"

# gcc ABP.c -o $name_exe -lm -O3 -Wall

name_all="pfaps_test48_wca_wall"
dt=0.001
# N=34000
Lx=200
Ly=200
rmax=1
right_limit=$(echo "scale=6; $Lx - $rmax"  | bc)
rho0=$1
N=$(echo "scale=0; $rho0 * $Ly * $Lx"  | bc)
# N=$(echo "scale=0; 0.8 / $rmax / $rmax * $Ly * $Lx"  | bc)
v=5
epsilon=100
# epsilon=$(echo "scale=0; 50 * $rmax"  | bc)
# Pe=$1
# Dr=$(printf %.3f $(echo "scale=4; $v / $Pe / 0.89 + 0.0002" | bc)) # 0.0002 is to make it round up
# rf=$1
# Pe=$(echo "scale=4; 5 / $rf"  | bc)
# Dr=$(echo "scale=4; $v / $Pe" | bc)
Dr=0.5
# v_min=5
# v_max=$1
# rho_m=10
# rho_large=$(echo "scale=1; (($v_max/$v_min-2)/10*1.5+1.5)*$rho_m"  | bc)
# rho_small=$(echo "scale=1; (1-($v_max/$v_min-2)/10)*$rho_m"  | bc)
# rho_large=$(echo "scale=3; 1 + 0.05 / $Dr"  | bc)
# rho_large=1.6
# rho_small=$(echo "scale=3; 1 - 0.10 / $Dr"  | bc)
# rho_small=0.2
# echo "rho_large = $rho_large"
# echo "rho_small = $rho_small"
# rho_0=0.9
# rho_large=$3
# rho_small=$2
# liquid_fraction=0.5
# liquid_fraction=$(echo "scale=3; ( $rho_0-$rho_small ) / ( $rho_large-$rho_small )" | bc)
final_time=1000
# linear wall force
wall_size=10
right_wall=$(echo "scale=0; $Lx - $wall_size"  | bc)
omega=20
density_box_size=5
# rho_rf2=0.4
# N=$(echo "scale=0; $rho_rf2 / $rmax / $rmax * $Ly * $Lx"  | bc)
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=10
data_store=$timestep
update_histo=2
histo_store=$timestep
start_time=0
resume="no"
terminal_x=1500
terminal_y=1500

name="$name_all"_"$rho0"
# {
#     # time ./$name_exe $dt $rho_small $rho_large $liquid_fraction $Lx $Ly $v $epsilon $rmax $Dr $final_time \
#     # $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 1234
#     time ./$name_exe $dt $N $Lx $Ly $v $epsilon $rmax $Dr $final_time $wall_size $omega \
#     $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 1234
# } \
# 2>> "$name"_param

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

# in HOMOGENEOUS simulations, prepare a file for measurements of v(rho)
file_out="$name"_v
rm $file_out
awk 'NF>2 {print $4,$5  >> "'"$file_out"'"}' "$file"

#Calcul of max density
rho=$(awk 'BEGIN{rho=0} $4>rho {rho=$4} END{print rho}' $file)
echo "rho max is $rho"

#Check if the walls are strong enough to prevent particles at 2 ends interacting with each other
{
echo "Times (if any) at which particles from 2 ends interact through walls (not supposed to happen):"
awk 'BEGIN{t=0;}
NF==1 {t+=10;}
NF>2 {if ($1<'"$rmax"' || $1>'"$right_limit"') print t}' $file
} >> "$name"_param

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

echo "simulation"

for i in "$dir"/data*;
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    # echo "file $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0.7:$rho]

    # set limits to x and y axes
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size $terminal_x,$terminal_y
    set output "$i.png"
    unset key
    size=$(echo "scale=2; $rmax / 2"  | bc)
    set style fill solid
    set term png font ",25"
    set object 1 rect from  0,0 to $wall_size,$Ly fc rgb "yellow" fillstyle solid 0.5
    set object 1 rect from  $right_wall,0 to $Lx,$Ly fc rgb "yellow" fillstyle solid 0.5

    # x:y:size:color
    # skip the first line which contains time
    # pt 7 gives you a filled circle and ps 10 is the size, lt -1 solid line
    # us 1:2 w p pt 7 ps .2 lt -1, "$i" 
    plot "$i" skip 1 using 1:2:(size):4 with circles palette
    
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
# ffmpeg makes nicer movies

ffmpeg -loglevel fatal -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$name"-gnuplot-density.mp4 

# Delete every file except the last one, for continuing simulation
last_file=$(printf "%s/data%0${pad}d" "$dir" "$last")
new_file="$dir"/"$name"_last_state
mv "$last_file" "$new_file"
rm "$dir"/data*

# in HOMOGENEOUS simulations, prepare a file for measurements of sigma_IK
file_out="$name"_density_data
rm $file_out
awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$name"_density

# sigmaIK
dir="$name"_sigmaIK

for file in "$dir"/"$dir"_??;
do
# in HOMOGENEOUS simulations, prepare a file for measurements of sigma_IK
file_out="$file"_data
rm $file_out
awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$file"

#Calcul of max sigma
sigma=$(awk 'BEGIN{max=0} 
NF>2 {for (i = 1; i <= NF; i++) {
    if ($i > max) {
        max = $i
    }
}
} END{print max}' $file)

# plot sigma
awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0 >> file_out}' "$file"

echo "sigmaIK"

for i in "$dir"/data*
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    # echo "sigma $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$sigma]

    # set limits to x and y axes
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size $terminal_x,$terminal_y
    set output "$i.png"
    unset key
    set term png font ",25"

    pl "$i" skip 1 matrix with image
    
    set output
EOF
done

ffmpeg -loglevel fatal -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$file".mp4 
rm "$dir"/data*
done

# vid of nematic, sigmaA, and total sigma (sigmaA + sigmaIK)
for file in "$name"_sigmaA "$name"_sigma;
do
# in HOMOGENEOUS simulations, prepare a file for measurements of sigma
file_out="$file"_data
rm $file_out
awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$file"

awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0 >> file_out}' "$file"

sigma=$(awk 'BEGIN{max=0} 
NF>2 {for (i = 1; i <= NF; i++) {
    if ($i > max) {
        max = $i
    }
}
} END{print max}' $file)

echo $file

for i in "$dir"/data*
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    # echo "sigma $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$sigma]

    # set limits to x and y axes
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size $terminal_x,$terminal_y
    set output "$i.png"
    unset key
    set term png font ",25"

    pl "$i" skip 1 matrix with image
    
    set output
EOF
done

ffmpeg -loglevel fatal -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$file".mp4 
rm "$dir"/data*
done