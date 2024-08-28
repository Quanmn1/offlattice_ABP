#!/bin/bash
{
name_exe="abp_pfaps_harmonic_qsaps_exp_slab"

# gcc ABP.c -o $name_exe -lm -O3 -Wall
name_all=$1
dt=0.0001
Lx=40
Ly=20
v=5
Dr=1
lp=$(echo "scale=0; $v / $Dr"  | bc)
lprf=$2
rmax_pfap=$(printf "%.2f" "$(echo "scale=4; $lp / $lprf" | bc)")
rho_large=$4
rho_small=$3
liquid_fraction=0.5
rho_m=25
# N=8000
# rmax_pfap=$1
lambda=1
phi=10
rmax_qsap=1
# epsilon=1
# epsilon=0.
epsilon=$(echo "scale=1; 100 * $rmax_pfap"  | bc)
epsilon=$(echo "scale=1; 100 * $rmax_pfap"  | bc)
# rho_small=20
# rho_large=$(echo "scale=1; 0.6 / $rmax_pfap / $rmax_pfap"  | bc)
# N=17500
# rho0=$1
# N=$(echo "scale=0; $rho0 * $Ly * $Lx"  | bc)
# rho_rf2=0.45
# N=$(echo "scale=0; $rho_rf2 * $Ly * $Lx / $rmax_pfap / $rmax_pfap"  | bc)
final_time=1000
density_box_size=2
ratio=$(echo "scale=1; $Ly / $Lx"  | bc)
timestep=10
data_store=$timestep
update_histo=2
histo_store=$timestep
start_time=0
resume="no"
terminal_x=2000
terminal_y=1000

name="$name_all"_"$rmax_pfap"
{
    # time ./$name_exe $dt $N $Lx $Ly $rho_m $v $lambda $phi $rmax_qsap $epsilon $rmax_pfap $Dr $final_time \
    # $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 1234
    time ./$name_exe $dt $rho_small $rho_large $liquid_fraction $Lx $Ly $rho_m $v $lambda $phi $rmax_qsap $epsilon $rmax_pfap $Dr $final_time \
    $density_box_size $start_time $update_histo $start_time $histo_store $start_time $data_store $name $resume 1234
} \
2>> "$name"_param

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

last=$(awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0  >> file_out}
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
    set palette defined ( 0 "orange", 1 "dark-orange" )

    # set limits to x and y axes
    set xr[0:$Lx]
    set yr[0:$Ly]
    set size ratio $ratio
    set terminal png size $terminal_x,$terminal_y
    set output "$i.png"
    unset key
    size=$(echo "scale=3; $rmax_pfap / 2"  | bc)
    # size=1
    set style fill border
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
# ffmpeg makes nicer movies

ffmpeg -loglevel fatal -y -r 10 -i "$dir"/data%0"$pad"d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p "$name"-gnuplot-density.mp4 

# Delete every file except the last one, for continuing simulation
last_file=$(printf "%s/data%0${pad}d" "$dir" "$last")
new_file="$dir"/"$name"_last_state
mv "$last_file" "$new_file"
rm "$dir"/data*

# in HOMOGENEOUS simulations, prepare a file for measurements of bulk density
file_out="$name"_density_data
rm $file_out
awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$name"_density

# threshold in pressure videos
max_sigma=600

# in HOMOGENEOUS simulations, prepare a file for measurements of bulk density
file_out="$name"_density_data
rm $file_out
awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$name"_density

# threshold in pressure videos
max_sigma=600

dir="$name"_sigmaIK
rm "$dir"/data*

for file in "$dir"/"$dir"_??;
do
# in HOMOGENEOUS simulations, prepare a file for measurements of sigma_IK
# file_out="$file"_data
# rm $file_out
# awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$file"

# plot sigmaIK
awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0 >> file_out}' "$file"

sigma=$(awk -v max_sigma="$max_sigma" 'BEGIN{max=0} 
sigma=$(awk -v max_sigma="$max_sigma" 'BEGIN{max=0} 
NF>2 {for (i = 1; i <= NF; i++) {
    if ($i > max) {
        max = $i
    }
}
} END{ {if (max > max_sigma) {print max_sigma} else {print max} } }' $file)

for i in "$dir"/data*
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    # echo "sigma $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$sigma]
    set palette defined ( 0 "white", 1 "dark-orange" )

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
for file in "$name"_sigmaAxx "$name"_sigma "$name"_Qxx;
do
echo "$file"
# in HOMOGENEOUS simulations, prepare a file for measurements of sigma
# file_out="$file"_data
# rm $file_out
# awk 'NF>2 {print $0  >> "'"$file_out"'"}' "$file"

awk 'BEGIN{iread=1;i=0;t=0;t_increment='"$timestep"';eps=0.000001;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i)}
NF==1  {if(iread==1) {i+=1;iread=0;t+=t_increment;file_out=sprintf("'"$dir"/'data%0'"$pad"'d",i);print $1 >> file_out}}
NF>2 {iread=1;print $0 >> file_out}' "$file"

sigma=$(awk -v max_sigma="$max_sigma" 'BEGIN{max=0} 
NF>2 {for (i = 1; i <= NF; i++) {
    if ($i > max) {
        max = $i
    }
}
} END{ {if (max > max_sigma) {print max_sigma} else {print max} } }' $file)

for i in "$dir"/data*
do
    # Read time from the first column of the first line
    Time=$(head -n 1 $i | awk '{printf("%04d",int($1))}')
    # echo "sigma $i $Time";
    gnuplot <<EOF
    set title 'Time $Time'
    set cbrange [0:$sigma]
    set palette defined ( 0 "white", 1 "dark-orange" )

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

exit
}