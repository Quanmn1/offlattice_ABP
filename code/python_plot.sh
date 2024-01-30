gcc ABP.c -o abp_pfaps_density -lm
# ./abp_pfaps_density dt N Lx Ly v epsilon Dr r_max final_time density_box_size next_histogram_update histogram_update_interval next_histogram_store histogram_store_interval next_store_time store_time_interval file_name seed 
name="pfaps_test7"
dt=0.05
N=1000
Lx=50
Ly=50
v=0.5
epsilon=0.125
Dr=0.025
final_time=1000
density_box_size=2

{ time ./abp_pfaps_density $dt $N $Lx $Ly $v $epsilon $Dr 1 $final_time $density_box_size 0 1 0 1 0 1 $name 10 ; } 2>> time.txt
{ time python3 pfap_visualize.py $name ; } 2>> time.txt
ffmpeg -r 10 -i ${name}_video/${name}_%d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p ${name}_video.mp4 