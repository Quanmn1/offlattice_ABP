gcc ABP_QSAPs.c -o abp_qsaps_hashing -lm
# dt N Lx Ly box_size rho_m v_min v_max interaction_range FileName FinalTime NextStoreTime StoreTimeInterval Dr seed
name="qsaps_coarsening1"
./abp_qsaps_hashing 0.05 10000 20 20 1 20 0.5 3 1 $name 1000 0 10 1 10
python3 visualize.py $name
ffmpeg -i ${name}_video/${name}_%d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p ${name}_video.mp4 