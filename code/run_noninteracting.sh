gcc ABP_QSAPs_noninteracting.c -o abp_qsaps_noninteracting -lm
# dt N Lx Ly v KernelName KernelWidth FileName FinalTime NextStoreTime StoreTimeInterval Dr seed
name="non_interacting_test0"
./abp_qsaps_noninteracting 0.01 500 5 5 2.5 exp 0.2 $name 1 0 0.5 0.1 10
python3 noninteracting_visualize.py $name
ffmpeg -i ${name}_video/${name}_%d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p ${name}_video.mp4 