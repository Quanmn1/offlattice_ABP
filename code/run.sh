gcc ABP_QSAPs.c -o abp_qsaps_interacting -lm
# dt N Lx Ly rho_m v_min v_max KernelName KernelWidth FileName FinalTime NextStoreTime StoreTimeInterval Dr seed
name="interacting_test_bash"
./abp_qsaps_interacting 0.01 500 5 5 20 0.1 2.5 exp 0.2 $name 1 0 0.5 0.1 10
python3 visualize.py $name
ffmpeg -i ${name}_video/${name}_%d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p ${name}_video.mp4 