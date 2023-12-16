
N bacteria that do run and tumble at some speed 
- define variables
- read parameters: store the input parameters
- make the bacteria evolve in time:
	- N bacteria evolve in parallel
	- If can solve explicitly for prob distribution: just sample
	- Nonlinear interaction: can't do exact -> time-stepping
	- Time-stepping: Euler, RK. Don't want to be stuck in one choice.
	- Keep track:
		- r(t)
		- $$\theta(t)$$
		- $$\rho[r(t)] \sim O(N^2)$$ bad
		- -> Spatial hasing
		- $$\rho(r_i)=\sum_j K(r_i-r_j)$$ with $$K$$ a window function of width $$w$$
		- Do that fast by creating boxes in space, for each particle only look at neighboring boxes and count particles in those boxes _> $$O(N)$$
	- pseudocode:
		- while (t<Tf) {
			- Move particles
			- t += dt
			- store data
		- }
 - N bacteria -> array of N particles

- Pointers:
int * titi
space not important
can read as titi is a pointer to an int, or * titi is an int
Can manipulate as array: * (titi+1) = titi[1]


ffmpeg -pattern_type glob -i "interacting_test_video/*.png" -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p interacting_text_video.mp4