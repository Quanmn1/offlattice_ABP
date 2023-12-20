
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

Calculate density:
For each pair of points, in periodic boundary conditions: since the kernel is much smaller than the system size, take the min of the possible distances. That is the only distance that can contribute significantly to the density. DONE
Calculate density function on a grid of points, write the matrix to a file with parameters of the grid. Python code read in the file and plot using color.
- Check steady state density with velocity landscape. Looks good qualitatively!
- Starting at 0: Fit density profile to a gaussian and plot the parameters as a function of time. Good!

Merged two function files, scrap the noninteracting one? Can't do that because the attributes in parameters are required to be defined, even if the function isn't used. Maybe make a common function file, and import from it?

Handle Keyboard interrupt so that save data