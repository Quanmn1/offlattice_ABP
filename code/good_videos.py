import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import sys

plt.rcParams.update({
    "font.size": 9,
    "text.usetex": True,            # Use LaTeX for text
    "font.family": "serif",         # Set font to serif
    "font.serif": ["Computer Modern"]  # Use the Computer Modern font (LaTeX default)
})

# Parameters (You need to define or input these values accordingly)
filename = sys.argv[1]
arguments = [float(sys.argv[i]) for i in range(2, len(sys.argv))]
Lx, Ly, rmax, rho_max, terminal_x, terminal_y, timestep, last = arguments
ratio = Ly/Lx
size = rmax / 2  # calculate size

lst = [f"{filename}_video/data{i:04}" for i in range(1,int(last)+1)]

for (i, data_file) in enumerate(lst):
    if i < int(last) -1:
        continue
    time = timestep*i
    # Load data from the file
    # Assuming data is in a space-separated text file with at least 4 columns
    data = np.loadtxt(data_file, skiprows=1)  # skip the first line which contains time

    # Prepare the figure
    fig, ax = plt.subplots()
    ax.set_title(fr'$r_f={rmax}$')

    # Set the axes limits
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect(ratio)
    ax.tick_params(axis='both', which='both', length=0, labelbottom=False, labelleft=False)
    M = ax.transData.get_matrix()
    xscale = M[0,0]
    yscale = M[1,1]

    # Define the color range and palette (orange to dark orange)
    norm = mcolors.Normalize(vmin=0, vmax=rho_max)
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["blanchedalmond", "darkorange"])

    # Extract columns from the data
    x = data[:, 0]
    y = data[:, 1]
    colors = data[:, 3]

    # Create the plot with filled circles
    scatter = ax.scatter(x, y, s=size**2 * xscale * yscale, c=colors, cmap=cmap, norm=norm, facecolors='none')

    # Add a colorbar
    plt.colorbar(scatter, label=r"$\tilde{\rho}$")

    # Save the output
    plt.savefig(f"{data_file}.png", dpi=500, bbox_inches="tight")
    plt.close()
