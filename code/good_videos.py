import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Parameters (You need to define or input these values accordingly)
rho = 200  # replace with your rho value
Lx, Ly = 20, 20  # replace with your x and y axis limits
ratio = 1  # set to the aspect ratio
terminal_x, terminal_y = 1500, 1500  # replace with your terminal size
rmax = 0.08  # replace with your rmax value
size = rmax / 2  # calculate size

lst = [f"pfaps_qsaps_test16_harmonic_largescaleeps_0.08_video/data{i:04}" for i in range(1,102)]

for (i, data_file) in enumerate(lst):
    time = 10*i
    # Load data from the file
    # Assuming data is in a space-separated text file with at least 4 columns
    data = np.loadtxt(data_file, skiprows=1)  # skip the first line which contains time

    # Prepare the figure
    plt.figure(figsize=(terminal_x / 100, terminal_y / 100))
    plt.title(f'Time {time}')

    # Set the axes limits
    plt.xlim(0, Lx)
    plt.ylim(0, Ly)
    plt.gca().set_aspect(ratio)

    # Define the color range and palette (orange to dark orange)
    norm = mcolors.Normalize(vmin=0, vmax=rho)
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["orange", "darkorange"])


    # Extract columns from the data
    x = data[:, 0]
    y = data[:, 1]
    colors = data[:, 3]

    # Create the plot with filled circles
    scatter = plt.scatter(x, y, s=size * 100, c=colors, cmap=cmap, norm=norm)

    # Add a colorbar
    plt.colorbar(scatter)

    # Save the output
    plt.savefig(f"{data_file}.png", dpi=100)
    plt.close()
