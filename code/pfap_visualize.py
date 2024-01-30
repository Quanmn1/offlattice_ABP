import numpy as np
from matplotlib import pyplot as plt
import os
import sys
from matplotlib.ticker import FixedLocator

def plot_arrow(ax, x, y, theta, scale=0.3):
    dx = np.cos(theta) * scale
    dy = np.sin(theta) * scale
    ax.arrow(x-dx/2, y-dy/2, dx, dy, width=0.02, color='C0')

def plot_point(ax, x, y, label='', color='C0'):
    ax.plot(x, y, marker='.', ms=4, label=label, color=color)

def get_edges(array):
    distance = (array[1] - array[0])/2
    return np.append(array-distance, array[-1] + distance)

name = sys.argv[1]
param_file = name + '_param'
data_file = name + '_data'
histogram_file = name + '_histogram'

# Read parameters into dict params, with string values
with open(param_file, 'r') as f:
    params_info = f.readline().split()
    params_values = f.readline().split()
    params = {var: value for (var, value) in zip(params_info[2:], params_values[1:])}

N = int(params['N'])
v = params['v']
r0 = float(params['r_max'])
Dr = params['Dr']
Lx = float(params['Lx'])
Ly = float(params['Ly'])
density_box_size = float(params['density_box_size'])
number_of_boxes_x = int(Lx/r0/density_box_size)
number_of_boxes_y = int(Ly/r0/density_box_size)
rho = N/Lx/Ly
rho_max = 10*rho
number_of_boxes = number_of_boxes_x * number_of_boxes_y
max_number = N // number_of_boxes * 10

densities = np.zeros(max_number)
counts = np.zeros(max_number)

# create a folder in which the images are stored
video_folder = name+'_video'
if not os.path.isdir(video_folder):
    os.mkdir(video_folder)

histogram_folder = name+'_histogram_video'
if not os.path.isdir(histogram_folder):
    os.mkdir(histogram_folder)

with open(histogram_file, 'r') as f:
    os.chdir(histogram_folder)
    counter = 0
    for line in f:
        data = line.split()
        if len(data) == 1:
            # start of a time
            counter += 1
            row = 0
            t = float(data[0])

            fig, ax = plt.subplots(figsize = (6,6))
            ax.set_title(fr'$N={N}, v={v}, \rho={rho:.2f}, D_r={Dr}, t={t:.2f}$')
            ax.set_xlabel(r'$\rho$')
            ax.set_ylabel('Averaged count')
            ax.set_xlim(0, rho_max)
            ax.set_ylim(0, number_of_boxes)

        elif len(data) == 0:
            # end of a time
            plt.stairs(counts, get_edges(densities), fill=True)
            plt.savefig(name + f'_{counter}.png', dpi=300, bbox_inches='tight')
            plt.close()

        else:
            densities[row] = float(data[0])
            counts[row] = float(data[1])
            row += 1

os.chdir('..')   
with open(data_file, 'r') as f:
    os.chdir(video_folder)
    counter = 0
    for line in f:
        data = line.split()
        if len(data) == 1:
            # start of a time
            counter += 1
            t = float(data[0])

            fig, ax = plt.subplots(figsize = (6,6))
            ax.set_title(fr'$N={N}, v={v}, \rho={rho:.2f}, D_r={Dr}, t={t:.2f}$')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_xlim(-Lx/2,Lx/2)
            ax.set_ylim(-Ly/2,Ly/2)

            ax.grid(which= 'major')
            x_ticks = np.linspace(-Lx/2, Lx/2, number_of_boxes_x+1)
            y_ticks = np.linspace(-Ly/2, Ly/2, number_of_boxes_y+1)
            ax.xaxis.set_major_locator(FixedLocator(x_ticks))
            ax.yaxis.set_major_locator(FixedLocator(y_ticks))

        elif len(data) == 0:
            # end of a time
            # plt.colorbar(label="Density")
            plt.savefig(name + f'_{counter}.png', dpi=300, bbox_inches='tight')
            plt.close()

        else:
            x = float(data[0])
            y = float(data[1])
            theta = float(data[2])
            rho = float(data[3])
            plot_point(ax, x, y)
os.chdir('..')