import numpy as np
from matplotlib import pyplot as plt
import os
import sys

def plot_arrow(ax, x, y, theta, scale=0.3):
    dx = np.cos(theta) * scale
    dy = np.sin(theta) * scale
    ax.arrow(x-dx/2, y-dy/2, dx, dy, width=0.02, color='C0')

def plot_point(ax, x, y):
    ax.plot(x, y, marker='.', ms=1, ls='')

name = sys.argv[1]
density_file = name + '_density'
param_file = name + '_param'

# Read parameters into dict params, with string values
with open(param_file, 'r') as f:
    params_info = f.readline().split()
    params_values = f.readline().split()
    params = {var: value for (var, value) in zip(params_info[2:], params_values[1:])}

N = int(params['N'])
v = params['v']
Dr = params['Dr']
Lx = float(params['Lx'])
Ly = float(params['Ly'])
next_store_time = float(params['NextStoreTime'])
store_time_interval = float(params['StoreTimeInterval'])
final_time = float(params['FinalTime'])

density_folder = name+'_density_folder'
if not os.path.exists(density_folder):
    os.mkdir(density_folder)

with open(density_file, 'r') as f:
    density_param = f.readline().split()

half_number_of_points_x = int(density_param[0])
half_number_of_points_y = int(density_param[1])
density_grid_spacing = float(density_param[2])

# One time
data = np.loadtxt(density_file, skiprows=1)
number_of_times = data.shape[0] // (1+2*half_number_of_points_x)
data = data.reshape(number_of_times, 1+2*half_number_of_points_x, 1+2*half_number_of_points_y)
times = np.arange(next_store_time, final_time+0.01, store_time_interval) # sus
assert number_of_times == len(times), str(number_of_times) + ' ' + str(len(times))

os.chdir(density_folder)
for i in range(data.shape[0]):
    t = times[i]
    fig, ax = plt.subplots(figsize = (6,6))
    ax.set_title(f'$N={N}, v={v}, Dr={Dr}, t={t:.2f}$')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(-Lx/2,Lx/2)
    ax.set_ylim(-Ly/2,Ly/2)
    cmap = plt.cm.get_cmap('viridis')
    cax = ax.imshow(data[i], cmap=cmap, origin='lower', interpolation='nearest', extent=(-Lx/2,Lx/2,-Ly/2,Ly/2))
    cbar = fig.colorbar(cax)
    cbar.set_label('Density')
    plt.savefig(name + f'_{i}.png', dpi=300, bbox_inches='tight')
    plt.close()