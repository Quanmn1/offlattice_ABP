import numpy as np
from matplotlib import pyplot as plt
import os
import sys

def plot_arrow(ax, x, y, theta, scale=0.3):
    dx = np.cos(theta) * scale
    dy = np.sin(theta) * scale
    ax.arrow(x-dx/2, y-dy/2, dx, dy, width=0.02, color='C0')

def plot_point(ax, x, y):
    ax.plot(x, y, marker='.', ms=1)

name = sys.argv[1]
param_file = name + '_param'
data_file = name + '_data'

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

# create a folder in which the images are stored
video_folder = name+'_video'
if not os.path.exists(video_folder):
    os.mkdir(video_folder)

with open(data_file, 'r') as f:
    os.chdir(video_folder)
    for line in f:
        # One time
        data = line.split()
        t = float(data[0])
        fig, ax = plt.subplots(figsize = (Lx, Ly))
        ax.set_title(f'$N={N}, v={v}, Dr={Dr}, t={t:.2f}$')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(-Lx/2,Lx/2)
        ax.set_ylim(-Ly/2,Ly/2)
        for i in range(0, len(data)//5):
            x = float(data[i*5+2])
            y = float(data[i*5+3])
            theta = float(data[i*5+4])
            plot_arrow(ax, x, y, theta)
            # plot_point(ax, x, y)
        assert i==N-1, str(i) + ' ' + str(N) # sanity check
        plt.savefig(name + f't={t:.2f}.png', dpi=300, bbox_inches='tight')
        plt.close()
