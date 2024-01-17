import numpy as np
from matplotlib import pyplot as plt
import os
import sys
from matplotlib.ticker import FixedLocator

def plot_arrow(ax, x, y, theta, scale=0.3):
    dx = np.cos(theta) * scale
    dy = np.sin(theta) * scale
    ax.arrow(x-dx/2, y-dy/2, dx, dy, width=0.02, color='C0')

def plot_point(ax, x, y, label=''):
    ax.plot(x, y, marker='.', ms=1, color='C0', label=label)

name = sys.argv[1]
param_file = name + '_param'
data_file = name + '_data'

# Read parameters into dict params, with string values
with open(param_file, 'r') as f:
    params_info = f.readline().split()
    params_values = f.readline().split()
    params = {var: value for (var, value) in zip(params_info[2:], params_values[1:])}

N = int(params['N'])
v_max = params['v_max']
v_min = params['v_min']
rho_m = params['rho_m']
r0 = params['interaction_range']
Dr = params['Dr']
Lx = float(params['Lx'])
Ly = float(params['Ly'])
NxBox = int(float(params['Lx'])/float(params['box_size']))
NyBox = int(float(params['Ly'])/float(params['box_size']))

# create a folder in which the images are stored
video_folder = name+'_video'
if not os.path.exists(video_folder):
    os.mkdir(video_folder)

with open(data_file, 'r') as f:
    os.chdir(video_folder)
    counter = 0
    for line in f:
        # One time
        data = line.split()
        t = float(data[0])
        fig, ax = plt.subplots(figsize = (6,6))
        ax.set_title(f'$N={N}, v_{{max}}={v_max}, v_{{min}}={v_min}, \\rho_m={rho_m}, r_0={r0}, Dr={Dr}, t={t:.2f}$')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(-Lx/2,Lx/2)
        ax.set_ylim(-Ly/2,Ly/2)
        # ax.grid(which= 'major')
        # x_ticks = np.linspace(-Lx/2, Lx/2, NxBox+1)
        # y_ticks = np.linspace(-Ly/2, Ly/2, NyBox+1)
        # ax.xaxis.set_major_locator(FixedLocator(x_ticks))
        # ax.yaxis.set_major_locator(FixedLocator(y_ticks))

        for i in range(0, len(data)//8):
            x = float(data[i*8+2])
            y = float(data[i*8+3])
            theta = float(data[i*8+4])
            rho = float(data[i*8+5])
            bi = int(data[i*8+6])
            bj = int(data[i*8+7])
            # plot_arrow(ax, x, y, theta)
            plot_point(ax, x, y, label=fr"$\rho$={rho:.2f}, bi={bi}, bj={bj}")
        assert i==N-1, str(i) + ' ' + str(N) # sanity check
        # plt.legend(loc=(1,0.3))
        plt.savefig(name + f'_{counter}.png', dpi=300, bbox_inches='tight')
        plt.close()
        counter += 1
