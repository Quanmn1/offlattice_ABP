import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 14})
marker_size = 10

mode = sys.argv[1]

file = "PhaseDiag-MIPS-Gianma.csv"
pad = 3
param_label = "Pe"

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel(param_label)
ax.set_xlabel('Density')

# test_name = sys.argv[2]

data = np.loadtxt(file, delimiter=',', skiprows=1)
rho_gases_prev = data[:, -2]
rho_liquids_prev = data[:,-1]
vars_prev = data[:,-3]
label=r"$L_x=200, L_y=200, \rho_0=0.9$"
if mode == 'pfap':
    vars_prev = data[:,-4] / vars_prev * np.float_power(2, 1/6)
# Plot binodal from previous
ax.errorbar(rho_gases_prev, vars_prev, color='grey', ms=marker_size, ls='-', marker='.', label=label)
ax.errorbar(rho_liquids_prev, vars_prev, color='grey', ms=marker_size, ls='-', marker='.')


# medium homo, dilute, aspect ratio 1:1
file = "pfaps_test26_aspect1_dilute_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=200, L_y=200, \rho_0=0.65$"
    color='C4'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='s', fillstyle='none', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='s', fillstyle='none')

# medium homo, dilute, aspect ratio 1:1
file = "pfaps_test34_scalefast_binodal_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=200, L_y=200, \rho_0=0.75, v=5, \epsilon=1.25$"
    color='C5'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='s', fillstyle='none', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='s', fillstyle='none')

file = "pfaps_test27_redo_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=200, L_y=200, \rho_0=0.9$"
    color='C6'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='d', fillstyle='none', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='d', fillstyle='none')

# file = test_name + '_histo_phase_diagram'
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_hist = data[:, -2]
#     rho_liquids_hist = data[:,-1]
#     vars_hist = data[:,-3]
#     # Plot binodal from histogram
#     ax.errorbar(rho_gases_hist, vars_hist, color='blue', ls='', marker='.', label="Gas (histogram)")
#     ax.errorbar(rho_liquids_hist, vars_hist, color='orange', ls='', marker='.', label = "Liquid (histogram)")

# file = 'pfaps_slab/pfaps_test17_diagram_slab_phase_diagram'
# # this one takes density too high, the gas phase is unstable.
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_slab = data[:, -2]
#     rho_liquids_slab = data[:,-1]
#     vars_slab = data[:,-3]
#     # Plot binodal from slab
#     ax.errorbar(rho_gases_slab, vars_slab, color='purple', ls='', marker='.', label="Gas (slab)")
#     ax.errorbar(rho_liquids_slab, vars_slab, color='red', ls='', marker='.', label="Liquid (slab)")

# small homo
file = "pfaps_test21_diagram_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=248, L_y=100, \rho_0=0.9$"
    color='C1'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='v', fillstyle='none', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='v', fillstyle='none')


# medium homo
file = "pfaps_test23_diagram_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label=r"$L_x=500, L_y=200, \rho_0=0.9$"
    color='C0'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='+', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='+')

# medium homo, dilute
file = "pfaps_test25_dilute_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=500, L_y=200, \rho_0=0.65$"
    color='C3'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='x', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='x')

# small slab
# file = "pfaps_test19_diagram_histo_phase_diagram"
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_hist = data[:, -2]
#     rho_liquids_hist = data[:,-1]
#     vars_hist = data[:,-3]
#     label = "Lx = 248, slab"
#     # Plot binodal from histogram
#     ax.errorbar(rho_gases_hist, vars_hist, ls='', marker='.', label=label)
#     ax.errorbar(rho_liquids_hist, vars_hist, ls='', marker='.', label = label)

# big slab
# file = "pfaps_test20_diagram_histo_phase_diagram"
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_slab = data[:, -2]
#     rho_liquids_slab = data[:,-1]
#     vars_slab = data[:,-3]
#     label = "Lx = 1000, slab"
#     # Plot binodal from slab
#     ax.errorbar(rho_gases_slab, vars_slab, ls='', marker='.', label=label)
#     ax.errorbar(rho_liquids_slab, vars_slab, ls='', marker='.', label=label)

# big homo
file = "pfaps_test22_diagram_histo_phase_diagram"
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_slab = data[:, -2]
    rho_liquids_slab = data[:,-1]
    vars_slab = data[:,-3]
    label = r"$L_x=1000, L_y=400, \rho_0=0.9$"
    color='C2'
    ax.errorbar(rho_gases_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='o', fillstyle='none', label=label)
    ax.errorbar(rho_liquids_slab, vars_slab, color=color, ms=marker_size, ls='-', marker='o', fillstyle='none')

ax.set_title("Homogeneous inital conditions")
ax.legend()
# ax.set_xlim(0.0, 1.2)
# ax.set_ylim(14, 33)
plt.savefig('pfap_test_finite_size_compare_phase_diagram_homo.png', dpi=300, bbox_inches='tight')
plt.close()