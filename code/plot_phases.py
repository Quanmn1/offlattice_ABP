import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

mode = sys.argv[1]

if mode == "qsap":
    pad = 1
    param_label = r"$\lambda$"
    file = "PhaseDiag-MIPS-QSAP.csv"
elif mode == "pfap":
    pad = 3
    param_label = r"$l_p/r_f$"
    # file = "PhaseDiag-MIPS-Gianma.csv"
    file = "PhaseDiag-MIPS-Yongfeng"
# file_anal = "phase_diag_qsap_rmap"

marker_size = 10

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel(param_label)
ax.set_xlabel('Density')

test_name_histo = sys.argv[2]
test_name_slab = sys.argv[3]

# data = np.loadtxt(file, delimiter=' ', skiprows=1)
# skip=0
# rho_gases_prev = data[:6-skip, 3]
# rho_liquids_prev = data[:6-skip,1]
# rho_gases_prev_std = data[:6-skip, 4]
# rho_liquid_prev_std = data[:6-skip, 2]
# Drs = data[:6-skip,0]
# vars_prev = 5/Drs
# # Plot binodal from previous
# ax.errorbar(rho_gases_prev, vars_prev, xerr=rho_gases_prev_std, color='green', ms=marker_size, marker='.', label="Yongfeng's")
# ax.errorbar(rho_liquids_prev, vars_prev, xerr=rho_liquid_prev_std, color='green', ms=marker_size, marker='.')

# data = np.loadtxt(file_anal, skiprows=1)
# rho_gases_prev = data[:, -2]*50
# rho_liquids_prev = data[:,-1]*50
# vars_prev = data[:,-3]
# if mode == 'pfap':
#     vars_prev = data[:,-4] / vars_prev * np.float_power(2, 1/6)
# # Plot binodal from r mapping
# ax.errorbar(rho_gases_prev, vars_prev, color='grey', ls='-', ms=marker_size, marker='', label="R mapping")
# ax.errorbar(rho_liquids_prev, vars_prev, color='grey', ls='-', ms=marker_size, marker='')


file = test_name_histo + '_histo_phase_diagram'
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases_hist = data[:, 1]
    rho_liquids_hist = data[:,2]
    vars_hist = data[:,0]
    # Plot binodal from histogram
    ax.errorbar(rho_gases_hist, vars_hist, color='brown', ls='', ms=marker_size, marker='+', label="Homo, "+r"$L_x=L_y=200$")
    ax.errorbar(rho_liquids_hist, vars_hist, color='brown', ls='', ms=marker_size, marker='+')

# file = test_name_histo + '_histo_phase_diagram_full_fit'
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_hist = data[:, -2]
#     rho_liquids_hist = data[:,-1]
#     vars_hist = data[:,-3]
#     # Plot binodal from histogram
#     ax.errorbar(rho_gases_hist, vars_hist, color='red', ls='', ms=marker_size, marker='.', label="Full fit, Homo, "+r"$L_x=L_y=200$")
#     ax.errorbar(rho_liquids_hist, vars_hist, color='red', ls='', ms=marker_size, marker='.')

file = test_name_slab + '_density_pressure'
if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    skip=5
    rho_gases_slab = data[:, 1]
    rho_liquids_slab = data[:,2]
    vars_slab = data[:,0]
    # Plot binodal from slab
    ax.errorbar(rho_gases_slab, vars_slab, color='purple', ls='', ms=marker_size, marker='x', label="Slab, "+r"$L_x=400,L_y=200$")
    ax.errorbar(rho_liquids_slab, vars_slab, color='purple', ls='', ms=marker_size, marker='x')

ax.set_title("Phase diagram for " + r"$\epsilon=100$")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend()
# ax.set_xlim(0.0, 1.2)
# ax.set_ylim(14, 33)
plt.savefig('pfap_harmonic_compare_phase_diagram_largeeps.png', dpi=300, bbox_inches='tight')
plt.close()