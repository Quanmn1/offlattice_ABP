import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

def area_fraction(density, rf):
    return density * np.pi/4 * rf**2

def ycoord(rf):
    if config == 0:
        lp = 1
        return lp/rf
    elif config == 1:
        return rf
    
def xcoord(rho, rf):
    if config == 0:
        return area_fraction(rho, rf)
    elif config == 1:
        return rho

config = 1

matplotlib.rcParams.update({'font.size': 16})
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
# mode = sys.argv[1]

# if mode == "qsap":
#     pad = 1
#     param_label = r"$\lambda$"
#     file = "PhaseDiag-MIPS-QSAP.csv"
# elif mode == "pfap":
#     pad = 3
#     param_label = r"$\ell_p/r_f$"
#     # param_label = r"$r_f$"
#     # file = "PhaseDiag-MIPS-Gianma.csv"
#     file = "PhaseDiag-MIPS-Yongfeng"

file_anal = "pfap_qsap_gradient_liquid"
file_anal2 = "pfap_new_theory"
# file_anal3 = "pfap_harmonic_xi1"
# file_anal4 = "pfap_harmonic_xi100"
# file_anal5 = "pfap_harmonic_rmap_symbolic"
# file_anal = "pfap_harmonic_stiff200_rmap"

marker_size = 5

fig, ax = plt.subplots(figsize = (8, 10))

# test_name_histo = "pfaps_qsaps_denser_stable"
# test_name_slab = "pfaps_qsaps_denser_stable"

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

data = np.loadtxt(file_anal, skiprows=1)
skip_high = 0
skip_low = 0
rho_gases_prev = data[skip_high:, 1]
rho_liquids_prev = data[skip_high:,2]
vars_prev = data[skip_high:,0]
# Plot binodal from r mapping
ax.errorbar(rho_gases_prev, vars_prev, color='grey', ls='-', ms=marker_size, marker='', label="old v(rho)")
ax.errorbar(rho_liquids_prev, vars_prev, color='grey', ls='-', ms=marker_size, marker='')

data = np.loadtxt(file_anal2, skiprows=1)
skip_high = 0
skip_low = 0
rho_gases_prev = data[skip_high:, 1]
rho_liquids_prev = data[skip_high:,2]
vars_prev = data[skip_high:,0]
# Plot binodal from r mapping
ax.errorbar(rho_gases_prev, vars_prev, color='orange', ls='-', ms=marker_size, marker='', label="new v(rho)")
ax.errorbar(rho_liquids_prev, vars_prev, color='orange', ls='-', ms=marker_size, marker='')

# data = np.loadtxt(file_anal3, skiprows=0)
# skip_high = 0
# skip_low = 0
# rho_gases_prev = data[skip_high:, 1]
# rho_liquids_prev = data[skip_high:,2]
# vars_prev = data[skip_high:,0]
# if mode == 'pfap':
#     vars_prev = 5 / vars_prev
# # Plot binodal from modification
# ax.errorbar(rho_gases_prev, vars_prev, color='blue', ls='-', ms=marker_size, marker='', label=r"$\xi=1$")
# ax.errorbar(rho_liquids_prev, vars_prev, color='blue', ls='-', ms=marker_size, marker='')


# data = np.loadtxt(file_anal4, skiprows=0)
# skip_high = 0
# skip_low = 0
# rho_gases_prev = data[skip_high:, 1]
# rho_liquids_prev = data[skip_high:,2]
# vars_prev = data[skip_high:,0]
# if mode == 'pfap':
#     vars_prev = 5 / vars_prev
# # Plot binodal from modification
# ax.errorbar(rho_gases_prev, vars_prev, color='pink', ls='-', ms=marker_size, marker='', label=r"$\xi=100$")
# ax.errorbar(rho_liquids_prev, vars_prev, color='pink', ls='-', ms=marker_size, marker='')

# data = np.loadtxt(file_anal5, skiprows=0)
# skip_high = 0
# skip_low = 0
# rho_gases_prev = data[skip_high:, 1]
# rho_liquids_prev = data[skip_high:,2]
# vars_prev = data[skip_high:,0]
# if mode == 'pfap':
#     vars_prev = 5 / vars_prev
# # Plot binodal from r mapping
# ax.errorbar(rho_gases_prev, vars_prev, color='black', ls='-', ms=marker_size, marker='', label="R mapping with EoS from symbolic regression")
# ax.errorbar(rho_liquids_prev, vars_prev, color='black', ls='-', ms=marker_size, marker='')


# file = test_name_histo
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     skip = 0
#     indices = (data[:,0] < 0.21) * (data[:,0] > 0.15)
#     rho_gases_hist = data[indices, 1]
#     rho_liquids_hist = data[indices,2]
#     rfs = data[indices,0]
    
#     # Plot binodal from histogram
#     ax.plot(xcoord(rho_gases_hist, rfs), ycoord(rfs), ls='--', color='C0', ms=marker_size, marker='o', zorder=10, clip_on=False)
#     ax.plot(xcoord(rho_liquids_hist,rfs), ycoord(rfs), ls='--', color='C0', ms=marker_size, marker='o', label="Histogram, Lx=Ly=20")
#     # ax.fill_betweenx(vars_hist, rho_gases_hist, rho_liquids_hist, alpha=0.3, color="purple")

# file = test_name_histo + '_histo_phase_diagram_full_fit'
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     rho_gases_hist = data[:, -2]
#     rho_liquids_hist = data[:,-1]
#     vars_hist = data[:,-3]
#     # Plot binodal from histogram
#     ax.errorbar(rho_gases_hist, vars_hist, color='red', ls='', ms=marker_size, marker='.', label="Full fit, Homo, "+r"$L_x=L_y=200$")
#     ax.errorbar(rho_liquids_hist, vars_hist, color='red', ls='', ms=marker_size, marker='.')

# file = test_name_slab + '_slab_density'
# if os.path.exists(file):
#     data = np.loadtxt(file, skiprows=1)
#     indices = (data[:,0] < 0.21) * (data[:,0] > 0.15)
#     rho_gases_slab = data[indices, 1]
#     rho_liquids_slab = data[indices,2]
#     rfs = data[indices,0]
#     # Plot binodal from slab
#     ax.plot(xcoord(rho_gases_slab, rfs), ycoord(rfs), ls='--', color='C1', ms=marker_size, marker='o', zorder=10, clip_on=False)
#     ax.plot(xcoord(rho_liquids_slab,rfs), ycoord(rfs), ls='--', color='C1', ms=marker_size, marker='o', label="Slab, Lx=70, Ly=40")

    # ax.errorbar(, vars_slab, color='purple', ls='', ms=marker_size, marker='x', label="Numerics")
    # ax.errorbar(rho_liquids_slab, vars_slab, color='purple', ls='', ms=marker_size, marker='x')
lp = 5
# ax.axhline(lp/0.154, color="red", label="My estimate of " + r"$r_{12}$")
# ax.axhline(0.152, color="C1", label="Prediction of phase boundary")
# ax.axhline(0.215, color="C1")
# ax.axvline(25, ls="--", color="C3")
plt.legend()
# ax.set_title("PFAPs: theory and numerics")    
if config == 0:
    ax.set_xlim(left=0, right=1.2)
    ax.set_ylim(bottom=0, top=150)
    plotname = 'pfqs_phase_diagram_absorbing.png'
    ax.set_ylabel(r"$\ell_p/r_f$")
    ax.set_xlabel(r'$\phi=\pi\rho r_f^2/4$')
elif config == 1:
    # ax.axvline(x=25,ls='--',color="C3")
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=0)
    plotname = 'pfaps_qsaps_denser_compare_big_theory.png'
    ax.set_ylabel(r"$r_f$")
    ax.set_xlabel(r'$\rho$')

plt.savefig(plotname, dpi=300, bbox_inches='tight')
plt.close()