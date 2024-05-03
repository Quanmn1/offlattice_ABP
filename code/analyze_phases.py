import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import os
import sys
from helper import *
from analyze_density import *
from analyze_histogram import *

start = float(sys.argv[3])
end = float(sys.argv[4])
space = float(sys.argv[5])
test_name = sys.argv[1]
mode = sys.argv[2]
num_segments = int(sys.argv[6])
init = sys.argv[7]
if len(sys.argv) > 8:
    fit = sys.argv[8]
else:
    fit = 'gauss'


vars_hist, rho_gases_hist, rho_liquids_hist, rho_gases_std_hist, rho_liquids_std_hist, param_label = analyze_histogram(test_name, mode, start, end, space, num_segments, fit=fit)
if init == "slab":
    vars_slab, rho_gases_slab, rho_liquids_slab, rho_gases_std_slab, rho_liquids_std_slab, param_label = analyze_slab(test_name, mode, start, end, space, num_segments)

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel(param_label)
ax.set_xlabel('Density')
# Plot binodal from histogram
ax.errorbar(rho_gases_hist, vars_hist, xerr=rho_gases_std_hist, color='blue', ls='', marker='.', label="Gas (histogram)")
ax.errorbar(rho_liquids_hist, vars_hist, xerr=rho_liquids_std_hist, color='orange', ls='', marker='.', label = "Liquid (histogram)")
if init == "slab":
    # Plot binodal from slab
    ax.errorbar(rho_gases_slab, vars_slab, xerr=rho_gases_std_slab, color='purple', ls='', marker='.', label="Gas (slab)")
    ax.errorbar(rho_liquids_slab, vars_slab, xerr=rho_liquids_std_slab, color='red', ls='', marker='.', label="Liquid (slab)")
ax.set_title("Phase diagram")
ax.legend()
# ax.set_xlim(0.0, 1.2)
# ax.set_ylim(14, 33)
plt.savefig(test_name + '_phase_diagram.png', dpi=300, bbox_inches='tight')
plt.close()