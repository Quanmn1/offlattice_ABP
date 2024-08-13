import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import os
import sys
from helper import *
from analyze_slab import *
from analyze_histogram import *

test_name = sys.argv[1]
mode = sys.argv[2]
vars = np.array(sys.argv[3].split(),dtype=float)
num_segments = int(sys.argv[4])
init = sys.argv[5]
homo_fit = sys.argv[6]
if init == "slab":
    slab_fit = sys.argv[7]

vars_hist, rho_gases_hist, rho_liquids_hist, rho_gases_std_hist, rho_liquids_std_hist, param_label = analyze_histogram(test_name, mode, vars, num_segments, fit=homo_fit)
if init == "slab":
    vars_slab, rho_gases_slab, rho_liquids_slab, rho_gases_std_slab, rho_liquids_std_slab, param_label = analyze_slab(test_name, mode, vars, num_segments)

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel(param_label)
ax.set_xlabel('Density')
# Plot binodal from histogram
ax.errorbar(rho_gases_hist, vars_hist, xerr=3*rho_gases_std_hist, color='blue', ls='', marker='.', label="Gas (histogram)")
ax.errorbar(rho_liquids_hist, vars_hist, xerr=3*rho_liquids_std_hist, color='orange', ls='', marker='.', label = "Liquid (histogram)")
if init == "slab":
    # Plot binodal from slab
    ax.errorbar(rho_gases_slab, vars_slab, xerr=3*rho_gases_std_slab, color='purple', ls='', marker='.', label="Gas (slab)")
    ax.errorbar(rho_liquids_slab, vars_slab, xerr=3*rho_liquids_std_slab, color='red', ls='', marker='.', label="Liquid (slab)")
ax.set_title("Phase diagram")
ax.legend()
ax.set_xlim(left=0)
# ax.set_ylim(14, 33)
plt.savefig(test_name + '_phase_diagram.png', dpi=300, bbox_inches='tight')
plt.close()

if init == "slab":
    # in pfaps_test44: pressures_gas_tot (read from _sigma files) is sum of actual pressures and Qxx. 
    # Qxx in test44 is NOT nematic tensor, it's the traceful nematic tensor.
    vars_data, pressures_gas, pressures_liquid, gas_std, liquid_std, param_label = analyze_slab(test_name, mode, vars, num_segments, data="sigma", fit=slab_fit)
    # vars_data, pressures_gas_Q, pressures_liquid_Q, std, std, param_label = analyze_slab(test_name, mode, vars, num_segments, data="Qxx")

    # pressures_gas = pressures_gas_tot - pressures_gas_Q
    # pressures_liquid = pressures_liquid_tot - pressures_liquid_Q

    with open(test_name + f'_density_pressure', 'w') as f:
        f.write(f"{param_label} rho_gas rho_liquid pressure_gas pressure_liquid pressure_gas_std pressure_liquid_std\n")
        for i in range(len(vars)):
            f.write(f"{vars_data[i]:.2f} {rho_gases_slab[i]:.4f} {rho_liquids_slab[i]:.4f} {pressures_gas[i]:.4f} {pressures_liquid[i]:.4f} {gas_std[i]:.4f} {liquid_std[i]:.4f}\n")
