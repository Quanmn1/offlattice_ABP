import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

"""
Plot scaled pressure vs measured pressure
"""
v0 = 5.0

def vstar(rho):
    return v0 * (1 - 0.96436447 * rho + 0.20893963 * rho**2) * (1 - np.tanh(6.5619715 * (rho - 1.06965412)))/2

def pA(rho, Dr, v0=v0):
    return v0 * vstar(rho) * rho / (2 * Dr)

def pIK(rho):
    return 6.71615706e-2 * (np.exp(4.32838715 * rho) - 1) + 5.91672633e-9 * (np.exp(14.8247273 * rho) - 1)

def p(rho, Dr):
    return pA(rho, Dr) + pIK(rho)

data = np.loadtxt("pfaps_test47_harmonic_slab_largeeps_density_pressure", skiprows=1)
Drs_measured = np.array([5])/(data[:, 0])
Drs_measured_targetted = [0.02, 0.05, 0.08, 0.10, 0.15]
rho_g = data[:, 1]
rho_l = data[:, 2]
rho_g_std = data[:, 3]
rho_l_std = data[:, 4]
p_g = data[:, 5]
p_l = data[:, 6]
p_g_std = data[:, 7]
p_l_std = data[:, 8]
p_scaled_g = p(rho_g, Drs_measured)
p_scaled_l = p(rho_l, Drs_measured)
p_scaled_g_std = (p(rho_g+rho_g_std, Drs_measured) - p(rho_g-rho_g_std, Drs_measured))/2
p_scaled_l_std = (p(rho_l+rho_l_std, Drs_measured) - p(rho_l-rho_l_std, Drs_measured))/2

p_avg = (p_g + p_l)/2

# plot the relative error of the scaled pressure relative to the measured pressure
separation = 0.001
plt.errorbar(Drs_measured-separation, p_scaled_g/p_avg, yerr=p_scaled_g_std/p_avg, ls='', label="Scaled gas pressure")
plt.errorbar(Drs_measured+separation, p_scaled_l/p_avg, yerr=p_scaled_l_std/p_avg, ls='', label="Scaled liquid pressure")

# plot the relative error of the binodals obtained by R-mapping
data_theory = np.loadtxt("pfap_harmonic_rmap", skiprows=1)
Drs = data_theory[:, 0]
indices = np.searchsorted(Drs, Drs_measured_targetted)
rho_theory_g = data_theory[indices, 1]
rho_theory_l = data_theory[indices, 2]
error_g = (rho_theory_g - rho_g)/rho_g
error_l = (rho_theory_l - rho_l)/rho_l
# plt.plot(Drs_measured-separation, 1+error_g, ls='', ms=5, marker='o', label="Gas density")
plt.plot(Drs_measured+separation, 1+error_l, ls='', ms=5, marker='o', label="Liquid density")

# density uncertainty
plt.errorbar(Drs_measured-separation/2, rho_g/rho_g, rho_g_std/rho_g, ls='', label="Gas density uncert")
plt.errorbar(Drs_measured+separation/2, rho_g/rho_g, rho_l_std/rho_l, ls='', label="Liquid density uncert")

# for pressures in (p_scaled_g, p_scaled_l):
#     plt.scatter(Drs_measured, pressures, s=size)
# plt.errorbar(Drs_measured, p_g, yerr=p_g_std, ls='', marker='o', ms=size)
# plt.errorbar(Drs_measured, p_l, yerr=p_l_std, ls='', marker='o', ms=size)

# plt.legend(["Scaled gas", "Scaled liquid", "Measured gas", "Measured liquid"])
plt.legend(loc=(1.05,0.2))
plt.savefig("pfap_harmonic_pressure_vs_binodal.png", dpi=300, bbox_inches="tight")
