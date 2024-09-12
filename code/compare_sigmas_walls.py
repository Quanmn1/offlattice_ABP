import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

def sigma_fit(rho, d1, d2, d3, d4):
    return d1 * (np.exp(d2*rho)-1) + d3 * (np.exp(d4*rho)-1)

# def sigma_fit(rho, d1, d2):
#     return d1 * (np.exp(d2*rho)-1)

def Pa_fit(rho, s1, s2, s3, s4):
    return 13.3/2 * rho * 24/2 * (1-s1*rho+s2*rho**2) * (1-np.tanh(s3*(rho-s4)))

rho_dense = np.linspace(0, 1.3, 1000)
fig, ax = plt.subplots()
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$p(\rho)$')

# compare = "pressure_pfaps_2018_sigma_data"
# data = np.loadtxt(compare, skiprows=1)
# rhos = data[:,0]
# sigmas = data[:,1]
# color = 'C0'
# ax.plot(rhos, sigmas, label=f"Comparison: lp=13.3", marker='.', ls='', c=color, ms=10)
# popt, pcov = optimize.curve_fit(sigma_fit, rhos, sigmas, p0=(0.5,4, 2e-9, 10))
# print(popt)
# ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

# compare = "pressure_pfaps_2018_Pa_data"
# lp = 13.33
# data = np.loadtxt(compare, skiprows=1)
# rhos = data[:,0]
# Pa = data[:,1]
# ax.plot(rhos, Pa, marker='.', ls='', ms=10, c=color)
# popt, pcov = optimize.curve_fit(Pa_fit, rhos, Pa, p0=(1,-0.2,0.4,0.3))
# print(popt)
# ax.plot(rho_dense, Pa_fit(rho_dense, *popt), c=color)

compare = "PhaseDiag-MIPS-Yongfeng"
data = np.loadtxt(compare, skiprows=1)
rhos = data[:,0]
sigmas = data[:,1]
color = 'C0'
ax.plot(rhos, sigmas, label=f"Comparison: lp=13.3", marker='.', ls='', c=color, ms=10)
popt, pcov = optimize.curve_fit(sigma_fit, rhos, sigmas, p0=(0.5,4, 2e-9, 10))
print(popt)
ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

lp=10
tests = sys.argv[1:]
i = 0
num_params = 2
while i < len(tests):
    test_name = tests[i]
    # epsilon = float(tests[i+1])
    omega = float(tests[i+1])
    # rf = float(tests[i+2])
    # delta = penetration(v, 0.125)
    test_lp = 10
    color = f'C{i//num_params+1}'
    i += num_params
    data = np.loadtxt(test_name+'_total_pressure', skiprows=1)
    rhos = data[:,0]
    sigmas = data[:,2]
    lefts = data[:,4]
    rights = data[:,6]
    marker = '.'
    indices = rhos<1.5
    ax.errorbar(rhos[indices], sigmas[indices], xerr=3*data[indices,1], yerr=3*data[indices,3], label=f"Bulk, stiffness = {omega}", marker=marker, ls='', ms=2, c=color)
    ax.errorbar(rhos[indices], lefts[indices], xerr=3*data[indices,1], yerr=3*data[indices,5], label=f"Left wall, stiffness = {omega}", marker=marker, ls='', ms=2, c=color)
    ax.errorbar(rhos[indices], rights[indices], xerr=3*data[indices,1], yerr=3*data[indices,7], label=f"Right wall, stiffness = {omega}", marker=marker, ls='', ms=2, c=color)

    popt, pcov = optimize.curve_fit(sigma_fit, rhos, sigmas, p0=(0.5,4, 2e-9, 10))
    # print(popt)
    ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

    data = np.loadtxt(test_name+'_veff_data', skiprows=1)
    rhos = data[:,0]
    veffs = data[:,1]
    Pas = lp*veffs*rhos/2
    ax.plot(rhos[indices], Pas[indices], marker=marker, ls='', ms=10, c=color)
    popt, pcov = optimize.curve_fit(Pa_fit, rhos, Pas, p0=(1,-0.2,0.4,0.3))
    print(popt)
    ax.plot(rho_dense, Pa_fit(rho_dense, *popt), c=color)


ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.set_title(f"Direct and active pressure for l_p = {test_lp}, epsilon = 100")
ax.legend()
plt.savefig('pfap_harmonic_wall_pressures.png',  dpi=300, bbox_inches='tight')
