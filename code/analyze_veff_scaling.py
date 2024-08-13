import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

def veff_fit(rho, coef, s1, s2, s3, s4):
    return coef/2 * rho * (1-s1*rho+s2*rho**2) * (1-np.tanh(s3*(rho-s4)))

fig, ax = plt.subplots()
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$p^{IK}(\rho)/(v*\delta)$')

compare = "pressure_pfaps_2018_Pa_data"
data = np.loadtxt(compare, skiprows=1)
rhos = data[:,0]
sigmas = data[:,1]
ax.plot(rhos, sigmas, label=f"Comparison: Pe=13.3", marker='.', ms=10)

i = 0
v = 24

test_name = "pfaps_test40_homo_compare_Pe13.3"
test_Pe = 40/3
# delta = penetration(v, 0.125)
data = np.loadtxt(test_name+'_veff_data', skiprows=1)
rhos = data[:,0]
vs = data[:,1]
# coef = v0*lp/2 (lp=Pe)
# popt, pcov = optimize.curve_fit(veff_fit, rhos, rhos*vs)
ax.plot(rhos, rhos*vs*test_Pe/2, label=f"Pe = {test_Pe:.1f}", marker='.', ms=10)

test_name = test_name = "pfaps_test40_homo_compare_Pe6.7"
test_Pe = 20/3
# delta = penetration(v, 0.125)
data = np.loadtxt(test_name+'_veff_data', skiprows=1)
rhos = data[:,0]
vs = data[:,1]
# coef = v0*lp/2 (lp=Pe)
# popt, pcov = optimize.curve_fit(veff_fit, rhos, rhos*vs)
ax.plot(rhos, 2*rhos*vs*test_Pe/2, label=f"Pe = {test_Pe:.1f}", marker='.', ms=10)

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.set_title("Active pressure (collapsed) as a function of density")
ax.legend()
plt.savefig('pfap_v24_eps1_v.png',  dpi=300, bbox_inches='tight')
