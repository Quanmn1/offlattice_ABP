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

# def sigma_fit(rho, d1, d2, d3, d4):
#     return (d2 * rho**2 + d3 * rho**3 + d4 * rho**4)/np.sqrt(1-rho/d1)

def Pa_fit(rho, s1, s2, s3, s4):
    return lp * rho * v/2 * (1-s1*rho+s2*rho**2) * (1 - np.tanh(s3*(rho-s4)))/2

# def smooth(x, zero, slope):
#     return 1 - np.exp(-1/slope/(x-zero)**2) * np.heaviside(x-zero, 0)
    
# def Pa_fit(rho, s1, s2, s3, s4):
#     lp = 10
#     v = 5
#     return lp * rho * v/2 * (1-s1*rho+s2*rho**2) * smooth(rho, s4, s3)

def veff_fit(rho, s1, s2, s3, s4):
    v = 5
    return v * (1-s1*rho+s2*rho**2) * (1-np.tanh(s3*(rho-s4)))/2

rho_dense = np.linspace(0, 1.32, 1000)
fig, ax = plt.subplots()
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$p(\rho)$')

# compare = "pressure_pfaps_2018_sigma_data"
# data = np.loadtxt(compare, skiprows=1)
# rhos = data[:,0]
# sigmas = data[:,1]
# color = 'C0'
# ax.plot(rhos, sigmas, label=f"2018 data", marker='.', ls='', c=color, ms=10)
# popt, pcov = optimize.curve_fit(sigma_fit, rhos, sigmas, p0=(0.5,4, 2e-9, 10))
# print("SigmaIK 2018")
# print(popt)
# ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

# compare = "pressure_pfaps_2018_Pa_data"
# lp = 13.33
# data = np.loadtxt(compare, skiprows=1)
# rhos = data[:,0]
# Pa = data[:,1]
# ax.plot(rhos, Pa, marker='.', ls='', ms=10, c=color)
# popt, pcov = optimize.curve_fit(Pa_fit, rhos, Pa, p0=(1,-0.2,0.4,0.3))
# print("veff 2018")
# print(popt)
# ax.plot(rho_dense, Pa_fit(rho_dense, *popt), c=color)

v=1.8393972058572117
Dr = 0.2
lp = v/Dr
tests = sys.argv[1:]
i = 0
num_params = 1
while i < len(tests):
    test_name = tests[i]
    # epsilon = float(tests[i+1])
    # dt = float(tests[i+1])
    # rf = float(tests[i+2])
    # delta = penetration(v, 0.125)
    test_lp = 10
    color = f'C{i//num_params}'
    data = np.loadtxt(test_name+'_sigmaIK_data', skiprows=1)
    rhos = data[:,0]
    sigmas = data[:,2]
    sigmas_std = data[:,3]
    marker = '.'
    indices = rhos<1.5
    ax.errorbar(rhos[indices], sigmas[indices], yerr=sigmas_std, label=r"$p_{IK}$", marker=marker, ls='', ms=10, c=color)
    try:
        popt, pcov, infodict, mesg, ier = optimize.curve_fit(sigma_fit, rhos, sigmas, sigma=sigmas_std, absolute_sigma=True, p0=(0.1, 4, 0.1, 2), full_output=True)
    except RuntimeError:
        print(f"Could not fit to SigmaIK of {test_name}")
    else:
        chisq = np.sum(np.square(infodict['fvec']))
        dof = len(rhos) - 4
        print("SigmaIK " + test_name)
        print(popt)
        print(pcov)
        print(f"chi^2={chisq}/dof={dof}\n")
        ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

    color = f'C{i//num_params+1}'
    data = np.loadtxt(test_name+'_sigmaA_data', skiprows=1)
    rhos = data[:,0]
    Pas = data[:,2]
    Pas_std = data[:,3]
    indices = rhos<1.5
    ax.errorbar(rhos[indices], Pas[indices], yerr=Pas_std, label=r"$p_A$", marker=marker, ls='', ms=10, c=color)
    try:
        popt, pcov, infodict, mesg, ier = optimize.curve_fit(Pa_fit, rhos, Pas, sigma=Pas_std, absolute_sigma=True, p0=(1.1, 0.21, 3, 0.8), full_output=True)
    except RuntimeError:
        print(f"Could not fit to SigmaA of {test_name}")
    else:
        chisq = np.sum(np.square(infodict['fvec']))
        dof = len(rhos) - 4
        print("SigmaA " + test_name)
        print(popt)
        print(pcov)
        print(f"chi^2={chisq}/dof={dof}\n")
        ax.plot(rho_dense, Pa_fit(rho_dense, *popt), c=color)

    # color = f'C{i//num_params+2}'
    # data = np.loadtxt(test_name+'_veff_data', skiprows=1)
    # rhos = data[:,0]
    # veffs = data[:,1]
    # veffs_std = data[:,2]
    # ax.errorbar(rhos[indices], veffs[indices], yerr=veffs_std, label=r"$v_{eff}$", marker=marker, ls='', ms=10, c=color)
    # try:
    #     popt, pcov, infodict, mesg, ier = optimize.curve_fit(veff_fit, rhos, veffs, sigma=veffs_std, absolute_sigma=True, p0=(1,0.4,5,1.3), full_output=True)
    # except RuntimeError:
    #     print(f"Could not fit to SigmaA of {test_name}")
    # else:
    #     chisq = np.sum(np.square(infodict['fvec']))
    #     dof = len(rhos) - 4
    #     print("veff " + test_name)
    #     print(popt)
    #     print(pcov)
    #     print(f"chi^2={chisq}/dof={dof}\n")
    #     ax.plot(rho_dense, veff_fit(rho_dense, *popt), c=color)

    i += num_params

ax.set_xlim(left=0)
# ax.set_ylim(bottom=0)
ax.set_title(f"Direct and active pressure for " + fr"$l_p = {test_lp}$")
ax.legend()
plt.savefig('pfap_harmonic_vlowqsap_pressures_fit.png',  dpi=300, bbox_inches='tight')
