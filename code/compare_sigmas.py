import numpy as np
from scipy import optimize, stats, signal
import math
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
    lp = 12.5
    v0 = 2
    return rho * v0/2 * (1-s1*rho+s2*rho**2) * (1 - np.tanh(s3*(rho-s4)))/2 * lp

# def smooth(x, zero, slope):
#     return 1 - np.exp(-1/slope/(x-zero)**2) * np.heaviside(x-zero, 0)
    
# def Pa_fit(rho, s1, s2, s3, s4):
#     lp = 10
#     v = 5
#     return lp * rho * v/2 * (1-s1*rho+s2*rho**2) * smooth(rho, s4, s3)

# def veff_fit(rho, s1, s2, s3, s4):
#     v = 2
#     return v * (1-s1*rho+s2*rho**2) * (1-np.tanh(s3*(rho-s4)))/2


"""
Scaled pIK and pA from PFAPs to PFQS
"""
def v(rho, v0, lamb):
    return v0 * np.exp(-lamb * np.tanh((rho - 25) / 10))

# def pIK1(rho, rf, v0, lamb): # measured with v0
#     return (6.71615706e-2 * (np.exp(4.32838715 * rho*rf**2) - 1) + 5.91672633e-9 * (np.exp(14.8247273 * rho*rf**2) - 1)) / rf

# def pIK2(rho, rf, v0, lamb): # scaled from p measured with v0
#     return v(rho, v0, lamb)/v0 * (6.71615706e-2 * (np.exp(4.32838715 * rho*rf**2) - 1) + 5.91672633e-9 * (np.exp(14.8247273 * rho*rf**2) - 1)) / rf

# def pIK3(rho, rf, v0, lamb): # measured with v(\infty)
#     return (v0/np.exp(lamb))/(5/np.exp(1)) * (rho*rf**2)**2 * (rho*rf**2/0.61 + np.exp((rho*rf**2)**6-1.81) + 0.0604) / rf

# def pIK4(rho, rf, v0, lamb): # scaled from p measured with v(\infty)
#     return v(rho, v0, lamb)/(5/np.exp(1)) * (rho*rf**2)**2 * (rho*rf**2/0.61 + np.exp((rho*rf**2)**6-1.81) + 0.0604) / rf

# def pA1(rho, rf, v0, lamb, Dr):
#     return v(rho, v0, lamb) * v(rho, v0, lamb) * (1-np.tanh((rho*rf**2)**6-rho*rf**2)) * (1-np.tanh(0.783*rho*rf**2))**2 * rho /(2*Dr)

def pIK1(rho, rf, v0, lamb):
    return v(rho, v0, lamb) / (5/math.e) * 0.07879 * (np.exp(6.153 * rho * rf**2) - 1)

def pIK2(rho, rf, v0, lamb):
    return 0.0680 * (np.exp(6.22 * rho*rf**2)-1)

def pIK3(rho, rf, v0, lamb):
    return v(rho, v0, lamb) / (5/math.e) / rf * 0.576* (rho*rf**2 * (np.exp(rho*rf**2)-1) + (np.exp((rho*rf**2)**5)-1))

# def veff1(rho, rf, v0, lamb):


rho_dense = np.linspace(0, 60, 1000)
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


tests = sys.argv[1:]
i = 0
num_params = 2
while i < len(tests):
    test_name = tests[i]
    rf = float(tests[i+1])
    lamb = 0.4
    v0 = 5*math.exp(lamb-1)
    Dr = v0/8.379/rf
    color = f'C{i//num_params}'
    data = np.loadtxt(test_name+'_sigmaIK_data', skiprows=1)
    rhos = data[:,0]
    phis = rhos * rf**2 * np.pi/4
    sigmas = data[:,2]
    sigmas_std = data[:,3]
    marker = '.'
    cutoff = 1.02
    indices = phis < cutoff
    ax.errorbar(phis[indices], sigmas[indices] * rf, yerr=sigmas_std[indices] * rf, label=fr"$r_f={rf}$", marker=marker, ls='', ms=10, c=color)
    # try:
    #     popt, pcov, infodict, mesg, ier = optimize.curve_fit(sigma_fit, rhos, sigmas, sigma=sigmas_std, absolute_sigma=True, p0=(0.1, 4, 0.001, 2), full_output=True)
    # except RuntimeError:
    #     print(f"Could not fit to SigmaIK of {test_name}")
    # else:
    #     chisq = np.sum(np.square(infodict['fvec']))
    #     dof = len(rhos) - 4
    #     print("SigmaIK " + test_name)
    #     print(popt)
    #     print(pcov)
    #     print(f"chi^2={chisq}/dof={dof}\n")
    #     ax.plot(rho_dense, sigma_fit(rho_dense, *popt), c=color)

    rhos_dense = np.linspace(0, cutoff, 1000)
    phis_dense = rhos_dense * rf**2 * np.pi/4
    pIK_scaled1 = pIK1(rhos_dense, rf, v0, lamb) * rf
    pIK_scaled2 = pIK2(rhos_dense, rf, v0, lamb) * rf
    pIK_scaled3 = pIK3(rhos_dense, rf, v0, lamb) * rf
    # pIK_scaled4 = pIK4(rhos_dense, rf, v0, lamb)

    # pA_scaled1 = pA1(rhos_dense, rf, v0, lamb, Dr)

    plt.plot(phis_dense, pIK_scaled1, label=r"$p_{IK,1}$", color="cyan")
    plt.plot(phis_dense, pIK_scaled2, label=r"$p_{IK,2}$", color="skyblue")
    plt.plot(phis_dense, pIK_scaled3, label=r"$p_{IK,3}$", color="purple")
    # plt.plot(rhos_dense, pIK_scaled1, label=r"$p_{IK}(\rho r_f^2,v_0)/r_f$", color="cyan")
    # plt.plot(rhos_dense, pIK_scaled2, label=r"$\frac{v(\rho)}{v_0}p_{IK}(\rho r_f^2,v_0)/r_f$", color="skyblue")
    # plt.plot(rhos_dense, pIK_scaled3, label=r"$p(\rho,v_{\infty})$")
    # plt.plot(rhos_dense, pIK_scaled4, label=r"$\frac{v(\rho)}{v_{\infty}}p_{IK}(\rho r_f^2,v_{\infty})/r_f$", color="blue")

    # color = f'C{i//num_params+1}'
    # data = np.loadtxt(test_name+'_sigmaA_data', skiprows=1)
    # rhos = data[:,0]
    # Pas = data[:,2]
    # Pas_std = data[:,3]
    # indices = rhos>0
    # ax.errorbar(rhos[indices], Pas[indices], yerr=Pas_std, label=r"$p_A$", marker=marker, ls='', ms=10, c=color)
    # try:
    #     popt, pcov, infodict, mesg, ier = optimize.curve_fit(Pa_fit, rhos, Pas, sigma=Pas_std, absolute_sigma=True, p0=(1.1, 0.21, 3, 0.8), full_output=True)
    # except RuntimeError:
    #     print(f"Could not fit to SigmaA of {test_name}")
    # else:
    #     chisq = np.sum(np.square(infodict['fvec']))
    #     dof = len(rhos) - 4
    #     print("SigmaA " + test_name)
    #     print(popt)
    #     print(pcov)
    #     print(f"chi^2={chisq}/dof={dof}\n")
    #     ax.plot(rho_dense, Pa_fit(rho_dense, *popt), c=color)

    # plt.plot(rhos_dense, pA_scaled1, label=r"$\rho v(\rho)v^*(\rho)/2D_r$", c=color)

    i += num_params

ax.set_xlim(left=0)
# ax.set_ylim(bottom=0)
ax.set_title(fr"Direct pressure")
ax.legend(loc=(1.05,0.2))
plt.savefig('pfap_qsap_pressure_scaling_measured.png',  dpi=300, bbox_inches='tight')
