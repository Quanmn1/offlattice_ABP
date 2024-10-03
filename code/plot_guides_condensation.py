import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

param_label = r"$l_p/r_f$"
test_name = sys.argv[1]

# pf_rmax1_name = sys.argv[1]

pf_binodal_test = sys.argv[2]

file = test_name + '_histo_phase_diagram'
pf_file = pf_binodal_test + '_histo_phase_diagram'
# pf_rmax1_file = pf_rmax1_name + '_histo_phase_diagram'

if os.path.exists(pf_file):
    data = np.loadtxt(pf_file, skiprows=1)
    skiphigh = 1
    skiplow = -3
    Pes_pfs = data[skiphigh:skiplow, 0]
    pf_gases_density = data[skiphigh:skiplow, 1]
    pf_liquids_density = data[skiphigh:skiplow, 2]
    # rf_pfs = 5/0.5/Pes_pfs

# Add lines of QSAP binodals and PFAP crystal densities to the phase diagram
# Improvement: add lines of PFAP binodals. 
# We know the crystall densities are not followed by PFAPs system alone

def area_fraction(density, rf):
    return density * rf**2 * np.pi/4

def scaled_rho(density, rf):
    return density * rf**2

if os.path.exists(file):
    # fig, ax = plt.subplots(figsize = (6,6))
    # ax.set_ylabel(param_label)
    # ax.set_xlabel('Density')

    data = np.loadtxt(file, skiprows=1)
    rho_gases = data[:,1]
    rho_liquids = data[:,2]
    Pes = data[:,0]
    rfs = 5/1/Pes
    rfs_dense = np.linspace(0.04, 0.14, 100)
    Pes_dense = 5/1/rfs_dense
    # # Plot binodal from histogram
    # ax.errorbar(rho_gases, rfs, color='blue', ls='', marker='.', label="Gas (histogram)")
    # ax.errorbar(rho_liquids, rfs, color='orange', ls='', marker='.', label = "Liquid (histogram)")
    
    # # ax.axvline(56.5, label="Liquid QS", color='orange')
    # # ax.axvline(7.4, label="Gas QS", color='blue')
    # # rfs_dense = np.linspace(0.04, 0.15, 100)[::-1]
    # # pf_rho = 2/np.sqrt(3)/rfs_dense**2
    # # cutoff = np.searchsorted(pf_rho, np.max(rho_liquids))
    # # ax.plot(pf_rho[:cutoff], 5/rfs_dense[:cutoff], label="PF crystals", color='black')
    # # ax.scatter(pf_gases_density, rf_pfs, label="PF gas", color='C0')
    # # ax.scatter(pf_liquids_density, rf_pfs, label="PF liquid", color='C1')
    # ax.legend()
    # ax.set_title("Density phase diagram")
    # plt.savefig(test_name + '_annotated_phase_diagram.png', dpi=300, bbox_inches='tight')
    # plt.close()

    fig, ax = plt.subplots(figsize = (6,6))
    ax.set_ylabel(param_label)
    ax.set_xlabel(r'$\phi=\pi \rho r_f^2 / 4$')

    # data = np.loadtxt(pf_rmax1_file, skiprows=1)
    # skip = 3
    # rho_gases_rmax1 = data[skip:, -2]
    # rho_liquids_rmax1 = data[skip:,-1]
    # Pes_pfs_rmax1 = data[skip:,-3]
    # ax.plot(rho_gases_rmax1, Pes_pfs_rmax1, label="Vary Dr", color='C1', marker='+', ls='')
    # ax.plot(rho_liquids_rmax1, Pes_pfs_rmax1, color='C1', marker='+', ls='')

    # Plot binodal from histogram
    area_fraction_gas = area_fraction(rho_gases, rfs)
    area_fraction_liquid = area_fraction(rho_liquids, rfs)
    ax.plot(area_fraction_gas, Pes, color='C0', ls='', marker='.', ms=12, label="PF+QS")
    ax.plot(area_fraction_liquid, Pes, color='C0', ls='', marker='.', ms=12)

    # QS phase diagrams
    qs_gases_areafrac = area_fraction(7.4, rfs_dense)
    qs_liquids_areafrac = area_fraction(56.5, rfs_dense)
    
    # init = area_fraction(19.44, rfs_dense)
    # rho_m = area_fraction(25, rfs_dense)
    ax.plot(qs_gases_areafrac, Pes_dense, label="QS", color='C1')
    ax.plot(qs_liquids_areafrac, Pes_dense, color='C1')

    # PF phase diagram
    # ax.axvline(np.pi/2/np.sqrt(3), label="PF crystals", color='black')
    # ax.plot(rho_m, rfs_dense, label="QS density scale", color='grey')
    # pf_gases_areafrac = area_fraction(pf_gases_density, rf_pfs)
    # pf_liquids_areafrac = area_fraction(pf_liquids_density, rf_pfs)
    pf_gases_areafrac = pf_gases_density * np.pi/4
    pf_liquids_areafrac = pf_liquids_density * np.pi/4
    ax.plot(pf_gases_areafrac, Pes_pfs, label="PF", color='C2', marker='v', ls='')
    ax.plot(pf_liquids_areafrac, Pes_pfs, color='C2', marker='v', ls='')
    # ax.plot(init, 5/rfs_dense, label="Initial density", color='purple', marker='', ls='--')    
    ax.legend()
    ax.set_title(fr"$\epsilon=100$")
    plt.savefig('pfap_qsap_largeeps_phase_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()

