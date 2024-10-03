import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

def area_fraction(density, rf):
    return density * rf * rf * np.pi/4

def scaled_rho(density, rf):
    return density * rf**2

matplotlib.rcParams.update({'font.size': 16})

fig, ax = plt.subplots(figsize = (10,6))

param_label = r"$l_p/r_f$"
test_name = sys.argv[1]

# test_name_slab = sys.argv[2]
# file_slab = test_name_slab + '_density_pressure'

# test_name_solids = sys.argv[3]
# file_solids = test_name_solids + '_histo_phase_diagram'
# pf_rmax1_name = sys.argv[1]

pf_binodal_test = sys.argv[2]

file = test_name + '_histo_phase_diagram'
pf_file = pf_binodal_test + '_histo_phase_diagram'
# pf_rmax1_file = pf_rmax1_name + '_histo_phase_diagram'

if os.path.exists(pf_file):
    data = np.loadtxt(pf_file, skiprows=1)
    skiphigh = 1
    skiplow = 2
    Pes_pfs = data[skiphigh:-skiplow, 0]
    pf_gases_density = data[skiphigh:-skiplow, 1]
    pf_liquids_density = data[skiphigh:-skiplow, 2]
    # rf_pfs = 5/0.5/Pes_pfs


# rfs_upper = np.array([0.06, 0.07, 0.07, 0.06])
# left_side_rhos = np.array([10.694419931568554, 9.931359544747822, 117.26031214966136, 92.57068619002868])
# right_side_rhos  = np.array([68.12138768578814,  52.26277572949672, 304.8448007099825, 423.5427770896261])
# ax.plot(area_fraction(left_side_rhos, rfs_upper), 5/rfs_upper, color='orange', label="Local theory: metastable")
# ax.plot(area_fraction(right_side_rhos, rfs_upper), 5/rfs_upper, color='purple', label="Local theory: metastable")


# theory_stable = "pfap_qsap_local_stable"
# data = np.loadtxt(theory_stable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory = data[:,1]
# liquid_theory = data[:,2]
# ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="blue", label="Stable coexistence")
# # ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='blue', label="Local theory: stable")
# # ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='blue')

# theory_stable = "pfap_qsap_local_metastable_liquid"
# data = np.loadtxt(theory_stable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory = data[:,1]
# liquid_theory = data[:,2]
# ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="orange", label="Metastable coexistence")

# theory_stable = "pfap_qsap_local_metastable_solid"
# data = np.loadtxt(theory_stable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory = data[:,1]
# liquid_theory = data[:,2]
# ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="purple", label="Metastable coexistence")


# theory_stable = "pfap_qsap_local_solid"
# data = np.loadtxt(theory_stable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory_2 = data[:,1]
# liquid_theory = data[:,2]
# ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory_2, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="blue")
# # ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='blue', label="Local theory: stable")
# # ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='blue')

# theory_stable = "pfap_qsap_local_liquid"
# data = np.loadtxt(theory_stable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory = data[:,1]
# liquid_theory_1 = data[:,2]
# ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory_1, rf_theory), alpha=0.3, facecolor="blue")
# # ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='blue', label="Local theory: stable")
# # ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='blue')

# ax.fill_betweenx(5/rf_theory, area_fraction(liquid_theory_1, rf_theory), area_fraction(gas_theory_2, rf_theory), alpha=0.3, facecolor="grey", label="Reentrant homogeneous")


# theory_metastable = "pfap_qsap_local_metastable"
# data = np.loadtxt(theory_metastable, skiprows=1)
# rf_theory = data[:,0]
# gas_theory = data[:,1]
# liquid_theory = data[:,2]
# ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='brown', label="Local theory: metastable")
# ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='brown')


if os.path.exists(file):
    # fig, ax = plt.subplots(figsize = (6,6))
    # ax.set_ylabel(param_label)
    # ax.set_xlabel('Density')

    data = np.loadtxt(file, skiprows=1)
    rho_gases = data[:,1]
    rho_liquids = data[:,2]
    Pes = data[:,0]
    rfs = 5/1/Pes
    rfs_dense = np.linspace(0.035, 0.14, 100)
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

    ax.set_ylabel(param_label)
    ax.set_xlabel(r'$\phi=\pi \rho r_f^2 /4$')

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
    # ax.plot(area_fraction_gas[:4], Pes[:4], color='C0', ls='', marker='.', ms=12, label="PF+QS metastable")
    # ax.plot(area_fraction_liquid[:4], Pes[:4], color='C0', ls='', marker='.', ms=12)
    skip = 3
    ax.plot(area_fraction_gas[skip:], Pes[skip:], color='C0', ls='--', marker='x', ms=12, label="PF+QS stable")
    ax.plot(area_fraction_liquid[skip:], Pes[skip:], color='C0', ls='--', marker='x', ms=12)

    # Plot binodal from slab
    # data_slab = np.loadtxt(file_slab, skiprows=1)
    # rho_gases = data_slab[:,1]
    # rho_liquids = data_slab[:,2]
    # rfs = data_slab[:,0]
    # area_fraction_gas = area_fraction(rho_gases, rfs)
    # area_fraction_liquid = area_fraction(rho_liquids, rfs)
    # ax.plot(area_fraction_gas, Pes, color='C1', ls='', marker='+', ms=12, label="PF+QS slab")
    # ax.plot(area_fraction_liquid, Pes, color='C1', ls='', marker='+', ms=12)

    # data_solid = np.loadtxt(file_solids, skiprows=1)
    # rho_gases = data_solid[:,1]
    # rho_liquids = data_solid[:,2]
    # Pes = data_solid[:,0]
    # rfs = 5/1/Pes
    # area_fraction_gas = area_fraction(rho_gases, rfs)
    # area_fraction_liquid = area_fraction(rho_liquids, rfs)
    # ax.plot(area_fraction_gas[2], Pes[2], color='C0', ls='', marker='x', ms=12)
    # ax.plot(area_fraction_liquid[2], Pes[2], color='C0', ls='', marker='x', ms=12)

    # ax.plot(area_fraction_gas[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)
    # ax.plot(area_fraction_liquid[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)

    # 2 branches
    rfs_upper = np.array([0.035, 0.04, 0.05, 0.06, 0.06, 0.05, 0.04, 0.035])
    left_side_rhos = np.array([7.852, 7.382, 7.119, 6.880, 51.597, 51.220, 51.159, 50.430])
    ax.plot(area_fraction(left_side_rhos, rfs_upper), 5/rfs_upper, ls='--', marker='x', color='C0', ms=12)

    rfs_upper = np.array([0.05, 0.06, 0.06, 0.05, 0.04])
    right_side_rhos  = np.array([69.379,  54.783, 319.813, 460.540, 740])
    ax.plot(area_fraction(right_side_rhos, rfs_upper), 5/rfs_upper, ls='--', marker='x', color='C0', ms=12)

    # metastable
    rfs_upper = np.array([0.07, 0.08, 0.08, 0.07])
    right_side_rhos  = np.array([6.550,  6.222, 54.689, 53.420])
    ax.plot(area_fraction(right_side_rhos, rfs_upper), 5/rfs_upper, ls='--', marker='.', color='C1', ms=12, label="PF+QS metastable")

    # 120.5, 0.035; 107,0.04; 80, 0.05; 53, 0.06: homogeneous
    # rfs_upper = np.array([0.035, 0.04, 0.06])
    # right_side_rhos  = np.array([120.5, 107.0, 53.0])
    # ax.plot(area_fraction(right_side_rhos, rfs_upper), 5/rfs_upper, ls='', marker='o', color='C2', ms=12, label="Reentrant homogeneous phase")

    # QS phase diagrams
    qs_gases_areafrac = area_fraction(7.4, rfs_dense)
    qs_liquids_areafrac = area_fraction(56.5, rfs_dense)
    # init = area_fraction(19.44, rfs_dense)
    # rho_m = area_fraction(25, rfs_dense)
    ax.plot(qs_gases_areafrac, Pes_dense, label="QS", color='grey', alpha=0.5)
    ax.plot(qs_liquids_areafrac, Pes_dense, color='grey', alpha=0.5)

    # PF phase diagram
    # ax.axvline(np.pi/2/np.sqrt(3), label="PF crystals", color='black')
    # ax.plot(rho_m, rfs_dense, label="QS density scale", color='grey')
    # pf_gases_areafrac = area_fraction(pf_gases_density, rf_pfs)
    # pf_liquids_areafrac = area_fraction(pf_liquids_density, rf_pfs)
    rfs = 5/Pes_pfs
    pf_gases_areafrac = pf_gases_density * np.pi/4
    pf_liquids_areafrac = pf_liquids_density * np.pi/4
    ax.plot(pf_gases_areafrac, Pes_pfs, label="PF", color='grey', marker='v', ls='', alpha=0.5)
    ax.plot(pf_liquids_areafrac, Pes_pfs, color='grey', marker='v', ls='', alpha=0.5)
    # ax.plot(init, 5/rfs_dense, label="Initial density", color='purple', marker='', ls='--')    
    ax.legend(loc=(1.05,0.25))
    ax.set_xlim(left=0)
    ax.set_title("PF+QS numerics")
    plt.savefig('pfap_qsap_largeeps_phase_diagram_numerics.png', dpi=300, bbox_inches='tight')
    plt.close()

