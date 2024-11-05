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

param_label = r"$\ell_p/r_f$"
test_name = sys.argv[1]

# test_name_slab = sys.argv[2]
# file_slab = test_name_slab + '_density_pressure'

# test_name_solids = sys.argv[3]
# file_solids = test_name_solids + '_histo_phase_diagram'
# pf_rmax1_name = sys.argv[1]

pf_binodal_test = sys.argv[2]

file = test_name + '_histo_phase_diagram'
pf_file = pf_binodal_test + '_histo_phase_diagram2'
# pf_rmax1_file = pf_rmax1_name + '_histo_phase_diagram'

if os.path.exists(pf_file):
    data = np.loadtxt(pf_file, skiprows=1)
    skiphigh = 0
    skiplow = 3
    Pes_pfs = data[skiphigh:-skiplow, 0]
    pf_gases_density = data[skiphigh:-skiplow, 1]
    pf_liquids_density = data[skiphigh:-skiplow, 2]
    rf_pfs = 5/0.5/Pes_pfs

theory_type = "alternate_"

theory_stable = f"pfap_qsap_local_{theory_type}vlow_stable"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.5, facecolor="C0", label="Stable coexistence")
# ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='blue', label="Local theory: stable")
# ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='blue')

theory_stable = f"pfap_qsap_local_{theory_type}vlow_metastable_liquid"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="C1", label="Metastable coexistence")

theory_stable = f"pfap_qsap_local_{theory_type}vlow_metastable_solid"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="C1")


theory_stable = f"pfap_qsap_local_{theory_type}vlow_solid"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory_2 = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory_2, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.5, facecolor="C0")
# ax.plot(area_fraction(gas_theory, rf_theory), 5/rf_theory, color='blue', label="Local theory: stable")
# ax.plot(area_fraction(liquid_theory, rf_theory), 5/rf_theory, color='blue')

theory_stable = f"pfap_qsap_local_{theory_type}vlow_liquid"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory_1 = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory_1, rf_theory), alpha=0.5, facecolor="C0")

# spinodal
theory_stable = f"pfap_qsap_spinodal_onebranch"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="C2", label="Spinodal")

theory_stable = f"pfap_qsap_spinodal_twobranches"
data = np.loadtxt(theory_stable, skiprows=1)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
liquid2_theory = data[:,3]
solid_theory = data[:,4]
ax.fill_betweenx(5/rf_theory, area_fraction(gas_theory, rf_theory), area_fraction(liquid_theory, rf_theory), alpha=0.3, facecolor="C2")
ax.fill_betweenx(5/rf_theory, area_fraction(liquid2_theory, rf_theory), area_fraction(solid_theory, rf_theory), alpha=0.3, facecolor="C2")


if os.path.exists(file):
    # fig, ax = plt.subplots(figsize = (6,6))
    # ax.set_ylabel(param_label)
    # ax.set_xlabel('Density')

    # data = np.loadtxt(file, skiprows=1)
    # rho_gases = data[:,1]
    # rho_liquids = data[:,2]
    # Pes = data[:,0]
    # rfs = 5/1/Pes
    # rfs_dense = np.linspace(0.035, 0.30, 100)
    # rfs_dense_qsap_display = np.linspace(0.035, 0.16, 100)

    # Pes_dense = 5/1/rfs_dense

    alpha = 0.8

    # # QS phase diagrams
    # qs_gases_areafrac = area_fraction(7.4, rfs_dense)
    # qs_liquids_areafrac = area_fraction(51.2, rfs_dense_qsap_display)
    # # ax.plot(qs_gases_areafrac, Pes_dense, color='grey', alpha=alpha, label="QSAPs binodal")
    # # ax.plot(qs_liquids_areafrac, 5/rfs_dense_qsap_display, color='grey', alpha=alpha)

    # # PF phase diagram
    # pf_gases_areafrac = pf_gases_density * np.pi/4
    # pf_liquids_areafrac = pf_liquids_density * np.pi/4
    # ax.plot(pf_gases_areafrac, Pes_pfs, color='grey', marker='v', ls='', alpha=alpha, label="PFAPs binodal")
    # ax.plot(pf_liquids_areafrac, Pes_pfs, color='grey', marker='v', ls='', alpha=alpha)
    # ax.plot(init, 5/rfs_dense, label="Initial density", color='purple', marker='', ls='--')    

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

    # # Plot binodal from histogram
    # area_fraction_gas = area_fraction(rho_gases, rfs)
    # area_fraction_liquid = area_fraction(rho_liquids, rfs)
    # # ax.plot(area_fraction_gas[:4], Pes[:4], color='C0', ls='', marker='.', ms=12, label="PF+QS metastable")
    # # ax.plot(area_fraction_liquid[:4], Pes[:4], color='C0', ls='', marker='.', ms=12)
    # skip = 3
    # ax.plot(area_fraction_gas[skip:], Pes[skip:], color='C0', ls='--', marker='o')
    # ax.plot(area_fraction_liquid[skip:], Pes[skip:], color='C0', ls='--', marker='o')
    # ax.fill_betweenx(Pes[skip:], area_fraction_gas[skip:], area_fraction_liquid[skip:], alpha=0.4, color="C0")

    # # Plot binodal from slab
    # # data_slab = np.loadtxt(file_slab, skiprows=1)
    # # rho_gases = data_slab[:,1]
    # # rho_liquids = data_slab[:,2]
    # # rfs = data_slab[:,0]
    # # area_fraction_gas = area_fraction(rho_gases, rfs)
    # # area_fraction_liquid = area_fraction(rho_liquids, rfs)
    # # ax.plot(area_fraction_gas, Pes, color='C1', ls='', marker='+', ms=12, label="PF+QS slab")
    # # ax.plot(area_fraction_liquid, Pes, color='C1', ls='', marker='+', ms=12)

    # # data_solid = np.loadtxt(file_solids, skiprows=1)
    # # rho_gases = data_solid[:,1]
    # # rho_liquids = data_solid[:,2]
    # # Pes = data_solid[:,0]
    # # rfs = 5/1/Pes
    # # area_fraction_gas = area_fraction(rho_gases, rfs)
    # # area_fraction_liquid = area_fraction(rho_liquids, rfs)
    # # ax.plot(area_fraction_gas[2], Pes[2], color='C0', ls='', marker='x', ms=12)
    # # ax.plot(area_fraction_liquid[2], Pes[2], color='C0', ls='', marker='x', ms=12)

    # # ax.plot(area_fraction_gas[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)
    # # ax.plot(area_fraction_liquid[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)

    # # 2 branches
    # rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    # left_side_rhos1 = np.array([7.852, 7.382, 7.119, 6.880])
    # ax.plot(area_fraction(left_side_rhos1, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C0')

    # rfs_upper = np.array([0.06, 0.05, 0.04, 0.035])[::-1]
    # left_side_rhos2 = np.array([51.597, 51.220, 51.159, 50.430])[::-1]
    # print(rfs_upper)
    # print(left_side_rhos2)
    # ax.plot(area_fraction(left_side_rhos2, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C0')

    # ax.fill_betweenx(5/rfs_upper, area_fraction(left_side_rhos1,rfs_upper), area_fraction(left_side_rhos2,rfs_upper), alpha=0.4, color="C0")

    # rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    # right_side_rhos1  = np.array([109.035, 91.571, 69.379,  54.783])
    # ax.plot(area_fraction(right_side_rhos1, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C0')

    # rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    # right_side_rhos2  = np.array([939.115, 720.674, 460.540, 319.813])
    # ax.plot(area_fraction(right_side_rhos2, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C0')

    # ax.fill_betweenx(5/rfs_upper, area_fraction(right_side_rhos1,rfs_upper), area_fraction(right_side_rhos2,rfs_upper), alpha=0.4, color="C0")

    # # metastable
    # rfs_upper = np.array([0.07, 0.08])
    # metastable_rhos1  = np.array([6.550,  6.222])
    # ax.plot(area_fraction(metastable_rhos1, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C1')

    # rfs_upper = np.array([0.07, 0.08])
    # metastable_rhos2  = np.array([53.420, 54.689])
    # ax.plot(area_fraction(metastable_rhos2, rfs_upper), 5/rfs_upper, ls='--', marker='o', color='C1')

    # ax.fill_betweenx(5/rfs_upper, area_fraction(metastable_rhos1,rfs_upper), area_fraction(metastable_rhos2,rfs_upper), alpha=0.4, color="C1")
    
    # 120.5, 0.035; 107,0.04; 80, 0.05; 53, 0.06: homogeneous
    # rfs_upper = np.array([0.035, 0.04, 0.06])
    # right_side_rhos  = np.array([120.5, 107.0, 53.0])
    # ax.plot(area_fraction(right_side_rhos, rfs_upper), 5/rfs_upper, ls='', marker='o', color='C2', ms=12, label="Reentrant homogeneous phase")

    ax.legend(loc=(1.05,0.4))
    ax.set_xlim(left=0, right=1.2)
    ax.set_ylim(bottom=0, top=150)
    # ax.set_title("PF+QS numerics")
    plt.savefig(f'pfap_qsap_largeeps_phase_diagram_spinodal.png', dpi=300, bbox_inches='tight')
    plt.close()

