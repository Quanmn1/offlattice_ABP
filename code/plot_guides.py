import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

def area_fraction(density, rf):
    return density * np.pi/4 * rf**2

def ycoord(rf):
    if config == 0:
        lp = 5
        return lp/rf
    elif config == 1:
        return rf
    
def xcoord(rho, rf):
    if config == 0:
        return area_fraction(rho, rf)
    elif config == 1:
        return rho

def scaled_rho(density, rf):
    return density * rf**2

matplotlib.rcParams.update({'font.size': 16})

# whether plot as (phi, lp/rf) or (rho, rf)
config = int(sys.argv[1])

fig, ax = plt.subplots(figsize = (10,6))

test_name = "pfaps_qsaps_test16_harmonic_largescaleeps"

pf_binodal_test = "pfaps_test49_vlowqsap"

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

theory_type = "gradient"
skip = 0

theory_stable = f"pfap_qsap_{theory_type}_stable"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.5, facecolor="C0", label="Stable coexistence")

theory_stable = f"pfap_qsap_{theory_type}_metastable_liquid"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.3, facecolor="C1", label="Metastable coexistence")

theory_stable = f"pfap_qsap_{theory_type}_metastable_solid"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.3, facecolor="C1")

theory_stable = f"pfap_qsap_{theory_type}_solid"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory_2 = data[:,1]
liquid_theory = data[:,2]
ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory_2, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.5, facecolor="C0")

theory_stable = f"pfap_qsap_{theory_type}_liquid"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory_1 = data[:,2]
ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory_1, rf_theory), alpha=0.5, facecolor="C0")

# spinodal
spinodal_theory_type = "alternate_correctpIK"
theory_stable = f"pfap_qsap_spinodal_{spinodal_theory_type}_onebranch"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
# ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.3, facecolor="C2", label="Spinodal")
ax.plot(xcoord(gas_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2", label="Spinodal")
ax.plot(xcoord(liquid_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2")

theory_stable = f"pfap_qsap_spinodal_{spinodal_theory_type}_twobranches"
data = np.loadtxt(theory_stable, skiprows=skip)
rf_theory = data[:,0]
gas_theory = data[:,1]
liquid_theory = data[:,2]
liquid2_theory = data[:,3]
solid_theory = data[:,4]
# ax.fill_betweenx(ycoord(rf_theory), xcoord(gas_theory, rf_theory), xcoord(liquid_theory, rf_theory), alpha=0.3, facecolor="C2")
# ax.fill_betweenx(ycoord(rf_theory), xcoord(liquid2_theory, rf_theory), xcoord(solid_theory, rf_theory), alpha=0.3, facecolor="C2")
ax.plot(xcoord(gas_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2")
ax.plot(xcoord(liquid_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2")
ax.plot(xcoord(solid_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2")
ax.plot(xcoord(liquid2_theory, rf_theory), ycoord(rf_theory), ls="--", color="C2")


if os.path.exists(file):
    data = np.loadtxt(file, skiprows=1)
    rho_gases = data[:,1]
    rho_liquids = data[:,2]
    Pes = data[:,0]
    rfs = 5/1/Pes
    rfs_dense = np.linspace(0.03, 0.32, 100)
    # rfs_dense_qsap_display = np.linspace(0.035, 0.16, 100)

    Pes_dense = 5/1/rfs_dense

    alpha = 0.8

    # # QS phase diagrams
    # qs_gases_areafrac = xcoord(7.4, rfs_dense)
    # qs_liquids_areafrac = xcoord(51.2, rfs_dense_qsap_display)
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


    # data = np.loadtxt(pf_rmax1_file, skiprows=1)
    # skip = 3
    # rho_gases_rmax1 = data[skip:, -2]
    # rho_liquids_rmax1 = data[skip:,-1]
    # Pes_pfs_rmax1 = data[skip:,-3]
    # ax.plot(rho_gases_rmax1, Pes_pfs_rmax1, label="Vary Dr", color='C1', marker='+', ls='')
    # ax.plot(rho_liquids_rmax1, Pes_pfs_rmax1, color='C1', marker='+', ls='')

    # Plot binodal from histogram
    area_fraction_gas = xcoord(rho_gases, rfs)
    area_fraction_liquid = xcoord(rho_liquids, rfs)
    # ax.plot(area_fraction_gas[:4], Pes[:4], color='C0', ls='', marker='.', ms=12, label="PF+QS metastable")
    # ax.plot(area_fraction_liquid[:4], Pes[:4], color='C0', ls='', marker='.', ms=12)
    skip = 3
    ax.plot(area_fraction_gas[skip:], ycoord(5/Pes[skip:]), color='C0', ls='--', marker='o')
    ax.plot(area_fraction_liquid[skip:], ycoord(5/Pes[skip:]), color='C0', ls='--', marker='o')
    # ax.fill_betweenx(Pes[skip:], area_fraction_gas[skip:], area_fraction_liquid[skip:], alpha=0.4, color="C0")

    # # Plot binodal from slab
    # # data_slab = np.loadtxt(file_slab, skiprows=1)
    # # rho_gases = data_slab[:,1]
    # # rho_liquids = data_slab[:,2]
    # # rfs = data_slab[:,0]
    # # area_fraction_gas = xcoord(rho_gases, rfs)
    # # area_fraction_liquid = xcoord(rho_liquids, rfs)
    # # ax.plot(area_fraction_gas, Pes, color='C1', ls='', marker='+', ms=12, label="PF+QS slab")
    # # ax.plot(area_fraction_liquid, Pes, color='C1', ls='', marker='+', ms=12)

    # # data_solid = np.loadtxt(file_solids, skiprows=1)
    # # rho_gases = data_solid[:,1]
    # # rho_liquids = data_solid[:,2]
    # # Pes = data_solid[:,0]
    # # rfs = 5/1/Pes
    # # area_fraction_gas = xcoord(rho_gases, rfs)
    # # area_fraction_liquid = xcoord(rho_liquids, rfs)
    # # ax.plot(area_fraction_gas[2], Pes[2], color='C0', ls='', marker='x', ms=12)
    # # ax.plot(area_fraction_liquid[2], Pes[2], color='C0', ls='', marker='x', ms=12)

    # # ax.plot(area_fraction_gas[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)
    # # ax.plot(area_fraction_liquid[:2], Pes[:2], color='C0', ls='', marker='.', ms=12)

    # 2 branches
    rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    left_side_rhos1 = np.array([7.852, 7.382, 7.119, 6.880])
    ax.plot(xcoord(left_side_rhos1, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C0')

    rfs_upper = np.array([0.06, 0.05, 0.04, 0.035])[::-1]
    left_side_rhos2 = np.array([51.597, 51.220, 51.159, 50.430])[::-1]
    print(rfs_upper)
    print(left_side_rhos2)
    ax.plot(xcoord(left_side_rhos2, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C0')

    # ax.fill_betweenx(ycoord(rfs_upper), xcoord(left_side_rhos1,rfs_upper), xcoord(left_side_rhos2,rfs_upper), alpha=0.4, color="C0")

    rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    right_side_rhos1  = np.array([109.035, 91.571, 69.379,  54.783])
    ax.plot(xcoord(right_side_rhos1, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C0')

    rfs_upper = np.array([0.035, 0.04, 0.05, 0.06])
    right_side_rhos2  = np.array([939.115, 720.674, 460.540, 319.813])
    ax.plot(xcoord(right_side_rhos2, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C0')

    # ax.fill_betweenx(ycoord(rfs_upper), xcoord(right_side_rhos1,rfs_upper), xcoord(right_side_rhos2,rfs_upper), alpha=0.4, color="C0")

    # metastable
    rfs_upper = np.array([0.065, 0.07, 0.08])
    metastable_rhos1  = np.array([6.515, 6.550,  6.222])
    ax.plot(xcoord(metastable_rhos1, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C1')

    rfs_upper = np.array([0.065, 0.07, 0.08])
    metastable_rhos2  = np.array([52.754, 53.420, 54.689])
    ax.plot(xcoord(metastable_rhos2, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C1')

    # ax.fill_betweenx(ycoord(rfs_upper), xcoord(metastable_rhos1,rfs_upper), xcoord(metastable_rhos2,rfs_upper), alpha=0.4, color="C1")

    rfs_upper = np.array([0.065])
    metastable_rhos2  = np.array([52.533])
    ax.plot(xcoord(metastable_rhos2, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C1')

    rfs_upper = np.array([0.065])
    metastable_rhos2  = np.array([266.721])
    ax.plot(xcoord(metastable_rhos2, rfs_upper), ycoord(rfs_upper), ls='--', marker='o', color='C1')


    # ax.fill_betweenx(ycoord(rfs_upper), xcoord(metastable_rhos1,rfs_upper), xcoord(metastable_rhos2,rfs_upper), alpha=0.4, color="C1")
    
    # spinodal:
    # 260, 0.035; 120.5, 0.035; 107,0.04; 80, 0.05; 53, 0.06: homogeneous
    # 60, 0.035-0.075 outside of spinodal. 0.085 and above is within spinodal
    # 50 60 70 80, 0.065 outside of spinodal. 40 within spinodal
    # 55 60, 0.08 outside of spinodal but within binodal. 50 65, 0.08 within spinodal.

    # rfs = np.array([0.035, 0.04, 0.04, 0.05, 0.05, 0.06, 0.065, 0.065, 0.065, 0.065, 0.07, 0.07, 0.08, 0.08])
    # rhos  = np.array([260, 107.0, 220.0, 80.0, 60.0, 60.0, 50.0, 60.0, 70.0, 80.0, 53.4, 60.0, 54.7, 60.0])
    # ax.plot(xcoord(rhos, rfs), 5/rfs, ls='', marker='o', color='C3', label="Outside of spinodal")

    # rfs = np.array([0.035, 0.05, 0.06, 0.065, 0.07, 0.08, 0.08, 0.085])
    # rhos  = np.array([320.0, 170.0, 140.0, 40.0, 110.0, 50.0, 65.0, 60.0])
    # ax.plot(xcoord(rhos, rfs), 5/rfs, ls='', marker='o', color='C2', label="Inside spinodal")

    # ax.legend(loc=(1.05,0.4))
    if config == 0:
        ax.set_xlim(left=0, right=1.2)
        ax.set_ylim(bottom=0, top=150)
        plotname = 'pfap_qsap_largeeps_phase_diagram_gradient.png'
        ax.set_ylabel(r"$\ell_p/r_f$")
        ax.set_xlabel(r'$\phi=\pi\rho r_f^2/4$')
    elif config == 1:
        ax.axvline(x=25,ls='--',color="C3")
        ax.set_xlim(left=0, right=80)
        ax.set_ylim(bottom=0)
        plotname = 'pfap_qsap_largeeps_phase_diagram_gradient_rhorf.png'
        ax.set_ylabel(r"$r_f$")
        ax.set_xlabel(r'$\rho$')

    # ax.set_title("PF+QS numerics")
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
    plt.close()

