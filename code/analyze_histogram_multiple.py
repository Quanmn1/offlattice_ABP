import numpy as np
from scipy import optimize, stats, signal
from scipy.optimize import OptimizeWarning
from matplotlib import pyplot as plt
import os
import sys
import math
from helper import *
import warnings
warnings.simplefilter("error", OptimizeWarning)

"""
Take the last histogram of multiple simulations and add them together.
Used for condensed phase, where histogram will not get averaged out because everything is stationary.
"""

def profile(x, rho_gas, rho_liquid, liquid_ratio, width_gas, width_liquid):
    return (1-liquid_ratio) * stats.norm.pdf(x, rho_gas, width_gas) + liquid_ratio * stats.norm.pdf(x, rho_liquid, width_liquid)

def profile_single_peak(x, rho, width, amplitude):
    return amplitude * stats.norm.pdf(x, rho, width)

def fit_gauss(densities_coarse, histogram_mult, histogram_mult_std, params, rho0, width0):
    popt, pcov, infodict, mesg, ier = optimize.curve_fit(profile_single_peak, densities_coarse, histogram_mult,
                                                         p0=(rho0, width0, 0.5), sigma=histogram_mult_std, absolute_sigma=True,
                                                         full_output=True)
    chisquare = np.sum(np.square(infodict['fvec']))
    dof = len(densities_coarse) - 2
    return popt, pcov, chisquare, dof

def fit_gausses(densities_coarse, histogram_mult, histogram_mult_std, params, rho_gas0=False, rho_liquid0=False):
    if not rho_gas0:
        rho_gas0 = params['rho'] * 0.6
    if not rho_liquid0:
        rho_liquid0 = params['rho'] * 1.6
    liquid_ratio0 = 0.5
    width_gas0 = params['rho'] * 0.1
    width_liquid0 = params['rho'] * 0.1
    # rho_mid0 = params['rho']
    # width_mid0 = params['rho'] * 0.1
    # mid_ratio0 = 0.05

    # histogram_std_fit, densities_fit, histogram_fit = nonzero(histogram_mult_std, densities_coarse, histogram_mult)
    popt, pcov, infodict, mesg, ier = optimize.curve_fit(profile, densities_coarse, histogram_mult,
                                                         p0=(rho_gas0,rho_liquid0,liquid_ratio0,width_gas0,width_liquid0),
                                                         full_output=True)
    chisquare = np.sum(np.square(infodict['fvec']))
    dof = len(densities_coarse) - 5
    return popt, pcov, chisquare, dof

def process_histograms(histograms):
    """
    Take average and standard deviation of multiple normalized histograms.
    """
    histogram_mults = np.array(histograms)
    histogram_mult = np.average(histogram_mults, axis=0)
    histogram_mult_std = np.std(histogram_mults, axis=0)
    return histogram_mult, histogram_mult_std

def analyze_histogram(vars, fit='gauss'):
    number = len(vars)
    param_label = r"$l_p/r_f$"
    num_combined = 1
    # vars = np.linspace(start, end, number)
    rho_gases = np.zeros(number)
    rho_liquids = np.zeros(number)
    rho_gases_std = np.full(number, np.nan)
    rho_liquids_std = np.full(number, np.nan)

    param_file = "pfqs_condensation_test10_0.04_param"
    try:
        params, N = read_param(param_file)
    except FileNotFoundError:
        print("File not found")
        print(var)

    Lx = params['Lx']
    Ly = params['Ly']
    Dr = params['Dr']
    v = params['v']
    density_box_size = params['density_box_size']
    box_size = params['r_max_qsap']
    density_box_raw_size = density_box_size * box_size
    density_box_area = density_box_raw_size**2
    number_of_boxes_x = Lx // density_box_raw_size
    number_of_boxes_y = Ly // density_box_raw_size
    threshold_density = 25

    histogram_folder = 'pfqs_condensation_test10_histogram_folder'
    if not os.path.isdir(histogram_folder):
        os.mkdir(histogram_folder)

    files = {}
    simulation_indices = [1,2,3,4,5,6]
    var_formatting = 2
    for var in vars:
        files[var] = [f'pfqs_condensation_test10.{seed}_{var:.{var_formatting}f}_histogram' for seed in simulation_indices]

    # files = {
    #     0.04: ['pfqs_condensation_test10.1_0.04_histogram','pfqs_condensation_test10.2_0.04_histogram','pfqs_condensation_test10.3_0.04_histogram'],
    #     0.08: ['pfqs_condensation_test10.1_0.08_histogram','pfqs_condensation_test10.2_0.08_histogram','pfqs_condensation_test10.3_0.08_histogram'],
    #     0.12: ['pfqs_condensation_test10.1_0.12_histogram','pfqs_condensation_test10.2_0.12_histogram','pfqs_condensation_test10.3_0.12_histogram'],
    #     0.13: ['pfqs_condensation_test10.1_0.13_histogram','pfqs_condensation_test10.2_0.13_histogram','pfqs_condensation_test10.3_0.13_histogram'],
    #     0.14: ['pfqs_condensation_test10.1_0.14_histogram','pfqs_condensation_test10.2_0.14_histogram','pfqs_condensation_test10.3_0.14_histogram'],
    #     0.15: ['pfqs_condensation_test10.1_0.15_histogram_coarsened','pfqs_condensation_test10.2_0.15_histogram_coarsened','pfqs_condensation_test10.3_0.15_histogram_coarsened'],
    # }

    for (test_num, var) in enumerate(vars):
        file_names = files[var]
        box_size_pfap = var
        if box_size_pfap > 0:
            max_number = math.ceil(density_box_area * 3 / box_size_pfap**2)
        else:
            max_number = params['rho_m'] * density_box_area * 10
        densities = np.arange(0, max_number+1, 1) / density_box_area
        histogram = np.zeros((max_number+1,len(file_names)))
        binwidth = num_combined / density_box_area
        densities_coarse = coarsen(densities, num_combined)

        for (ind, file_name) in enumerate(file_names):
            with open(file_name, 'r') as f:
                histogram_text = f.read().split("\n\n")[-2] # take the histogram of the last snapshots
                histogram_lines = histogram_text.split("\n") # first elt: time. following elements: "density count"
                for i in range(max_number+1):
                    histogram[i, ind] = float(histogram_lines[i+1].split()[1])
        
        histogram_avg = np.mean(histogram, axis=-1)
        histogram_std = np.std(histogram, axis=-1)/np.sqrt(2)

        fig, ax = plt.subplots(figsize = (6,6))
        ax.set_ylabel('Probability')
        ax.set_xlabel('Density')

        indices = (densities < 65) * (densities > 0)
        # here is the normalization factor for histogram, disregarding bin rho=0. num_combined should be 1.
        norm = np.sum(histogram_avg[indices]) * binwidth
        histogram_plot = histogram_avg[indices] / norm
        densities_plot = densities[indices]
        histogram_std_plot = histogram_std[indices] / norm

        ax.errorbar(densities_plot, 
                    histogram_plot, 
                    yerr=histogram_std_plot, color='C0', label="Initialization-averaged histogram")
        
        histogram_fit, histogram_std_fit, densities_fit = nonzero(histogram_plot, histogram_std_plot, densities_plot)
        
        if fit == 'gauss':
            # fit and plot the profile
            try:
                divider = np.searchsorted(densities_fit, threshold_density)
                gas = densities_fit[np.argmax(histogram_fit[:divider])]
                liquid = densities_fit[divider+np.argmax(histogram_fit[divider:])]
                max_density = densities_fit[-1]
                width = (max_density - liquid)*0.6
                # width = 10
                # liquid = 44
                # width = 7

                # gas_indices = (densities_fit >= gas - width) & (densities_fit <= gas + width)
                liquid_indices = (densities_fit >= liquid - width) & (densities_fit <= liquid + width)
                liquid_popt, liquid_pcov, chisquare, dof = fit_gauss(densities_fit[liquid_indices], histogram_fit[liquid_indices], histogram_std_fit[liquid_indices], params, liquid, width)
                # Gas density: know to be 0
                # gas_popt, gas_pcov, chisquare, dof = fit_gauss(densities_fit[gas_indices], histogram_fit[gas_indices], histogram_std_fit[gas_indices], params, gas, width)
                # print(chisquare)
                # print(dof)
                
                # if gas_popt[0] < 0:
                #     raise ValueError("Negative gas density!")
            except (RuntimeError, OptimizeWarning, ValueError) as e:
                print(f"could not fit gauss: {e=} for r={var}. Taking maximum values instead.")
                divider = np.searchsorted(densities_fit, threshold_density)
                rho_gases[test_num] = densities_fit[np.argmax(histogram_fit[:divider])]
                rho_liquids[test_num] = densities_fit[divider+np.argmax(histogram_fit[divider:])]
                ax.set_title(f'$Peaks: \\rho_g = {rho_gases[test_num]:.3f}, \\rho_l = {rho_liquids[test_num]:.3f}$')
                ax.text(0.4,0.9, f'Invalid fit', transform=ax.transAxes)
            else:
                densities_dense = np.linspace(0, 65, 1000)
                ax.plot(densities_dense, profile_single_peak(densities_dense, *liquid_popt), c='C1', label="Gaussian fits")
            # if np.all(np.isfinite(gas_pcov)) and np.all(np.isfinite(liquid_pcov)):
                rho_gases[test_num] = 0
                rho_liquids[test_num] = liquid_popt[0]
                rho_gases_std[test_num] = 0
                rho_liquids_std[test_num] = np.sqrt(np.diag(liquid_pcov))[0]
                ax.set_title(f'$\\rho_g = {rho_gases[test_num]:.3f}\pm {rho_gases_std[test_num]:.3f}, \\rho_l = {rho_liquids[test_num]:.3f}\pm {rho_liquids_std[test_num]:.3f}$')

                            # else:
                            #     rho_gases[test_num] = np.nan
                            #     rho_liquids[test_num] = np.nan
                            #     rho_gases_std[test_num] = np.inf
                            #     rho_liquids_std[test_num] = np.inf
                            #     ax.set_title(f'$t={t:.2f}, \\rho_g = {gas_popt[0]:.3f}, \\rho_l = {liquid_popt[0]:.3f}$')
                            #     ax.text(0.4,0.9, 'Invalid fit', transform=ax.transAxes)
                                # ax.text(0.4,0.8, f'$\chi^2 / dof = {chisquare:.2f}/{dof}$', transform=ax.transAxes)
        elif fit == 'max':
            divider = np.searchsorted(densities_fit, threshold_density)
            rho_gases[test_num] = densities_fit[np.argmax(histogram_fit[:divider])]
            rho_liquids[test_num] = densities_fit[divider+np.argmax(histogram_fit[divider:])]
            ax.set_title(f'$Peaks: \\rho_g = {rho_gases[test_num]:.3f}, \\rho_l = {rho_liquids[test_num]:.3f}$')
                    
        ax.legend()
        ax.grid()
        plt.savefig(histogram_folder + '/' + f'{var}.png', dpi=100, bbox_inches='tight')
        plt.close()


    rho_liquids, rho_gases, vars, rho_gases_std, rho_liquids_std = nonzero(rho_liquids, rho_gases, vars, rho_gases_std, rho_liquids_std)

    # print(rho_gases)
    # print(rho_liquids)

    with open('pfqs_condensation_test10_histo_phase_diagram_frozen', 'w') as f:
        f.write(f"{param_label} \t rho_gas \t rho_liquid \t rho_gas_std \t rho_liquid_std \n")
        for i in range(len(vars)):
            f.write(f"{vars[i]:.3f}\t{rho_gases[i]:.3f}\t{rho_liquids[i]:.3f}\t{rho_gases_std[i]:.3f}\t{rho_liquids_std[i]:.3f}\n")

    return vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label


if __name__ == "__main__":
    vars = np.array(sys.argv[1].split(),dtype=float)
    fit = sys.argv[2]
    vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label = analyze_histogram(vars, fit=fit)

    fig, ax = plt.subplots(figsize = (6,6))
    ax.set_ylabel(param_label)
    ax.set_xlabel('Density')
    ax.errorbar(rho_gases, vars, xerr=rho_gases_std, color='C0', ls='', marker='.')
    ax.errorbar(rho_liquids, vars, xerr=rho_liquids_std, color='C1', ls='', marker='.')
    ax.set_title("Histogram phase diagram")
    plt.savefig('pfqs_condensation_test10_histo_phase_diagram_frozen.png', dpi=300, bbox_inches='tight')
    plt.close()