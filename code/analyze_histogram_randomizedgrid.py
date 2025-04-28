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
Take the last snapshot of one simulation. Randomize the position of the grid N times, calculate the histogram, add them together.
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
    Take the mean and the standard deviation of the mean of multiple histograms.
    """
    histogram_mult = np.average(histograms, axis=0)
    histogram_mult_std = np.std(histograms, axis=0) / np.sqrt(histograms.shape[0])
    return histogram_mult, histogram_mult_std

def histogram_randomizedgrid(snapshot, density_box_size, density_box_origin, sizes, max_number):
    """
    For the grid specified by density_box_size and density_box_origin, take the histogram of the snapshot.
    """
    origin_x, origin_y = density_box_origin
    Lx, Ly = sizes
    number_box_x = int(Lx / density_box_size)
    number_box_y = int(Ly / density_box_size)
    box = np.zeros((number_box_x, number_box_y), dtype=int)
    for i in range(len(snapshot)):
        x = snapshot[i, 0]
        y = snapshot[i, 1]
        # calculate the adjusted distance to the new origin
        x_adjusted = (x - origin_x) % Lx
        y_adjusted = (y - origin_y) % Ly
        # determine the box that particle i is in
        box_x = int(x_adjusted / density_box_size)
        box_y = int(y_adjusted / density_box_size)
        box[box_x, box_y] += 1
    
    # add to the histogram of the density in the box
    histogram = np.zeros(max_number)
    for i in range(number_box_x):
        for j in range(number_box_y):
            histogram[box[i, j]] += 1

    return histogram


def analyze_histogram(test_name, mode, vars, pad, fit='gauss'):
    rng = np.random.default_rng()
    number = len(vars)
    param_label = r"$l_p/r_f$"
    num_combined = 1
    # vars = np.linspace(start, end, number)
    rho_gases = np.zeros(number)
    rho_liquids = np.zeros(number)
    rho_gases_std = np.full(number, np.nan)
    rho_liquids_std = np.full(number, np.nan)

    for (test_num, var) in enumerate(vars):
        name = test_name + f'_{var:.{pad}f}'
        param_file = name + '_param'
        snapshot_file = f'{name}_video/{name}_last_state'

        # Read parameters into dict params, with string values
        try:
            params, N = read_param(param_file)
        except FileNotFoundError:
            print("File not found")
            print(var)
            continue

        # N = params['N']
        Lx = params['Lx']
        Ly = params['Ly']
        sizes = (Lx, Ly)
        Dr = params['Dr']
        v = params['v']
        # density_box_size = params['density_box_size']
        density_box_size = 1

        try:
            if mode == "pfap":
                box_size = np.ceil(params['r_max_pfap'])
                rmax = params['r_max_pfap']
            elif mode == "qsap":
                box_size = params['r_max_qsap']
            elif mode == "pfqs":
                box_size = params['r_max_qsap']
                box_size_pfap = params['r_max_pfap']
        except KeyError:
            box_size = params['r_max'] # old version of code
        tf = params['final_time']
        try:
            N = params['N']
        except KeyError:
            lf = params['liquid_fraction']
            rsmall = params['rho_small']
            rlarge = params['rho_large']
            N = int(Lx*Ly * rsmall*(1-lf)) + int(Lx*Ly*rlarge*lf)
        density_box_raw_size = density_box_size * box_size
        density_box_area = density_box_raw_size**2
        number_of_boxes_x = Lx // density_box_raw_size
        number_of_boxes_y = Ly // density_box_raw_size
        params['rho'] = N/Lx/Ly
        rho_m = params['rho_m']
        max_density = rho_m * 3
        threshold_density = rho_m
        max_number = int(max_density*density_box_area)

        densities = np.arange(0, max_number, 1) / density_box_area
        binwidth = num_combined / density_box_area

        histogram_folder = name+'_histogram_folder'
        if not os.path.isdir(histogram_folder):
            os.mkdir(histogram_folder)

        # Read the snapshot
        snapshot = np.loadtxt(snapshot_file, skiprows=1)

        # Number of random histograms
        num_histograms = 1600

        histograms = np.zeros((num_histograms, max_number))

        for i in range(num_histograms):
            # Generate random grid origins
            density_box_origin = rng.random(size=2) * density_box_size
            histograms[i, :] = histogram_randomizedgrid(snapshot, density_box_raw_size, density_box_origin, sizes, max_number)

        # Take the average and the std
        histogram_avg, histogram_std = process_histograms(histograms)

        # print(densities)
        # print(histogram_avg)
        # print(histogram_std)

        fig, ax = plt.subplots(figsize = (6,6))
        ax.set_ylabel('Probability')
        ax.set_xlabel('Density')

        indices = (densities < 65) * (densities > 0)
        # here is the normalization factor for histogram, disregarding bin rho=0. num_combined should be 1.
        norm = np.sum(histogram_avg[indices]) * binwidth
        histogram_plot = histogram_avg[indices] / norm
        densities_plot = densities[indices]
        histogram_std_plot = histogram_std[indices] / norm

        histogram_fit, histogram_std_fit, densities_fit = nonzero(histogram_plot, histogram_std_plot, densities_plot)

        ax.errorbar(densities_fit, 
                    histogram_fit, 
                    yerr=histogram_std_fit, color='C0', label="Grid-averaged histogram")
        
        
        if fit == 'gauss':
            # fit and plot the profile
            try:
                divider = np.searchsorted(densities_fit, threshold_density)
                gas = densities_fit[np.argmax(histogram_fit[:divider])]
                liquid = densities_fit[divider+np.argmax(histogram_fit[divider:])]
                max_density = densities_fit[-1]
                width = (max_density - liquid)*0.5
                # gas_indices = (densities_fit >= gas - width) & (densities_fit <= gas + width)
                liquid_indices = (densities_fit >= liquid - width) & (densities_fit <= liquid + width)
                liquid_popt, liquid_pcov, chisquare, dof = fit_gauss(densities_fit[liquid_indices], histogram_fit[liquid_indices], histogram_std_fit[liquid_indices], params, liquid, width)
                # Gas density: know to be 0
                # gas_popt, gas_pcov, chisquare, dof = fit_gauss(densities_fit[gas_indices], histogram_fit[gas_indices], histogram_std_fit[gas_indices], params, gas, width)
                
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
                densities_dense = np.linspace(0, max_density, 1000)
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

    with open(f'{test_name}_histo_phase_diagram_frozen', 'w') as f:
        f.write(f"{param_label} \t rho_gas \t rho_liquid \t rho_gas_std \t rho_liquid_std \n")
        for i in range(len(vars)):
            f.write(f"{vars[i]:.3f}\t{rho_gases[i]:.3f}\t{rho_liquids[i]:.3f}\t{rho_gases_std[i]:.3f}\t{rho_liquids_std[i]:.3f}\n")

    return vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label


if __name__ == "__main__":
    """
    usage: python analyze_histogram_randomizedgrid.py test_name "r1 r2 ... rN" pad
    """
    test_name = sys.argv[1]
    vars = np.array(sys.argv[2].split(),dtype=float)
    pad = int(sys.argv[3])
    mode = 'pfqs'
    fit = 'gauss'
    vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label = analyze_histogram(test_name, mode, vars, pad, fit=fit)

    fig, ax = plt.subplots(figsize = (6,6))
    ax.set_ylabel(param_label)
    ax.set_xlabel('Density')
    ax.errorbar(rho_gases, vars, xerr=rho_gases_std, color='C0', ls='', marker='.')
    ax.errorbar(rho_liquids, vars, xerr=rho_liquids_std, color='C1', ls='', marker='.')
    ax.set_title("Histogram phase diagram")
    plt.savefig(f'{test_name}_histo_phase_diagram_frozen.png', dpi=300, bbox_inches='tight')
    plt.close()