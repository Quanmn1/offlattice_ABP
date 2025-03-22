import numpy as np
from scipy import optimize, stats, signal
from scipy.optimize import OptimizeWarning
from matplotlib import pyplot as plt
import os
import sys
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

def analyze_histogram(var, coarsen_number, fit='gauss'):
    test_name = f"pfqs_condensation_test10.1_{var}"

    param_file = f"{test_name}_param"
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
    box_size_pfap = params['r_max_pfap']
    density_box_raw_size = density_box_size * box_size
    density_box_area = density_box_raw_size**2
    number_of_boxes_x = Lx // density_box_raw_size // coarsen_number
    number_of_boxes_y = Ly // density_box_raw_size // coarsen_number
    threshold_density = 25

    # Coarsened histograms:
    new_density_box_area = density_box_area * coarsen_number**2
    max_number = int(new_density_box_area * 3 / box_size_pfap**2)
    densities = np.arange(0, max_number+1, 1) / new_density_box_area
    histogram = np.zeros(max_number+1)
    heatmap = np.zeros((number_of_boxes_y, number_of_boxes_x))

    for i in range(1,4):
        data_file = f"pfqs_condensation_test10.{i}_{var}_density"
        output_file = f"pfqs_condensation_test10.{i}_{var}_histogram_coarsened"
        with open(output_file, 'w') as f_output:
            # if this function is used to analyze sigmas profile: the data read in need to be coarsened, 
            # since the average is over cell of size 1
            count = 0 # index of the current segment
            row = 0
            t=0

            with open(data_file, 'r') as f:
                for raw_line in f:
                    line = [float(a) for a in raw_line.split()]
                    if len(line) == 1:
                        t = line[0]
                        row = 0
                        heatmap.fill(0)
                        histogram = np.zeros(max_number+1)
                        
                    elif len(line) > 1:
                        try:
                            heatmap[int(row), :] += np.array(coarsen(line, coarsen_number))/coarsen_number
                        except ValueError:
                            print(f"invalid density row {row} at var={var}, t={t}")
                            print(np.array(coarsen(line, coarsen_number))/coarsen_number)
                            raise
                        row += 1/coarsen_number

                    else:
                        if int(row) != number_of_boxes_y: # f"#rows = {row}, #boxes = {number_of_boxes_x}"
                            continue
                        
                        # print coarsened histogram to file
                        f_output.write(f"{t}\n")
                        for j in range(number_of_boxes_y):
                            for k in range(number_of_boxes_x):
                                index = np.searchsorted(densities, heatmap[j,k])
                                histogram[index] += 1

                        for j in range(max_number):
                            f_output.write(f"{densities[j]} {histogram[j]}\n")
                        f_output.write("\n")


    # rho_liquids, rho_gases, vars, rho_gases_std, rho_liquids_std = nonzero(rho_liquids, rho_gases, vars, rho_gases_std, rho_liquids_std)
    # Pes = v/Dr/vars
    # # vars = Pes

    # # print(rho_gases)
    # # print(rho_liquids)

    # with open('pfqs_condensation_test10_histo_phase_diagram_frozen', 'w') as f:
    #     f.write(f"{param_label} \t rho_gas \t rho_liquid \t rho_gas_std \t rho_liquid_std \n")
    #     for i in range(len(vars)):
    #         f.write(f"{vars[i]:.3f}\t{rho_gases[i]:.3f}\t{rho_liquids[i]:.3f}\t{rho_gases_std[i]:.3f}\t{rho_liquids_std[i]:.3f}\n")

    # return vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label


if __name__ == "__main__":
    # vars = np.array(sys.argv[1].split(),dtype=float)
    # fit = sys.argv[2]
    analyze_histogram(0.15, 2)

    # fig, ax = plt.subplots(figsize = (6,6))
    # ax.set_ylabel(param_label)
    # ax.set_xlabel('Density')
    # ax.errorbar(rho_gases, vars, xerr=rho_gases_std, color='C0', ls='', marker='.')
    # ax.errorbar(rho_liquids, vars, xerr=rho_liquids_std, color='C1', ls='', marker='.')
    # ax.set_title("Histogram phase diagram")
    # plt.savefig('pfqs_condensation_test10_histo_phase_diagram_frozen.png', dpi=300, bbox_inches='tight')
    # plt.close()