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
Coarsen the density histogram into larger boxes
"""

def process_histograms(histograms):
    """
    Take average and standard deviation of multiple normalized histograms.
    """
    histogram_mults = np.array(histograms)
    histogram_mult = np.average(histogram_mults, axis=0)
    histogram_mult_std = np.std(histogram_mults, axis=0)
    return histogram_mult, histogram_mult_std

def analyze_histogram_coarsen(test_name, vars, coarsen_number, pad):
    param_file = f"{test_name}_{vars[0]:.{pad}f}_param"
    try:
        params, N = read_param(param_file)
    except FileNotFoundError:
        print("File not found")
        print(vars[0])

    Lx = params['Lx']
    Ly = params['Ly']
    density_box_size = params['density_box_size']
    box_size = params['r_max_qsap']
    density_box_raw_size = density_box_size * box_size
    density_box_area = density_box_raw_size**2
    number_of_boxes_x = Lx // density_box_raw_size // coarsen_number
    number_of_boxes_y = Ly // density_box_raw_size // coarsen_number

    # Coarsened histograms:
    new_density_box_area = density_box_area * coarsen_number**2
    heatmap = np.zeros((number_of_boxes_y, number_of_boxes_x))

    for var in vars:
        param_file = f"{test_name}_{var:.{pad}f}_param"
        try:
            params, N = read_param(param_file)
        except FileNotFoundError:
            print("File not found")
            print(var)
        box_size_pfap = params['r_max_pfap']
        max_number = int(new_density_box_area * 3 / box_size_pfap**2)
        densities = np.arange(0, max_number+1, 1) / new_density_box_area
        histogram = np.zeros(max_number+1)

        data_file = f"{test_name}_{var:.{pad}f}_density"
        output_file = f"{test_name}_{var:.{pad}f}_histogram_coarsened_{coarsen_number}"
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
                            heatmap[int(row), :] += np.array(coarsen(line, coarsen_number))
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
                                histogram[int(heatmap[j,k])] += 1

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
    test_name = sys.argv[1]
    vars = np.array(sys.argv[2].split(),dtype=float)
    coarsen_number = int(sys.argv[3])
    pad = int(sys.argv[4])
    # vars = np.array(sys.argv[1].split(),dtype=float)
    # fit = sys.argv[2]
    analyze_histogram_coarsen(test_name, vars, coarsen_number, pad)

    # fig, ax = plt.subplots(figsize = (6,6))
    # ax.set_ylabel(param_label)
    # ax.set_xlabel('Density')
    # ax.errorbar(rho_gases, vars, xerr=rho_gases_std, color='C0', ls='', marker='.')
    # ax.errorbar(rho_liquids, vars, xerr=rho_liquids_std, color='C1', ls='', marker='.')
    # ax.set_title("Histogram phase diagram")
    # plt.savefig('pfqs_condensation_test10_histo_phase_diagram_frozen.png', dpi=300, bbox_inches='tight')
    # plt.close()