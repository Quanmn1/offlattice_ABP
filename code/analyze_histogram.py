import numpy as np
from scipy import optimize, stats
from matplotlib import pyplot as plt
import os
import sys
from math_helper import *

def profile(x, rho_gas, rho_liquid, liquid_ratio, width_gas, width_liquid):
    return (1-liquid_ratio) * stats.norm.pdf(x, rho_gas, width_gas) + liquid_ratio * stats.norm.pdf(x, rho_liquid, width_liquid)

def fit_gausses(densities_coarse, histogram_mult, histogram_mult_std, params):
    rho_gas0 = params['rho'] * 0.6
    rho_liquid0 = params['rho'] * 1.6
    liquid_ratio0 = 0.5
    width_gas0 = params['rho'] * 0.1
    width_liquid0 = params['rho'] * 0.2

    histogram_std_fit, densities_fit, histogram_fit = nonzero(histogram_mult_std, densities_coarse, histogram_mult)

    popt, pcov, infodict, mesg, ier = optimize.curve_fit(profile, densities_fit, histogram_fit,
                                                         p0=(rho_gas0,rho_liquid0,liquid_ratio0,width_gas0,width_liquid0),
                                                         full_output=True)
    chisquare = np.sum(np.square(infodict['fvec']))
    dof = len(densities_coarse) - 5
    return popt, pcov, chisquare, dof

def process_histograms(histograms, num_combined, target_length):
    """
    Add multiple histograms together and coarse-grain the bins. Normalize the result.
    """
    histogram_mults = np.array(histograms)
    histogram_mult = np.average(histogram_mults, axis=0)
    histogram_mult_std = np.std(histogram_mults, axis=0)
    return histogram_mult, histogram_mult_std

start = float(sys.argv[2])
end = float(sys.argv[3])
space = float(sys.argv[4])
number = int((end-start)/space)+1
lambs = np.linspace(start, end, number)
test_name = sys.argv[1]
rho_gases = np.zeros(number)
rho_liquids = np.zeros(number)
rho_gases_std = np.zeros(number)
rho_liquids_std = np.zeros(number)

for (test_num, lamb) in enumerate(lambs):
    print(lamb)
    name = test_name + f'_{lamb:.1f}'
    param_file = name + '_param'
    data_file = name + '_data'
    histogram_file = name + '_histogram'
    density_file = name + '_density'

    # Read parameters into dict params, with string values
    with open(param_file, 'r') as f:
        params_info = f.readline().split()
        params_values = f.readline().split()
        params = {var:convert_num(value) for (var, value) in zip(params_info[2:], params_values[1:])}
        N = convert_num(f.readline().split()[0])

    N = params['N']
    # v = params['v']
    Lx = params['Lx']
    Ly = params['Ly']
    density_box_size = params['density_box_size']
    box_size = params['r_max']
    tf = params['final_time']
    density_box_raw_size = density_box_size * box_size
    density_box_area = density_box_raw_size**2
    number_of_boxes_x = Lx // density_box_raw_size
    number_of_boxes_y = Ly // density_box_raw_size
    params['rho'] = N/Lx/Ly

    max_number = int(N/(number_of_boxes_x*number_of_boxes_y) * 4)

    num_segments = 5
    segments = np.linspace(0,tf,num_segments+1)

    histogram_folder = name+'_histogram_folder'
    if not os.path.isdir(histogram_folder):
        os.mkdir(histogram_folder)

    densities = np.arange(0, max_number+1, 1) / density_box_area
    histogram = np.zeros(max_number+1)
    num_combined = 50
    binwidth = num_combined / density_box_area
    densities_coarse = coarsen(densities, num_combined)
    count = 0 # index of the current segment
    histograms = []

    with open(histogram_file, 'r') as f:
        for raw_line in f:
            line = [float(a) for a in raw_line.split()]
            if len(line) == 1:
                t = line[0]
                row = 0
                if t >= segments[count+1]:
                    count += 1

                    # fit and plot the profile
                    histogram_mult, histogram_mult_std = process_histograms(histograms, num_combined, len(densities_coarse))
                    popt, pcov, chisquare, dof = fit_gausses(densities_coarse, histogram_mult, histogram_mult_std, params)                 
                    rho_gases[test_num] = popt[0]
                    rho_liquids[test_num] = popt[1]
                    fig, ax = plt.subplots(figsize = (6,6))
                    ax.set_ylabel('Probability')
                    ax.set_xlabel('Density')
                    ax.bar(densities_coarse, histogram_mult, color='C0')
                    # ax.errorbar(densities_coarse, histogram_mult, yerr=histogram_mult_std, ls='')
                    ax.plot(densities, profile(densities, *popt), c='C1')
                    if np.all(np.isfinite(pcov)):
                        perr = np.sqrt(np.diag(pcov))
                        rho_gases_std[test_num] = perr[0]
                        rho_liquids_std[test_num] = perr[1]
                        ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}\pm {perr[0]:.3f}, \\rho_l = {popt[1]:.3f}\pm {perr[1]:.3f}$')
                    else:
                        ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}, \\rho_l = {popt[1]:.3f}$')
                        ax.text(1.1,0.6, 'Invalid fit', transform=ax.transAxes)
                    ax.text(0.85,0.95, f'$\chi^2 / dof = {chisquare:.2f}/{dof}$', transform=ax.transAxes)
                    plt.savefig(histogram_folder + '/' + name + f'_{count}.png', dpi=100, bbox_inches='tight')
                    plt.close()

                    histograms = []

            elif len(line) > 1:
                histogram[row] = line[1]
                row += 1
            else:
                histogram_coarsened = coarsen(histogram, num_combined, normalize=True)
                histograms.append(histogram_coarsened)
                histogram.fill(0)

rho_gases, rho_liquids, lambs, rho_gases_std, rho_liquids_std = nonzero(rho_gases, rho_liquids, lambs, rho_gases_std, rho_liquids_std)
# ind = np.nonzero(rho_gases)
# rho_gases = rho_gases[ind]
# rho_liquids = rho_liquids[ind]
# lambs = lambs[ind]
# rho_gases_std = rho_gases_std[ind]
# rho_liquids_std = rho_liquids_std[ind]

with open(test_name + '_phase_diagram', 'w') as f:
    f.write("Lambda \t rho_gas \t rho_liquid\n")
    for i in range(len(lambs)):
        f.write(f"{lambs[i]:.2f} \t {rho_gases[i]:.3f} \t {rho_liquids[i]:.3f}\n")

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel('$\lambda$')
ax.set_xlabel('Density')
ax.errorbar(rho_gases, lambs, xerr=rho_gases_std, color='C0', ls='', marker='.')
ax.errorbar(rho_liquids, lambs, xerr=rho_liquids_std, color='C1', ls='', marker='.')
ax.set_title("Phase diagram")
plt.savefig(test_name + '_histo_phase_diagram.png', dpi=300, bbox_inches='tight')
plt.close()