import numpy as np
from scipy import optimize, stats, signal
from matplotlib import pyplot as plt
import os
import sys
from helper import *

def profile_trimodal(x, rho_gas, rho_liquid, liquid_ratio, width_gas, width_liquid, rho_mid, width_mid, mid_ratio):
    return (1-mid_ratio)*(1-liquid_ratio) * stats.norm.pdf(x, rho_gas, width_gas) + (1-mid_ratio)*liquid_ratio * stats.norm.pdf(x, rho_liquid, width_liquid) \
            + mid_ratio * stats.norm.pdf(x, rho_mid, width_mid)

def profile(x, rho_gas, rho_liquid, liquid_ratio, width_gas, width_liquid):
    return (1-liquid_ratio) * stats.norm.pdf(x, rho_gas, width_gas) + liquid_ratio * stats.norm.pdf(x, rho_liquid, width_liquid)

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

def process_histograms(histograms, num_combined, target_length):
    """
    Take average and standard deviation of multiple normalized histograms.
    """
    histogram_mults = np.array(histograms)
    histogram_mult = np.average(histogram_mults, axis=0)
    histogram_mult_std = np.std(histogram_mults, axis=0)
    return histogram_mult, histogram_mult_std

def analyze_histogram(test_name, mode, start, end, space, num_segments, fit='gauss'):
    number = round((end-start)/space)+1
    if mode == "qsap":
        pad = 1
        param_label = "lambda"
        num_combined = 25 # number of densities binned together
    elif mode == "pfap":
        pad = 3
        param_label = r"$l_p/r_f$"
        num_combined = 1
    elif mode == "pfqs":
        pad = 2
        param_label = r"$l_p/r_f$"
        num_combined = 4
    vars = np.linspace(start, end, number)
    rho_gases = np.zeros(number)
    rho_liquids = np.zeros(number)
    rho_gases_std = np.zeros(number)
    rho_liquids_std = np.zeros(number)

    for (test_num, var) in enumerate(vars):
        name = test_name + '_' + f'{var:0.{pad}f}'
        param_file = name + '_param'
        histogram_file = name + '_histogram'

        # Read parameters into dict params, with string values
        try:
            params, N = read_param(param_file)
        except FileNotFoundError:
            continue

        # N = params['N']
        Lx = params['Lx']
        Ly = params['Ly']
        density_box_size = params['density_box_size']
        try:
            if mode == "pfap":
                box_size = params['r_max_pfap']
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

        if mode == "pfap":
            v = params['v']
            max_number = int(density_box_area * 2 / box_size**2)
        elif mode == "qsap":
            max_number = int(N/(number_of_boxes_x*number_of_boxes_y) * 4)
        elif mode == "pfqs":
            max_number = int(density_box_area * 2 / box_size_pfap**2)
        segments = np.linspace(0,tf,num_segments+1)

        histogram_folder = name+'_histogram_folder'
        if not os.path.isdir(histogram_folder):
            os.mkdir(histogram_folder)

        densities = np.arange(0, max_number+1, 1) / density_box_area
        histogram = np.zeros(max_number*2)
        binwidth = num_combined / density_box_area
        densities_coarse = coarsen(densities, num_combined)
        count = 0 # index of the current segment
        histograms = []
        rho_liquid0 = params['rho'] * 1.6
        rho_gas0 = params['rho'] * 0.6

        with open(histogram_file, 'r') as f:
            for raw_line in f:
                line = [float(a) for a in raw_line.split()]
                if len(line) == 1:
                    t = line[0]
                    row = 0
                    if t >= segments[count+1]:
                        count += 1
                        fig, ax = plt.subplots(figsize = (6,6))
                        ax.set_ylabel('Probability')
                        ax.set_xlabel('Density')

                        histogram_mult, histogram_mult_std = process_histograms(histograms, num_combined, len(densities_coarse))
                        histogram_std_fit, densities_fit, histogram_fit = nonzero(histogram_mult_std, densities_coarse, histogram_mult)
                        ax.errorbar(densities_fit, histogram_fit, yerr=histogram_std_fit, color='C0', label="Time-averaged histogram")
                        
                        if fit == 'gauss':
                            # fit and plot the profile
                            try:
                                popt, pcov, chisquare, dof = fit_gausses(densities_fit, histogram_fit, histogram_std_fit, params, rho_gas0, rho_liquid0)                 
                            except RuntimeError:
                                histograms = []
                                continue
                            rho_gases[test_num] = popt[0]
                            rho_liquids[test_num] = popt[1]
                            rho_gas0, rho_liquid0 = popt[0], popt[1] # initial guess for next time segment
                            densities_dense = np.linspace(np.min(densities_fit), np.max(densities_fit), 1000)
                            ax.plot(densities_dense, profile(densities_dense, *popt), c='C1', label="Bimodal Gaussian fit")
                            if np.all(np.isfinite(pcov)):
                                perr = np.sqrt(np.diag(pcov))
                                rho_gases_std[test_num] = perr[0]
                                rho_liquids_std[test_num] = perr[1]
                                ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}\pm {perr[0]:.3f}, \\rho_l = {popt[1]:.3f}\pm {perr[1]:.3f}$')
                            else:
                                ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}, \\rho_l = {popt[1]:.3f}$')
                                ax.text(1.1,0.6, 'Invalid fit', transform=ax.transAxes)
                            ax.text(0.85,0.95, f'$\chi^2 / dof = {chisquare:.2f}/{dof}$', transform=ax.transAxes)
                            # put on a grid
                        elif fit == 'max':
                            divider = np.searchsorted(densities_fit, params['rho']*1.5)
                            rho_gases[test_num] = densities_fit[np.argmax(histogram_fit[:divider])]
                            rho_liquids[test_num] = densities_fit[divider+np.argmax(histogram_fit[divider:])]
                            ax.set_title(f'$t={t:.2f}, peaks: \\rho_g = {rho_gases[test_num]:.3f}, \\rho_l = {rho_liquids[test_num]:.3f}$')
                        ax.legend()
                        ax.grid()
                        plt.savefig(histogram_folder + '/' + name + f'_{count}.png', dpi=100, bbox_inches='tight')
                        plt.close()

                        histograms = []

                elif len(line) > 1:
                    histogram[row] = line[1]
                    row += 1
                else:
                    histogram_coarsened = coarsen(histogram, num_combined, binwidth, normalize=True)
                    histograms.append(histogram_coarsened)
                    histogram.fill(0)

    rho_gases, rho_liquids, vars, rho_gases_std, rho_liquids_std = nonzero(rho_gases, rho_liquids, vars, rho_gases_std, rho_liquids_std)
    if mode == "pfap":
        Pes = v/vars*np.float_power(2, 1/6)
        vars = Pes
    elif mode == "qsap":
        pass
    elif mode == "pfqs":
        Pes = 5/vars
        vars = Pes


    with open(test_name + '_histo_phase_diagram', 'w') as f:
        f.write(f"{param_label} \t rho_gas \t rho_liquid\n")
        for i in range(len(vars)):
            f.write(f"{vars[i]:.2f} \t {rho_gases[i]:.3f} \t {rho_liquids[i]:.3f}\n")

    return vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label


if __name__ == "__main__":
    start = float(sys.argv[3])
    end = float(sys.argv[4])
    space = float(sys.argv[5])
    test_name = sys.argv[1]
    mode = sys.argv[2]
    num_segments = int(sys.argv[6])
    if len(sys.argv) > 7:
        fit = sys.argv[7]

    vars, rho_gases, rho_liquids, rho_gases_std, rho_liquids_std, param_label = analyze_histogram(test_name, mode, start, end, space, num_segments, fit=fit)

    fig, ax = plt.subplots(figsize = (6,6))
    ax.set_ylabel(param_label)
    ax.set_xlabel('Density')
    ax.errorbar(rho_gases, vars, xerr=rho_gases_std, color='C0', ls='', marker='.')
    ax.errorbar(rho_liquids, vars, xerr=rho_liquids_std, color='C1', ls='', marker='.')
    ax.set_title("Histogram phase diagram")
    plt.savefig(test_name + '_histo_phase_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()