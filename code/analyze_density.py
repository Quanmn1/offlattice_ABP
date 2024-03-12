import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import os
import sys
from math_helper import *

def profile(x, rho_gas, rho_liquid, x1, x2, slope1, slope2):
    # use a PBC-friendly function?
    return np.tanh(slope1*(x-x1)) * np.tanh(slope2*(x2-x)) * (rho_liquid-rho_gas)/2 + (rho_liquid+rho_gas)/2

def fit_tanh(x_lattice, one_dim_profile, one_dim_profile_std, params):
    rho_gas0 = params['rho_small']
    rho_liquid0 = params['rho_large']
    x10 = params['Lx']*(0.5-params['liquid_fraction']/2)
    x20 = params['Lx']*(0.5+params['liquid_fraction']/2)
    # slopes: might be very large/undefined bc the density resolution is low
    # and the slope is of order interaction range
    slope1 = 1/params['r_max']/params['density_box_size']
    slope2 = slope1
    sigmas = one_dim_profile_std + 1/(params['r_max']*params['density_box_size'])**2

    popt, pcov, infodict, mesg, ier = optimize.curve_fit(profile, x_lattice, one_dim_profile,
                                                         p0=(rho_gas0,rho_liquid0,x10,x20,slope1,slope2),
                                                         sigma=sigmas,absolute_sigma=True,full_output=True)
    chisquare = np.sum(np.square(infodict['fvec']))
    dof = len(x_lattice) - 6
    return popt, pcov, chisquare, dof

def combine_profiles(x_lattice, one_dim_profiles, one_dim_profile_stds):
    """
    Combine multiple profiles with same x_lattices into a fit.
    one_dim_profiles is a list of 1D arrays.
    Calculate the CoM of each profile, translate so that they overlap, 
    and combine into a single fit.
    Good fit if x10+x20 == Lx, slope1 == slope2
    """
    # 
    Lx = params['Lx']
    x_lattice_multiple = []
    for profile in one_dim_profiles:
        center_of_mass = center_of_mass_pbc(x_lattice, Lx, profile)
        translated = translate_pbc(x_lattice, Lx, Lx/2-center_of_mass)
        x_lattice_multiple.append(translated)
    x_lattice_mult = np.concatenate(x_lattice_multiple)
    one_dim_profile_mult = np.concatenate(one_dim_profiles)
    one_dim_profile_std_mult = np.concatenate(one_dim_profile_stds)
    return x_lattice_mult, one_dim_profile_mult, one_dim_profile_std_mult


start = float(sys.argv[2])
end = float(sys.argv[3])
space = float(sys.argv[4])
number = int((end-start)/space)+1
Drs = np.linspace(start, end, number)
test_name = sys.argv[1]
rho_gases = np.zeros(number)
rho_liquids = np.zeros(number)
rho_gases_std = np.zeros(number)
rho_liquids_std = np.zeros(number)

for (test_num, Dr) in enumerate(Drs):
    name = test_name + f'_{Dr:.3f}'
    param_file = name + '_param'
    data_file = name + '_data'
    histogram_file = name + '_histogram'
    density_file = name + '_density'

    # Read parameters into dict params, with string values
    with open(param_file, 'r') as f:
        params_info = f.readline().split()
        params_values = f.readline().split()
        params = {var:convert_num(value) for (var, value) in zip(params_info[2:], params_values[1:])}

    # N = params['N']
    # v = params['v']
    Lx = params['Lx']
    Ly = params['Ly']
    v = params['v']
    density_box_size = params['density_box_size']
    box_size = params['r_max']
    tf = params['final_time']
    density_box_raw_size = density_box_size * box_size
    number_of_boxes_x = Lx // density_box_raw_size
    number_of_boxes_y = Ly // density_box_raw_size

    num_segments = 5
    segments = np.linspace(0,tf,num_segments+1)

    density_folder = name+'_avg_density_folder'
    if not os.path.exists(density_folder):
        os.mkdir(density_folder)

    heatmap = np.zeros((number_of_boxes_y, number_of_boxes_x))
    x_lattice = np.linspace(density_box_raw_size/2, Lx-density_box_raw_size/2, number_of_boxes_x)
    fine_lattice = np.linspace(0, Lx, 100)
    count = 0 # index of the current segment
    row = 0
    one_dim_profiles = []
    one_dim_profile_stds = []

    with open(density_file, 'r') as f:
        for raw_line in f:
            line = [float(a) for a in raw_line.split()]
            if len(line) == 1:
                t = line[0]
                row = 0
                if t >= segments[count+1]:
                    count += 1

                    # fit and plot the profile
                    x_lattice_mult, one_dim_profile_mult, one_dim_profile_std_mult = combine_profiles(x_lattice, one_dim_profiles, one_dim_profile_stds)
                    popt, pcov, chisquare, dof = fit_tanh(x_lattice_mult, one_dim_profile_mult, one_dim_profile_std_mult, params)                 
                    rho_gases[test_num] = popt[0]
                    rho_liquids[test_num] = popt[1]
                    fig, ax = plt.subplots(figsize = (6,6))
                    ax.set_xlabel('x')
                    ax.set_ylabel('Density')
                    ax.set_xlim(0, Lx)
                    ax.scatter(x_lattice_mult, one_dim_profile_mult, s=0.1, c='C0')
                    # ax.errorbar(x_lattice_mult, one_dim_profile_mult, yerr=one_dim_profile_std_mult, ls='')
                    ax.plot(fine_lattice, profile(fine_lattice, *popt), c='C1')
                    if np.all(np.isfinite(pcov)):
                        perr = np.sqrt(np.diag(pcov))
                        rho_gases_std[test_num] = perr[0]
                        rho_liquids_std[test_num] = perr[1]
                        ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}\pm {perr[0]:.3f}, \\rho_l = {popt[1]:.3f}\pm {perr[1]:.3f}$')
                    else:
                        ax.set_title(f'$t={t:.2f}, \\rho_g = {popt[0]:.3f}, \\rho_l = {popt[1]:.3f}$')
                        ax.text(1.1,0.6, 'Invalid fit', transform=ax.transAxes)
                    ax.text(0.85,0.95, f'$\chi^2 / dof = {chisquare:.2f}/{dof}$', transform=ax.transAxes)
                    plt.savefig(density_folder + '/' + name + f'_{count}.png', dpi=100, bbox_inches='tight')
                    plt.close()

                    one_dim_profiles = []
                    one_dim_profile_stds = []

            elif len(line) > 1:
                heatmap[row, :] = line
                row += 1
            else:
                if row != number_of_boxes_y: # f"#rows = {row}, #boxes = {number_of_boxes_x}"
                    continue
                one_dim_profile = np.average(heatmap, axis=0)
                one_dim_profile_std = np.std(heatmap, axis=0) # Or std of the mean?
                one_dim_profiles.append(one_dim_profile)
                one_dim_profile_stds.append(one_dim_profile_std)
                heatmap.fill(0)

ind = np.nonzero(rho_gases)
rho_gases = rho_gases[ind]
rho_liquids = rho_liquids[ind]
Drs = Drs[ind]
rho_gases_std = rho_gases_std[ind]
rho_liquids_std = rho_liquids_std[ind]
Pes = v/Drs*np.float_power(2, 1/6)

with open(test_name + '_phase_diagram', 'w') as f:
    f.write("Pe \t rho_gas \t rho_liquid\n")
    for i in range(len(Pes)):
        f.write(f"{Pes[i]:.2f} \t {rho_gases[i]:.3f} \t {rho_liquids[i]:.3f}\n")

fig, ax = plt.subplots(figsize = (6,6))
ax.set_ylabel('Pe')
ax.set_xlabel('Density')
ax.errorbar(rho_gases, Pes, xerr=rho_gases_std, color='C0', ls='', marker='.')
ax.errorbar(rho_liquids, Pes, xerr=rho_liquids_std, color='C1', ls='', marker='.')
ax.set_title("Phase diagram")
plt.savefig(test_name + '_phase_diagram.png', dpi=300, bbox_inches='tight')
plt.close()