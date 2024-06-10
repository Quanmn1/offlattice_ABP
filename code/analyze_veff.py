import numpy as np
from scipy import optimize, stats, signal
from matplotlib import pyplot as plt
import os
import sys
from helper import *

def analyze_veff(test_name, mode, start, end, space):
    # input: sequence of rhos.
    # output: plot v(rho)
    number = round((end-start)/space)+1
    rhos = np.linspace(start, end, number)
    # for each rho, get last state data, average the v's
    vs = np.zeros_like(rhos)
    v_stds = np.zeros_like(rhos)
    rho_stds = np.zeros_like(rhos)
    rho_measureds = np.zeros_like(rhos)

    for (ind,rho) in enumerate(rhos):
        if mode == "pfap":
            file = test_name + f'_{rho:.1f}'
            datafile = file + '_video/' + test_name + f'_{rho:.1f}_last_state'
        else:
            file = test_name + f'_{rho:.0f}'
            datafile = file + '_video/' + test_name + f'_{rho:.0f}_last_state'        
        
        try:
            data=np.loadtxt(datafile, skiprows=1, usecols=(3,4))
        except OSError:
            continue
        params, N = read_param(file + '_param')
        rho_is = data[:,0]
        v_is = data[:,1]
        N = len(v_is)
        # print(np.max(v_is))
        # print(np.min(v_is))
        
        # plot the spread of rho_i and v_i
        # plt.hist(v_is, bins=100)
        plt.figure(figsize=(6,12))
        indices = (v_is < 1) & (v_is > -1.5)
        plt.scatter(rho_is[indices], v_is[indices])
        plt.xlabel(r"$\rho$")
        plt.ylabel("v")
        plt.title(fr"$\rho={rho}$")
        plt.savefig(file + 'vspread.png', dpi=100, bbox_inches='tight')
        plt.close()
        rho_measureds[ind] = np.average(rho_is)
        rho_stds[ind] = np.std(rho_is)/np.sqrt(N)
        # obtain v average
        vs[ind] = np.average(v_is)
        v_stds[ind] = np.std(v_is)/np.sqrt(N)

    # check consistency
    # print(rhos)
    # print(rho_measureds)
    # print(rho_stds)
    # print((rhos-rho_measureds)/rho_stds)

    # plot v vs rho
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel('v')

    # put functional form in label, or maybe in title
    # put parameters
    # put as much info in figures as possible
    # do latex
    # make the expected smoother
    # make sure to switch off QS when measure PF
        
    if mode == "pfap":
        v = params['v']*np.exp(-0.4*np.tanh((0.78-25)/10))
        r_pf = params['r_max_pfap']
        rho_max = 2/np.sqrt(3)/r_pf/r_pf
        rho_dense = np.linspace(0,1.1,100)
        # v_expecteds = v*(1-rhos/rho_max)
        # ax.plot(rhos, v_expecteds, label=r"$v_0 (1-\rho / \rho_m)$")
        popt, pcov = optimize.curve_fit(linear, rhos, vs, p0=(-v/rho_max, v), sigma=v_stds, absolute_sigma=True)
        ax.plot(rho_dense, linear(rho_dense, *popt), label="Linear fit")
        perr = np.sqrt(np.diag(pcov))
        ax.text(0.7, 0.5, fr"$v = {popt[1]:.2f}\times (1 - \rho/{-popt[1]/popt[0]:.2f})$")
    elif mode == "qsap":
        v = params['v']
        phi = params['phi']
        rho_m = params['rho_m']
        lamb = params['lambda']
        rho_dense = np.linspace(0, 60, 1000)
        v_expecteds = v*np.exp(-lamb*np.tanh((rho_dense-rho_m)/phi))
        ax.plot(rho_dense, v_expecteds, label='Expected veff')           
    elif mode == "pfqs":
        v = params['v']
        phi = params['phi']
        rho_m = params['rho_m']
        lamb = params['lambda']
        r_pf = params['r_max_pfap']
        v_expecteds = v*np.exp(-lamb*np.tanh((rho_measureds-rho_m)/phi)) - v + v*(1-rhos/(1.29/r_pf**2))
        ax.plot(rhos, v_expecteds, label='Expected veff')     

    ax.errorbar(rhos, vs, yerr=v_stds, xerr=rho_stds, ls='', marker='.', label="Measured veff")
    ax.legend()
    ax.set_xlim(left=0)
    plt.savefig(test_name + '_result.png',  dpi=100, bbox_inches='tight')

    output = test_name + '_veff_data'
    np.savetxt(output, (rhos-rho_measureds)/rho_stds)

    with open(output, 'a') as f:
        f.write(f"rho \t v\n")
        for i in range(len(rhos)):
            f.write(f"{rho_measureds[i]:.2f} \t {vs[i]:.4f}\n")


if __name__ == "__main__":
    test_name = sys.argv[1]
    mode = sys.argv[2]
    start = float(sys.argv[3])
    end = float(sys.argv[4])
    space = float(sys.argv[5])

    analyze_veff(test_name, mode, start, end, space)