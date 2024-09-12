import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 14})


def veff(rho, v, phi):
    return v*(1-rho/phi)

def Pa_fit(rho, s1, s2, s3, s4):
    lp = 10
    v = 5
    return v * (1-s1*rho+s2*rho**2) * (1-np.tanh(s3*(rho-s4)))/2

def analyze_veff(test_name, mode, vars, num_segments):
    # input: sequence of rhos.
    # output: plot v(rho)
    number = len(vars)
    rhos = vars
    # for each rho, get last state data, average the v's

    # Average over EVERY SNAPSHOTS

    vs = np.zeros((number, num_segments))
    v_stds = np.zeros((number, num_segments))
    # rho_stds = np.zeros_like(rhos)
    # rho_measureds = np.zeros_like(rhos)

    for (ind,rho) in enumerate(rhos):
        if mode == "pfap":
            file = test_name + f'_{rho:.1f}'
        else:
            file = test_name + f'_{rho:.0f}'
        
        datafile = file + '_v'
        try:
            data=np.loadtxt(datafile, usecols=(0,1))
        except OSError:
            print(rho)
            continue
        
        # CONVERGENCE?
        N = data.shape[0]
        for i in range(num_segments):
            data_i = data[N*i//num_segments:N*(i+1)//num_segments, :]
            vs[ind, i] = np.mean(data_i[:,1])
            v_stds[ind, i] = np.std(data_i[:,1])/np.sqrt(N/num_segments)

        # rho_is = data[:,0]
        # v_is = data[:,1]
        # N = len(v_is)
        print(N)        
        # rho_measureds[ind] = np.average(rho_is)
        # rho_stds[ind] = np.std(rho_is)/np.sqrt(N)
        # obtain v average
        # vs[ind] = np.average(v_is)
        # v_stds[ind] = np.std(v_is)/np.sqrt(N)

    # check consistency
    # print(rhos)
    # print(rho_measureds)
    # print(rho_stds)
    # print((rhos-rho_measureds)/rho_stds)

    # plot v vs rho
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$v_{eff}(\rho)$')
    params, N = read_param(file + '_param')
    v = params['v']
    r_pf = params['r_max_pfap']
    Dr = params['Dr']
    vs_avg = np.mean(vs, axis=-1)
    vs_stds = np.std(vs, axis=-1)/np.sqrt(num_segments-1)

    if mode == "pfap":
        v = params['v']
        r_pf = params['r_max_pfap']
        rho_max = 2/np.sqrt(3)/r_pf/r_pf
        rho_dense = np.linspace(0,1.2,100)
        # v_expecteds = v*(1-rhos/rho_max)
        # ax.plot(rhos, v_expecteds, label=r"$v_0 (1-\rho / \rho_m)$")
        popt, pcov = optimize.curve_fit(veff, rhos[:-2], vs_avg[:-2], p0=(-v, rho_max), sigma=vs_stds[:-2], absolute_sigma=True)
        ax.plot(rho_dense,veff(rho_dense, *popt), label="Linear fit")
        print(pcov)
        ax.text(0.6, 0.5, fr"$v = {popt[0]:.2f}\times (1 - \rho/{popt[1]:.2f})$")
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
        rho_dense = np.linspace(0, 60, 1000)
        def veff_fit(rho, rho_star):
            return veff_pfqs(rho, rho_star, r_pf, lamb, v, rho_m, phi)
        popt, pcov = optimize.curve_fit(veff_fit, rhos, vs, p0=(1.3), sigma=v_stds, absolute_sigma=True)     
        # ax.plot(rho_dense,veff_fit(rho_dense, *popt), label="Fit " + fr'$v(\rho)(1-\rho r_f^2/{popt[0]:.2f})$')   
        # v_expecteds = v*np.exp(-lamb*np.tanh((rho_dense-rho_m)/phi))*(1-rho_dense/(1.29/r_pf**2))
        # ax.plot(rho_dense, v_expecteds, label=r'$v(\rho)(1-\rho r_f^2/1.29)$')     

    # ax.errorbar(rhos, vs, yerr=v_stds, xerr=rho_stds, ls='', marker='.')
    ax.plot(rhos, vs, ls='', marker='.', label=r"$v^*(\rho)/v(\rho)$")
    ax.legend()
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(r"$v_{eff}=\langle \dot{\vec{r}}\cdot \vec{u}\rangle$")
    final_time = params["final_time"]
    start_time = params["next_store_time"]
    legends = [f"Time {(final_time-start_time)*i//num_segments+start_time}-{(final_time-start_time)*(i+1)//num_segments+start_time}" for i in range(num_segments)]
    ax.legend(legends)
    plt.savefig(test_name + '_v.png',  dpi=300, bbox_inches='tight')

    output = test_name + '_veff_data'
    with open(output, 'w') as f:
        f.write(f"rho \t v\n")
        for i in range(len(rhos)):
            f.write(f"{rhos[i]:.2f} \t {vs_avg[i]:.4f} \t {vs_stds[i]:.4f}\n")

if __name__ == "__main__":
    test_name = sys.argv[1]
    mode = sys.argv[2]
    vars = np.array(sys.argv[3].split(),dtype=float)
    num_segments = int(sys.argv[4])

    analyze_veff(test_name, mode, vars, num_segments)