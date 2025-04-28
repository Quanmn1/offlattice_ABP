import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 14})

def v_QS(rho, v0, lamb, rhom, phi):
    return v0*np.exp(-lamb*np.tanh((rho-rhom)/phi))

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
    if mode == "pfap":
        file = test_name + f'_{rhos[0]:.0f}'
        params, num_particles = read_param(file + '_param')
    else:
        file = test_name + f'_{rhos[0]:.0f}'
        params, num_particles = read_param(file + '_param')

    # expectation: v(rho)*(PFAPS reduction)
    v_qs = v_QS(rhos, params['v'], params['lambda'], params['rho_m'], params['phi'])

    # Average over EVERY SNAPSHOTS that equilibriated
    vs = np.zeros((number, num_segments))
    v_stds = np.zeros((number, num_segments))
    vs_avg = np.zeros_like(rhos)
    vs_stds = np.zeros_like(rhos)
    rho_stds = np.zeros_like(rhos)
    rho_measureds = np.zeros_like(rhos)
    final_time = params['final_time']
    time_step = params['store_time_interval']
    N = num_particles * final_time // time_step
    equi_index = 20

    for (ind,rho) in enumerate(rhos):
        if mode == "pfap":
            file = test_name + f'_{rho:.0f}'
            params, num_particles = read_param(file + '_param')
        else:
            file = test_name + f'_{rho:.0f}'
            params, num_particles = read_param(file + '_param')
        
        datafile = file + '_v'
        try:
            data=np.loadtxt(datafile, usecols=(0,1))[num_particles:,:]
        except OSError:
            print(rho)
            continue
        
        # CONVERGENCE?
        for i in range(num_segments):
            data_i = data[N*i//num_segments:N*(i+1)//num_segments, :]
            # if simulation has not finished, those segments will has 0
            vs[ind, i] = np.mean(data_i[:,1])
            # v_stds[ind, i] = np.std(data_i[:,1])/np.sqrt(N/num_segments) # will not be used
        
        vs_rho, = nonzero(vs[ind,:])
        vs_avg[ind] = np.mean(vs_rho[equi_index:], axis=-1)
        vs_stds[ind] = np.std(vs_rho[equi_index:], axis=-1)/np.sqrt(len(vs_rho)-equi_index-1)
    
    # plot v vs rho
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$v_{eff}(\rho)$')
    # params, N = read_param(file + '_param')
    # v = params['v']
    # r_pf = params['r_max_pfap']
    # Dr = params['Dr']

    # if mode == "pfap":
    #     v = params['v']
    #     r_pf = params['r_max_pfap']
    #     rho_max = 2/np.sqrt(3)/r_pf/r_pf
    #     rho_dense = np.linspace(0,1.2,100)
    #     # v_expecteds = v*(1-rhos/rho_max)
    #     # ax.plot(rhos, v_expecteds, label=r"$v_0 (1-\rho / \rho_m)$")
    #     popt, pcov = optimize.curve_fit(veff, rhos[:-2], vs_avg[:-2], p0=(-v, rho_max), sigma=vs_stds[:-2], absolute_sigma=True)
    #     ax.plot(rho_dense,veff(rho_dense, *popt), label="Linear fit")
    #     print(pcov)
    #     ax.text(0.6, 0.5, fr"$v = {popt[0]:.2f}\times (1 - \rho/{popt[1]:.2f})$")
    # elif mode == "qsap":
    #     v = params['v']
    #     phi = params['phi']
    #     rho_m = params['rho_m']
    #     lamb = params['lambda']
    #     rho_dense = np.linspace(0, 60, 1000)
    #     v_expecteds = v*np.exp(-lamb*np.tanh((rho_dense-rho_m)/phi))
    #     ax.plot(rho_dense, v_expecteds, label='Expected veff')           
    # elif mode == "pfqs":
    #     v = params['v']
    #     phi = params['phi']
    #     rho_m = params['rho_m']
    #     lamb = params['lambda']
    #     r_pf = params['r_max_pfap']
    #     rho_dense = np.linspace(0, 60, 1000)
    #     def veff_fit(rho, rho_star):
    #         return veff_pfqs(rho, rho_star, r_pf, lamb, v, rho_m, phi)
    #     popt, pcov = optimize.curve_fit(veff_fit, rhos, vs, p0=(1.3), sigma=v_stds, absolute_sigma=True)     
    #     # ax.plot(rho_dense,veff_fit(rho_dense, *popt), label="Fit " + fr'$v(\rho)(1-\rho r_f^2/{popt[0]:.2f})$')   
    #     # v_expecteds = v*np.exp(-lamb*np.tanh((rho_dense-rho_m)/phi))*(1-rho_dense/(1.29/r_pf**2))
    #     # ax.plot(rho_dense, v_expecteds, label=r'$v(\rho)(1-\rho r_f^2/1.29)$')     

    # Check for convergence
    # ax.plot(rhos, vs, ls='', marker='.')
    # filename = test_name + '_v_convergence.png'

    # After converged, plot the reduced v
    ax.plot(rhos, vs_avg/v_qs, ls='', marker='.', label=r"$v^*(\rho)/v(\rho)$")
    filename = test_name + '_v_reduced.png'

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(r"$v_{eff}=\langle \dot{\vec{r}}\cdot \vec{u}\rangle$")
    final_time = params["final_time"]
    start_time = params["next_store_time"]
    # legends = [f"Time {(final_time-start_time)*i//num_segments+start_time}-{(final_time-start_time)*(i+1)//num_segments+start_time}" for i in range(num_segments)]
    ax.legend()
    plt.savefig(filename,  dpi=300, bbox_inches='tight')

    output = test_name + '_veff_reduced_data'
    with open(output, 'w') as f:
        f.write(f"rho \t v\n")
        for i in range(len(rhos)):
            f.write(f"{rhos[i]:.2f} \t {vs_avg[i]/v_qs[i]:.4f} \t {vs_stds[i]/v_qs[i]:.4f}\n")

if __name__ == "__main__":
    test_name = sys.argv[1]
    mode = sys.argv[2]
    vars = np.array(sys.argv[3].split(),dtype=float)
    num_segments = int(sys.argv[4])

    analyze_veff(test_name, mode, vars, num_segments)