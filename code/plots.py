import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *

matplotlib.rcParams.update({'font.size': 16})

def v(rho):
    return 5*np.exp(-1*np.tanh((rho-50)/20))

def F(r):
    rf = 1
    return np.where(100/rf*(1-r/rf)>0,100/rf*(1-r/rf),0)

rhos = np.linspace(0, 150, 1000)
rs = np.linspace(0,2.0,1000)
fig, ax = plt.subplots()
# ax.plot(rhos, v(rhos))
ax.plot(rs, F(rs))
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$F(r)$")
ax.set_title("Repulsive force")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
# ax.text(50, 10, r"$v(\rho)=v_0 exp\left(-\lambda \tanh \left(\frac{\rho-\rho_m}{\phi}\right)\right)$")
ax.text(0.8, 60, r"$F(r)=\frac{\epsilon}{r_f}\left(1-\frac{r}{r_f}\right)$")
plt.savefig("QSPF_visualize.png", dpi=300, bbox_inches='tight')


