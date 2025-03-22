import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import os

matplotlib.rcParams.update({'font.size': 20})
lw = 3
rhos = np.linspace(0, 50, 1000)
vs = 5 * np.exp(-1*np.tanh((rhos-25)/1)) - 5
# plt.subplots(figsize = (6,6))
plt.xlabel(r"$\rho$")
plt.ylabel(r"$v$")
# plt.text(23, 10, r"$v(\rho) = v_0 \exp \left[-\lambda \tanh \left(\frac{\rho-\rho_m}{\phi} \right) \right]$")
plt.title(r"$v_2(\rho)$")
plt.plot(rhos, vs, c='C2', lw=lw)
plt.xlim(left=0, right=50)
plt.ylim(bottom=0, top=15)
plt.savefig("Quorum sensing funcion zero.png", dpi=300, bbox_inches='tight')
plt.close()

rs = np.linspace(0.0, 2.0, 1000)
# fs = np.where(rs <= 1, -12*0.125*(1/rs**7 - 1/rs**13), 0)
Vs = np.where(rs <= 1, 100*(1-rs)**2, 0)
# plt.text(1.2, 50.0, r"$V(r)=\epsilon \left(1-\frac{r}{r_f}\right)^2$")
plt.xlabel(r"$r$")
plt.ylabel(r"$V(r)$")
plt.title(r"Steric potential $V(r)$")
plt.plot(rs, Vs, c='C0', lw=lw)
plt.axvline(1.0, ls='--', c='C1', lw=lw)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig("Harmonic potential.png", dpi=300, bbox_inches='tight')
plt.close()

rs = np.linspace(0.0, 2.0, 1000)
# fs = np.where(rs <= 1, -12*0.125*(1/rs**7 - 1/rs**13), 0)
Vs = np.where(rs <= 1, np.exp(-1/(1-rs**2))/0.46651239317833015, 0)
# plt.text(1.2, 1.0, r"$\rho(r) = \int K(r-r') \rho_{micro}(r') = \sum_i \frac{1}{Z}\exp\left(-\frac{r_0^2}{r_0^2-r_i^2}\right)$")
plt.xlabel(r"$r$")
plt.ylabel(r"$K(r)$")
plt.title(r"Quorum sensing kernel $K(r)$")
plt.plot(rs, Vs, c='C0', lw=lw)
plt.axvline(1.0, ls='--', c='C1', lw=lw)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig("Kernel.png", dpi=300, bbox_inches='tight')
plt.close()

rhos = np.linspace(0, 50, 1000)
def vzero(rhos):
    v = np.zeros_like(rhos)
    for i, rho in enumerate(rhos):
        if rho < 24:
            v[i]=5
        elif rho < 25:
            v[i] = 5 - (rho-24)*5
        else:
            break
    return v

# plt.xlabel(r"$\rho$")
# plt.ylabel(r"$v$")
# plt.title(r"$v_1(\rho)$")
# plt.plot(rhos, vzero(rhos), c='C2', lw=lw)
# plt.axvline(25.0, ls='--', c="C3", lw=lw)
# plt.xlim(left=0)
# plt.ylim(bottom=0, top=6)
# plt.savefig("Quorum sensing funcion zero.png", dpi=300, bbox_inches='tight')
# plt.close()
