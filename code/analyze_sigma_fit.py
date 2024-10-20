import numpy as np
from scipy import optimize, stats, signal
import matplotlib
from matplotlib import pyplot as plt
import os
import sys
from helper import *
from pysr import PySRRegressor
import sympy

"""
Use symbolic regression to fit an expression to pressure data.
"""

matplotlib.rcParams.update({'font.size': 16})

model_pC = PySRRegressor(
    niterations=200,  # < Increase me for better results
    population_size=100,
    # ^ Slightly larger populations, for greater diversity.
    maxsize=30,
    # ^ Allow greater complexity.
    binary_operators=["+", "*", "-", "/"],
    unary_operators=[
        "exp",
        "square",
        "cube",
        "invsqrt(x::T) where {T} = x>0 ? 1/sqrt(x) : T(NaN)"
        # Custom operator (julia syntax)
    ],
    # Define operator for SymPy as well
    extra_sympy_mappings={"invsqrt": lambda x: 1/sympy.sqrt(x)},
    # elementwise_loss="loss(prediction, target) = (prediction - target)^2",
    # ^ Custom loss function (julia syntax)
    nested_constraints={
        "square": {"square": 1, "cube": 0, "exp": 1, "invsqrt": 0},
        "cube": {"square": 1, "cube": 0, "exp": 1, "invsqrt": 0},
        "exp": {"square": 1, "cube": 1, "exp": 0, "invsqrt": 0},
        "invsqrt": {"square": 1, "cube": 1, "exp": 0, "invsqrt": 0},
    },
)

model_pA = PySRRegressor(
    niterations=100,  # < Increase me for better results
    population_size=60,
    # ^ Slightly larger populations, for greater diversity.
    maxsize=25,
    # ^ Allow greater complexity.
    binary_operators=["+", "*", "-", "/"],
    unary_operators=[
        "mexp(x)=exp(-x)",
        "square",
        "cube",
        "mtanh(x)=1-tanh(x)"
        # Custom operator (julia syntax)
    ],
    # Define operator for SymPy as well
    extra_sympy_mappings={"mexp": lambda x: sympy.exp(-x), "mtanh": lambda x: 1-sympy.tanh(x)},
    # elementwise_loss="loss(prediction, target) = (prediction - target)^2",
    # ^ Custom loss function (julia syntax)
    nested_constraints={
        "square": {"square": 1, "cube": 1, "mexp": 1, "mtanh": 2},
        "cube": {"square": 0, "cube": 0, "mexp": 1, "mtanh": 2},
        "mexp": {"square": 1, "cube": 1, "mexp": 0, "mtanh": 1},
        "mtanh": {"square": 2, "cube": 2, "mexp": 1, "mtanh": 0},
    },
)


rho_dense1D = np.linspace(0, 1.25, 1000)
rho_dense = rho_dense1D.reshape(-1, 1)
fig, ax = plt.subplots()
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$p(\rho)$')

v0=1.8393972058572117
Dr = 0.2
lp = v0/Dr
# test_name = "pfaps_test44_harmonic_measure_largeeps"
test_name = "pfaps_test50_vlowqsap_smalldt_measure"

# Direct pressure
print("Direct pressure")
data = np.loadtxt(test_name+'_sigmaIK_data', skiprows=1)
rhos = data[:,0]
sigmas = data[:,2]
sigmas_std = data[:,3]
marker = '.'
ax.errorbar(rhos, sigmas, yerr=sigmas_std, label=r"$p_{IK}$", marker=marker, ls='', ms=10, c='C0')

# perform symbolic regression fit
try:
    model_pC.fit(rhos.reshape(-1, 1) , sigmas, weights=1/np.square(sigmas_std))
    print(model_pC)
    equation = model_pC.latex()
    plt.plot(rho_dense1D, model_pC.predict(rho_dense), label=fr"${equation}$")
except RuntimeError:
    print(f"Could not fit to SigmaIK of {test_name}")

# Active pressure, which is converted to effective speed for symbolic fitting
print("Active pressure")
data = np.loadtxt(test_name+'_sigmaA_data', skiprows=1)
rhos = data[:,0]
indices = rhos<1.5
Pas = data[:,2]
Pas_std = data[:,3]
veffs = np.where(rhos > 0, Pas * 2 * Dr / rhos / v0, v0)
veffs_std = np.where(rhos > 0, Pas_std * 2 * Dr / rhos / v0, 0.00001)
ax.errorbar(rhos[indices], Pas[indices], yerr=Pas_std[indices], label=r"$p_A$", marker=marker, ls='', ms=10, c='C1')

# perform symbolic regression fit
try:
    model_pA.fit(rhos.reshape(-1, 1), veffs, weights=1/np.square(veffs_std))
    print(model_pA)
    equation = model_pA.latex()
    plt.plot(rho_dense1D, rho_dense1D*v0*model_pA.predict(rho_dense)/2/Dr, label=fr"$\rho * v_0 * ({equation}) / (2D_r)$", c='C1')
except RuntimeError:
    print(f"Could not fit to SigmaA of {test_name}")

ax.set_xlim(left=0)
# ax.set_ylim(bottom=0)
ax.set_title(f"Direct and active pressure for " + fr"$l_p = {lp}$")
ax.legend()
plt.savefig('pfap_harmonic_pressures_fit92_vlowqsap.png',  dpi=300, bbox_inches='tight')
