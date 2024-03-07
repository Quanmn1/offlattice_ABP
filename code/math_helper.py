import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import os
import sys

def convert_num(s):
    """
    If s is string of an integer (contains only numeric chars), return int(s)
    Else if s can be converted to a float, return float(s)
    Else, return s
    """
    if s.isnumeric():
        return int(s)
    try:
        return float(s)
    except ValueError:
        return s

def center_of_mass_pbc(x_lattice, Lx, one_dim_profile):
    """
    Find center of mass of a 1D system with PBC.
    Calculate CoM of the ring in yz plane, then project.
    """
    thetas = x_lattice/Lx*2*np.pi
    zs = np.cos(thetas)
    ys = np.sin(thetas)
    z_cm = np.average(one_dim_profile * zs)
    y_cm = np.average(one_dim_profile * ys)
    theta_cm = np.arctan2(y_cm, z_cm)
    return theta_cm/2/np.pi * Lx

def translate_pbc(x_lattice, Lx, distance):
    """
    The lattice is translated over 'distance'. This is used to move CoM into middle.
    The x_lattice array is reassigned, the density profile array remains the same.
    Same index, different position.
    """
    return (x_lattice+distance) % Lx