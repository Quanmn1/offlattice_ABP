import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import os
import sys

def read_param(param_file):
    with open(param_file, 'r') as f:
        params_info = f.readline().split()
        params_values = f.readline().split()
        params = {var:convert_num(value) for (var, value) in zip(params_info[2:], params_values[1:])}
        N = convert_num(f.readline().split()[0])
    return params, N

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

def coarsen(array, num_combined, binwidth=1, normalize=False):
    """
    Each element of the returned array is average of the next num_combined elements of array.
    If normalize=True, normalize the array to 1/binwidth (to fit with continuous distributions).
    """
    if num_combined == 1:
        if normalize:
            return array/np.sum(array)/binwidth
        else:
            return array
    original_len = len(array)
    result = np.zeros(original_len//num_combined)
    for i in range(len(result)):
        result[i] = np.average(array[i*num_combined:(i+1)*num_combined])
    if normalize:
        return result/np.sum(result)/binwidth
    else:
        return result
    
def nonzero(first, *arrs):
    ind = np.nonzero(first)
    return [first[ind]] + [arr[ind] for arr in arrs]

def bump(x, height, slope, mean, width):
    return np.where(np.abs(x-mean) < width/2, height * np.exp(slope * (1 - 1/(1 - 4/width**2 * (x-mean)**2)) ), 0)

def linear(x, a, b):
    return a*x + b