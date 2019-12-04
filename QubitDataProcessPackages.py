import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
import QubitDecayFunc
import ExtractDataFunc as edf
import h5py
import time
import calendar

e = 1.60217662e-19 # C
hbar = 1.0545718e-34 # m^2 * kg / s
h = 6.62607004e-34 # m^2 * kg / s
R_K = 25812.807557 # Ohm
Phi0 = 2.067833831e-15 # Wb
phi0 = Phi0 / (2. * np.pi) # Wb
k_B = 1.38064852e-23 # m^2 * kg / (s^2 * K)

C_J_unit_area = 45.e-15 # F / um^2
C_J_unit_area_fF = 45. # fF / um^2

Delta_Al_eV = 0.18e-3 # eV
Delta_Al = Delta_Al_eV * e