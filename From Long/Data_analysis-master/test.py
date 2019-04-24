import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from qutip import*
plt.figure(figsize=(10, 10))
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')

directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_38.56to38.66mA_Qubit_4.2to5.1GHz_-6dBm_Cav_10.3045GHz_5dBm_IF_0.05GHz_measTime_500ns_avg_20000'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1:-1] - 0.037
freq = np.genfromtxt(path + '_FREQ.dat')
freq = freq[1::]
data = np.genfromtxt(path + '_PHASEMAG.dat')
phase = data[1::,0] #phase is recorded in rad
phase = phase#
mag = data[1::,0]

Z = np.zeros((len(current),len(freq)))
for idx in range(len(current)):
    temp = np.unwrap(phase[idx*len(freq):(idx+1)*len(freq)])
    Z[idx,:] = temp - np.average(temp)
Z = Z*180/(np.pi)
Z = Z.transpose()[1:len(freq)-1]
X,Y = np.meshgrid(current,freq[1:len(freq)-1])
plt.figure(1)
plt.pcolormesh(X,Y,Z, cmap= 'GnBu_r', vmin = -4, vmax=-0.5)
#####################################################################################################################################################################################
#####################################################################################################################################################################################

# Define constants
e = 1.602e-19  # Fundamental charge
h = 6.62e-34  # Placnk's constant
phi_o = h / (2 * e)  # Flux quantum

"""
First section of the script attempts to plot the energies vs external flux
"""


# Hamiltonian definition
def Ho(N, E_l, E_c, E_j_sum, d, phi_squid, phi_ext):
    E_j1 = 0.5 * E_j_sum * (1 + d)
    E_j2 = 0.5 * E_j_sum * (1 - d)
    a = tensor(destroy(N))
    mass = 1.0 / (8.0 * E_c)
    w = sqrt(8.0 * E_c * E_l)
    phi = (a + a.dag()) * (8 * E_c / E_l) ** (0.25) / np.sqrt(2)
    na = 1j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2)
    ope1 = 1j * (-phi + phi_ext)
    ope2 = 1j * (
    phi - phi_ext + phi_squid)  # phi_squid and phi_ext here are the external phases, or normalized flux, = flux*2pi/phi_o
    H = 4.0 * E_c * na ** 2 + 0.5 * E_l * (phi) ** 2 - 0.5 * E_j1 * (ope1.expm() + (-ope1).expm()) - 0.5 * E_j2 * (
    ope2.expm() + (-ope2).expm())
    return H.eigenenergies()


def coupled_H(Na, E_l, E_c, E_j_sum, d, phi_squid, phi_ext, Nr, wr, g):
    E_j1 = 0.5 * E_j_sum * (1 + d)
    E_j2 = 0.5 * E_j_sum * (1 - d)
    a = tensor(destroy(Na), qeye(Nr))
    b = tensor(qeye(Na), destroy(Nr))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
    ope1 = 1.0j * (phi_ext - phi)
    ope2 = 1.0j * (phi + phi_squid - phi_ext)
    H_f = 4.0 * E_c * na ** 2 + 0.5 * E_l * (phi) ** 2 - 0.5 * E_j1 * (ope1.expm() + (-ope1).expm()) - 0.5 * E_j2 * (
    ope2.expm() + (-ope2).expm())
    H_r = wr * (b.dag() * b + 1.0 / 2)
    H_c = -2 * g * na * (b.dag + b)
    H = H_f + H_r + H_c
    return H.eigenenergies()


def trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState):
    B_field = current * B_coeff * 1e-4  # in T, this depends on a seperate measurement
    phi_squid = B_field * A_j  # these are flux, not normalized
    phi_ext = B_field * A_c
    trans_energy = np.zeros((level_num - iState, len(phi_ext)))
    for idx in range(len(phi_ext)):
        energies = Ho(N, E_l, E_c, E_j_sum, d, 2 * np.pi * (phi_squid[idx] / phi_o - beta_squid),
                      2 * np.pi * (phi_ext[idx] / phi_o - beta_ext))  # normalize the flux -> phase here
        for level in range(iState + 1, level_num):
            trans_energy[level - iState, idx] = energies[level] - energies[iState]
    return trans_energy


def coupled_trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState,
                           Nr, wr, g):
    B_field = current * B_coeff * 1e-4  # in T, this depends on a seperate measurement
    phi_squid = B_field * A_j  # these are flux, not normalized
    phi_ext = B_field * A_c
    trans_energy = np.zeros((level_num - iState, len(phi_ext)))
    for idx in range(len(phi_ext)):
        energies = coupled_H(N, E_l, E_c, E_j_sum, d, 2 * np.pi * (phi_squid[idx] / phi_o - beta_squid),
                             2 * np.pi * (phi_ext[idx] / phi_o - beta_ext), Nr, wr,
                             g)  # normalize the flux -> phase here
        for level in range(iState + 1, level_num):
            trans_energy[level - iState, idx] = energies[level] - energies[iState]
    return trans_energy


########################################################################
# Fitting for bottom spectrum
N = 50
E_l=0.722729827116
E_c=0.552669197076
E_j_sum=17.61374383
A_j=4.76321410213e-12
A_c=1.50075181762e-10
d=0.125005274368
beta_squid=0.129912406349
beta_ext=0.356925557542

current = np.linspace(0.0385, 0.0387, 100)
B_coeff = 60
level_num = 5

iState = 0
spectrum = trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState)
for idx in range(iState, level_num):
    line = plt.plot(current * 1e3, spectrum[idx, :])  # transition from state (iState)
    plt.setp(line, linewidth=3.0, linestyle='--', color="black")#, alpha=0.5)

ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
plt.xlim([38.523, 38.606])
plt.ylim([4.2,5.1])
plt.xticks([38.55, 38.6])
plt.yticks([4.2, 4.6, 5.0])
plt.tick_params(labelsize=26)
plt.show()