import numpy as np
from matplotlib import pyplot as plt
from qutip import*
import QubitSpectrumFunc as qsf

from Single_small_junction import coupled_hamiltonian

directory = "E:\Projects\Fluxonium\data_process"
simulation = "Fluxonium042619"
path = directory + "\\" + simulation

#Qubit and computation parameters
Na = 30
Nr = 6
E_l = 1.016
E_c = 1.421
E_j = 14.858
wr = 7.45
g = 0.
I0 = 0
I_period = 1

# E_l = 1.125
# E_c = 0.847
# E_j = 4.79
# wr = 7.3375
# g = 0.1


phi_ext = np.linspace(0.297, 0.298, 51)
# phi_ext = np.linspace(0.58,0.584,11)
level_num = 10
energies = np.zeros((len(phi_ext),level_num))

########################################################################################
for idx, phi in enumerate(phi_ext):
    H = coupled_hamiltonian(Na, E_l, E_c, E_j, phi*2*np.pi, Nr, wr, g)
    for idy in range(level_num):
        energies[idx,idy] = H.eigenenergies()[idy]
    qsf.printPercent(idx, len(phi_ext))
np.savetxt(path + '_energies.txt', energies)
########################################################################################
energies = np.genfromtxt(path+'_energies.txt')

#Plot eigensnergies
fig1 = plt.figure(1)
for idx in range(level_num):
    if idx == 0:
        plt.plot(phi_ext, energies[:, idx], linewidth='2', color = 'k')
    else:
        plt.plot(phi_ext, energies[:,idx], linewidth = '2')
plt.xlabel(r'$\varphi_\mathrm{ext}/2\pi$')
plt.ylabel(r'Energy (GHz)')
plt.ylim(top=30)

#Plot transition energies
fig2 = plt.figure(2)
for idx in range(1, 3):
    plt.plot(I0 + I_period * phi_ext, energies[:,idx]-energies[:, 0], linewidth = '2', color = 'b')
plt.xlabel(r'$\varphi_\mathrm{ext}/2\pi$')
plt.ylabel(r'$\mathrm{E_i} - \mathrm{E_0}$')
# plt.ylim([0,30])


plt.show()