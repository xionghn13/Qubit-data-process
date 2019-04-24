import numpy as np
from matplotlib import pyplot as plt
from qutip import*
import QubitSpectrumFunc as qsf

from Single_small_junction import coupled_hamiltonian

directory = "E:\Projects\Fluxonium\data_process"
simulation = "Fluxonium022319"
path = directory + "\\" + simulation

#Qubit and computation parameters
Na = 30
Nr = 6
E_l = 0.30316786766768816
E_c = 1.3707449371055807
E_j = 5.081608341619772
wr = 7.332
g = 0.01
I0 = 2.41
I_period = 3.43

# E_l = 1.125
# E_c = 0.847
# E_j = 4.79
# wr = 7.3375
# g = 0.1


phi_ext = np.linspace(-0.2, -0.15, 101)
# phi_ext = np.linspace(0.58,0.584,11)
level_num = 10
energies = np.zeros((len(phi_ext),level_num))

########################################################################################
for idx, phi in enumerate(phi_ext):
    H = coupled_hamiltonian(Na, E_l, E_c, E_j, phi*2*np.pi, Nr, wr, g)
    for idy in range(level_num):
        energies[idx,idy] = abs(H.eigenenergies()[idy])
    qsf.printPercent(idx, len(phi_ext))
np.savetxt(path + '_energies.txt', energies)
########################################################################################
energies = np.genfromtxt(path+'_energies.txt')

#Plot eigensnergies
# fig1 = plt.figure(1)
# for idx in range(level_num):
#     plt.plot(phi_ext, energies[:,idx], linewidth = '2')
# plt.xlabel(r'$\varphi_\mathrm{ext}/2\pi$')
# plt.ylabel(r'Energy (GHz)')
# plt.ylim(top=30)

#Plot transition energies
# fig2 = plt.figure(2)
for idx in range(3, 5):
    plt.plot(I0 + I_period * phi_ext, energies[:,idx]-energies[:,0], linewidth = '2', color = 'b')
# plt.xlabel(r'$\varphi_\mathrm{ext}/2\pi$')
# plt.ylabel(r'$\mathrm{E_i} - \mathrm{E_0}$')
# plt.ylim([0,30])


plt.show()