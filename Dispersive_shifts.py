import numpy as np
from matplotlib import pyplot as plt
import QubitSpectrumFunc as qsf

from Single_small_junction import charge_dispersive_shift as nChi

N = 50
E_l = 1
E_c = 2.3
E_j = 13
level_num = 30
g = 0.146
I0 = 5.27
hI = 2.57 # mA
w = 7.3


I_period = hI * 2

iState = 0
fState = 1

SingleCurrentPoint = 0.4
# plot dispersive shift as a function of flux
phi_ext1 = np.array([(SingleCurrentPoint - I0) / hI / 2])
phi_ext = np.linspace(0., 1, 51)
phi_ext = np.sort(np.concatenate((phi_ext1, phi_ext)))

chi = np.zeros(len(phi_ext))
# chi1 = np.zeros(len(phi_ext))
# chi2 = np.zeros(len(phi_ext))
# chi3 = np.zeros(len(phi_ext))
# kappa = 8.45 * 2  # MHz
for idx, phi in enumerate(phi_ext):
    chi[idx] = nChi(N, level_num, E_l, E_c, E_j, phi * 2 * np.pi, 0, 1, w, g)
    # chi1[idx] = pChi(N, level_num, E_l, E_c, E_j, phi*2*np.pi, 2, 0, w, g)
    # chi2[idx] = nChi(N, level_num, E_l, E_c, E_j, phi * 2 * np.pi, 2, 0, w, g)
    # chi3[idx] = nChi(N, level_num, E_l, E_c, E_j, phi * 2 * np.pi, 3, 0, w, g)
    qsf.printPercent(idx, len(phi_ext))
# chi is in GHz
# chi_angle is in degree
Ind0 = phi_ext == 0
chi0 = chi[Ind0] * 1e3
Indh = phi_ext == 0.5
chih = chi[Indh] * 1e3
IndS = phi_ext == phi_ext1[0]
chiS = chi[IndS] * 1e3
print('At 0 flux, chi/2pi=%.3GMHz. At half flux quanta, chi/2pi=%.3GMHz.' % (chi0, chih))
print('At %.3GmA, chi/2pi=%.3GMHz.' % (SingleCurrentPoint, chiS))
# chi_angle = chi*1e3/(kappa/2) *180/np.pi
# chi1_angle = chi1*1e3/(kappa/2) *180/np.pi
# chi2_angle = chi2*1e3/(kappa/2) *180/np.pi
# chi3_angle = chi3*1e3/(kappa/2) *180/np.pi

plt.figure(1)
plt.plot(phi_ext * I_period + I0, chi * 1e3, 'k-')
# plt.plot(phi_ext, chi1*1e3 , 'b-')
# plt.plot(phi_ext, chi2*1e3 , 'r-')
# plt.plot(phi_ext, chi3*1e3 , 'g-')
# plt.grid()
# plt.ylim([-2,2])
# plt.tick_params(labelsize = 18.0)

# plt.figure(2)
# plt.plot(phi_ext, chi_angle , 'k-')
# plt.plot(phi_ext, chi1_angle , 'b-')
# plt.plot(phi_ext, chi2_angle , 'r-')
# plt.plot(phi_ext, chi3_angle , 'g-')
# plt.ylim([-20,20])
# plt.tick_params(labelsize = 18.0)
# chi_angle = chi1*1e3/(kappa/2) *180/np.pi
# plt.plot(phi_ext, chi_angle , 'r.')
# chi_angle = chi2*1e3/(kappa/2) *180/np.pi
# plt.plot(phi_ext, chi_angle , 'g.')

# plot dispersive shift as a function of cavity frequency
# phi_ext = 0.5
# w = np.linspace(3,12,901)
# chi = np.zeros(len(w))
# kappa = 5 #MHz
# #
iState = 0
fState = 1
# for idx, freq in enumerate(w):
#     chi[idx]= nChi(N, level_num, E_l, E_c, E_j, phi_ext*2*np.pi, iState, fState, freq, g)

# chi_angle = chi*1e3/(kappa/2) *180/np.pi
# plt.plot(w, chi*1e3, '.')
#
# iState = 0
# fState = 2
# for idx, freq in enumerate(w):
#     chi[idx]= nChi(N, level_num, E_l, E_c, E_j, phi_ext*2*np.pi, iState, fState, freq, g)
#
# chi_angle = chi*1e3/(kappa/2) *180/np.pi
# plt.plot(w, chi*1e3, '-')
# plt.ylim([-2,2])
# iState = 1
# fState = 2
# for idx, freq in enumerate(w):
#     chi[idx]= nChi(N, level_num, E_l, E_c, E_j, phi_ext*2*np.pi, iState, fState, freq, g)
#
# chi_angle = chi*1e3/(kappa/2) *180/np.pi
# plt.plot(w, chi_angle, '.')
#
# plt.ylim([-10,10])
plt.xlabel('Current (mA)')
plt.ylabel('Chi_01/2pi (MHz)')
plt.grid()
plt.show()
