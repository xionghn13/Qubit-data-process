from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt
from qutip import *
import SubtractBackgroundFunc as sbf
import pickle

#####################################################################################
######################################Data###########################################
#####################################################################################
DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium032619/'
PickleFile = "SubtractBackgroundPickleDump.dat"

# Read data
with open(DataPath + PickleFile, "rb") as fp:  # Unpickling
    NoBackgroundDataList = pickle.load(fp)
    [OneCurrentUniqList, OneFreqUniqList, OneComplex3List] = NoBackgroundDataList

[[Imin, Imax], [freqmin, freqmax]] = sbf.dataListRange([OneCurrentUniqList, OneFreqUniqList, OneComplex3List])
NumFile = len(OneCurrentUniqList)

fig, ax = plt.subplots()
# for i in range(NumFile):
#     OnePhase = np.angle(OneComplex3List[i])
#     OneCurrentUniq = OneCurrentUniqList[i]
#     OneFreqUniq = OneFreqUniqList[i]
#     plt.plot(OneCurrentUniq, OneFreqUniq[np.argmax(abs(OnePhase), axis=0)], 'b')
# plt.xlabel('Current/mA', fontsize='x-large')
# plt.ylabel('freq/GHz', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')

clicked_data1 = np.array([
    [0.000978, 6.649072],
    [0.101601, 6.497885],
    [0.373282, 5.772188],
    [4.096316, 5.681476],
    [4.569241, 6.649072],
    [4.891234, 7.042157],
    [5.243412, 6.649072],
    [5.464782, 6.255986],
])

clicked_data2 = np.array([
    [2.285110, 7.737616],
    [2.486355, 7.949278],
    [2.878783, 8.281889],
    [3.160526, 8.372601],
    [4.106378, 8.130702],
    [4.901296, 7.525955],
    [5.374222, 7.888803],
])
current1 = clicked_data1[:, 0]  # In mA
freq1 = clicked_data1[:, 1]  # in GHz

if len(clicked_data2) != 0:
    current2 = clicked_data2[:, 0]  # In mA
    freq2 = clicked_data2[:, 1]  # in GHz

    current = np.concatenate([current1, current2], axis=0)
    freq = np.concatenate([freq1, freq2], axis=0)
else:
    current = current1
    freq = freq1
plt.plot(current, freq, 'o')  # plot mA

#####################################################################################
######################################Fit###########################################
#####################################################################################
# Define constants
e = 1.602e-19  # Fundamental charge
h = 6.62e-34  # Placnk's constant
phi_o = h / (2 * e)  # Flux quantum

N = 50
I0_guess = -0.36
I_period_guess = 5.28
E_l_guess = 0.45418260142466865
E_c_guess = 2.9
E_j_guess = 6.057730215694281

guess = ([E_l_guess, E_c_guess, E_j_guess, I0_guess, I_period_guess])


def single_trans_energy(current_single, E_l, E_c, E_j, I0, I_period, level_ind1, level_ind2):
    energy1 = np.zeros(len(current_single))  # predefined variable

    phi_ext1 = (current_single - I0) / I_period * 2 * np.pi
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
    for idx in range(len(current_single)):
        ope1 = 1.0j * (phi - phi_ext1[idx])
        H1 = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope1.expm() + (-ope1).expm())
        energy_levels1 = H1.eigenenergies()
        energy1[idx] = energy_levels1[level_ind2] - energy_levels1[level_ind1]

    return energy1


##################Change levels here############################
def trans_energy_all(current, E_l, E_c, E_j, I0, I_period):
    energy1 = single_trans_energy(current1, E_l, E_c, E_j, I0, I_period, 0, 1)
    if len(clicked_data2) != 0:
        energy2 = single_trans_energy(current2, E_l, E_c, E_j, I0, I_period, 0, 2)
        energy = np.concatenate([energy1, energy2], axis=0)
    else:
        energy = energy1
    print('trans_energy_all called ')
    print('EL=' + str(E_l) + '\nEC=' + str(E_c) + '\nEJ=' + str(E_j) +
          '\nI0=' + str(I0) + '\nI_period=' + str(I_period))
    return energy


# opt, cov = curve_fit(trans_energy_all, current, freq, guess)
# E_l_fit, E_c_fit, E_j_fit, I0_fit, I_period_fit = opt
# print('EL = ' + str(E_l_fit) + '\nEC = ' + str(E_c_fit) + '\nEJ = ' + str(E_j_fit) +
#       '\nI0 = ' + str(I0_fit) + '\nI_period = ' + str(I_period_fit))

# current_nice = np.linspace(Imin, Imax, 100)  # In mA

plt.plot(current, trans_energy_all(1, E_l_guess, E_c_guess, E_j_guess, I0_guess, I_period_guess), 'r')
# plt.plot(current, trans_energy_all(1, E_l_fit, E_c_fit, E_j_fit, I0_fit, I_period_fit), 'b')

plt.show()
