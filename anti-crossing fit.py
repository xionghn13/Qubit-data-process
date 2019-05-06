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

# fig, ax = plt.subplots()
# for i in range(NumFile):
#     OnePhase = np.angle(OneComplex3List[i])
#     OneCurrentUniq = OneCurrentUniqList[i]
#     OneFreqUniq = OneFreqUniqList[i]
#     plt.plot(OneCurrentUniq, OneFreqUniq[np.argmax(abs(OnePhase), axis=0)], 'b')
# plt.xlabel('Current/mA', fontsize='x-large')
# plt.ylabel('freq/GHz', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')

clicked_data1 = np.array([
    [2.28, 7.227],
    [2.21, 7.177],
    [2.17, 7.1251],

])

clicked_data2 = np.array([

    [2.28, 7.303],
    [2.23, 7.2777],
    [2.13, 7.2607],
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

Na = 30
Nr = 5
I0_guess = 2.28
I_period_guess = 5.28
E_l = 0.4373
E_c = 2.8894
E_j = 6.9756
wr_guess = 7.2582
g_guess = .1622
FitData = True

guess = ([I0_guess, I_period_guess, g_guess, wr_guess])
bounds = (
    (I0_guess * 0.5, I_period_guess * 0.5, g_guess * 0.1, wr_guess * 0.9),
    (I0_guess * 1.5, I_period_guess * 1.5, g_guess * 10, wr_guess * 1.1)
)


def single_trans_energy(current_single, E_l, E_c, E_j, I0, I_period, g, wr, level_ind1, level_ind2):
    energy1 = np.zeros(len(current_single))  # predefined variable

    phi_ext1 = (current_single - I0) / I_period * 2 * np.pi
    a = tensor(destroy(Na), qeye(Nr))
    b = tensor(qeye(Na), destroy(Nr))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
    H_r = wr * (b.dag() * b + 1.0 / 2)
    H_c = g * na * (b.dag() + b)

    for idx in range(len(current_single)):
        ope1 = 1.0j * (phi - phi_ext1[idx])
        H_f = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope1.expm() + (-ope1).expm())
        H = H_f + H_r + H_c
        energy_levels1 = H.eigenenergies()
        energy1[idx] = energy_levels1[level_ind2] - energy_levels1[level_ind1]

    return energy1


##################Change levels here############################
def trans_energy_all(current, I0, I_period, g, wr):
    energy1 = single_trans_energy(current1, E_l, E_c, E_j, I0, I_period, g, wr, 0, 1)
    if len(clicked_data2) != 0:
        energy2 = single_trans_energy(current2, E_l, E_c, E_j, I0, I_period, g, wr, 0, 2)
        energy = np.concatenate([energy1, energy2], axis=0)
    else:
        energy = energy1
    print('trans_energy_all called ')
    print('I0=' + str(I0) + '\nI_period=' + str(I_period) + '\ng = ' + str(g) + '\nwr = ' + str(wr) + '\nE1 = ' + str(energy1[0]) + '\nE2 = ' + str(energy2[-1]))
    return energy


if FitData:
    opt, cov = curve_fit(trans_energy_all, current, freq, guess, bounds=bounds)
    I0_fit, I_period_fit, g_fit, wr_fit = opt
    print('I0 = ' + str(I0_fit) + '\nI_period = ' + str(I_period_fit) + '\ng = ' + str(g_fit) + '\nwr = ' + str(wr_fit))

# current_nice = np.linspace(Imin, Imax, 100)  # In mA

plt.plot(current, trans_energy_all(1, I0_guess, I_period_guess, g_guess, wr_guess), 'r')
if FitData:
    plt.plot(current, trans_energy_all(1, I0_fit, I_period_fit, g_fit, wr_fit), 'b')

plt.show()
