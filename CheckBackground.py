import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
import ExtractDataFunc as edf

DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium123118/'
# BackgroundFile = 'one_tone_7GHz_to_8GHz_5dBm_2.5mA_10us integration_50Kavg_50KHz step_121918.dat'
BackgroundFileList = [
    'one_tone_3.9GHz_to_4.2GHz_-15dBm_0.8mA_10us integration_100Kavg_50KHz step_021319.dat',
    '021519_one_tone_3.9GHz_to_4.2GHz_0.89mA_to_0.89mA_pumped_1.312GHz_20dBm30kavg.dat',
    '021519_one_tone_3.9GHz_to_4.2GHz_0.89mA_to_0.89mA_pumped_1.312GHz_20dBm10kavg.dat',
    '021519_one_tone_3.9GHz_to_4.2GHz_0.89mA_to_0.89mA_pumped_1.312GHz_-20dBm10kavg.dat',
]
LegendList = [
    'no pump qubit at 3.94GHz',
    '20dBm longer time',
    '20dBm',
    '-20dBm',
]

UseBackgroundRange = True

NumFile = len(BackgroundFileList)
[BackFreqList, BackAmpList, BackAmpNormalizedList, BackPhaseList] = [[[] for i in range(NumFile)] for j in range(4)]
FreqLimit = [7.6, 7.9]
PhaseSlope = 0
for i in range(NumFile):
    if BackgroundFileList[i].startswith('one_tone'):
        Background = np.loadtxt(DataPath + BackgroundFileList[i], skiprows=1)
        BackFreqList[i] = Background[:, 0]
        BackFreq = BackFreqList[i]
        BackPhase = np.unwrap(Background[:, 1])
        BackPowerStr = BackgroundFileList[i].split('_')[5][:-3]
        BackPower = float(BackPowerStr)
        BackAmpList[i] = Background[:, 2]
    else:
        [OneFreq, OneCurrent, OneComplex] = edf.readFISweepDat(DataPath + BackgroundFileList[i])
        BackFreq = OneFreq[:, 0]
        BackFreqList[i] = BackFreq
        # BackPowerStr = BackgroundFileList[i].split('_')[11][:-7]
        BackPower = float(-15)
        BackAmpList[i] = np.abs(OneComplex[:, 0])
        BackPhase = np.unwrap(np.angle(OneComplex[:, 0]))

    BackAmpNormalizedList[i] = BackAmpList[i] * 10 ** (- BackPower / 20)
    if PhaseSlope == 0:
        PhaseSlope = (BackPhase[-1] - BackPhase[0]) / (BackFreq[-1] - BackFreq[0])
    # BackPhaseList[i] = BackPhase
    BackPhaseList[i] = (BackPhase - PhaseSlope * (BackFreq - 4.105)) % (2 * np.pi)

fig, ax = plt.subplots()
for i in range(NumFile):
    # if i == 1:
    #     BackFreqList[i] = BackFreqList[i] - 0.1
    plt.plot(BackFreqList[i], BackAmpNormalizedList[i])
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Abs', fontsize='x-large')
plt.legend(LegendList)
if not UseBackgroundRange:
    plt.xlim(FreqLimit[0], FreqLimit[1])
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(BackFreqList[i], BackPhaseList[i])
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Phase', fontsize='x-large')
plt.legend(LegendList)
if not UseBackgroundRange:
    plt.xlim(FreqLimit[0], FreqLimit[1])
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
