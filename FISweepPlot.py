import numpy as np
import ExtractDataFunc as edf
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
import pickle
import pprofile


# Click on the points on screen to define an approximation line for interpolation
def onclick(event):
    print('[%f, %f],' % (event.xdata, event.ydata))


profiler = pprofile.Profile()

# DataPath = 'E:/Projects\Fluxonium\data_process/test_data/'
# BackgroundFile = '081717_one_tone_4GHz_to_12GHz_-20dBm.dat'
# OneToneFile = '081817_one_tone_6.4GHz_to_9.0GHz_-2.5mA_to_2.5mA_1.dat'
DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium032619/'
BackgroundFile = 'one_tone_6.45GHz_to_6.6GHz_5dBm_0mA_10us integration_200Kavg_100KHz step_091918.dat'
# BackgroundFile = 'one_tone_3.5GHz_to_12GHz_5dBm_3mA_10us integration_5Kavg_500KHz step_121018.dat'
# OneToneFile = '110518_one_tone_4GHz_to_9GHz_-2mA_to_3mA_10us integration_3Kavg.dat'
# OneToneFileList = ['110518_one_tone_4GHz_to_9GHz_-2mA_to_3mA_10us integration_3Kavg.dat',
#                    # OneToneFile = '110618_one_tone_4GHz_to_9GHz_3mA_to_10mA_10us integration_3Kavg.dat'
#                    '110618_one_tone_4GHz_to_9GHz_3mA_to_10mA_10us integration_3Kavg.dat',
#                    '110718_one_tone_7.5GHz_to_9GHz_4.5mA_to_5.5mA_10us integration_3Kavg.dat',
#                    '110718_one_tone_4GHz_to_9GHz_10mA_to_15mA_10us integration_3Kavg.dat']
# OneToneFileList = ['111918_one_tone_8GHz_to_9GHz_-1mA_to_1mA_10us integration_3Kavg.dat',
#                    '111918_one_tone_8GHz_to_9GHz_1mA_to_2mA_10us integration_3Kavg.dat',
#                    '111918_one_tone_7GHz_to_8GHz_1.5mA_to_3mA_10us integration_3Kavg.dat',
#                    '111918_one_tone_6GHz_to_7GHz_2.5mA_to_3.4mA_10us integration_3Kavg.dat',
#                    '112018_one_tone_5GHz_to_6GHz_3.1mA_to_3.7mA_10us integration_3Kavg.dat',
#                    '112018_one_tone_4GHz_to_5GHz_3.5mA_to_4.1mA_10us integration_3Kavg.dat',
#                    '112618_one_tone_3GHz_to_4GHz_4mA_to_4.5mA_10us integration_3Kavg.dat',
#                    '112518_one_tone_4.6GHz_to_6.2GHz_4mA_to_6mA_10us integration_5Kavg.dat']
# OneToneFileList = [
#     '010319_one_tone_3.5GHz_to_6GHz_0mA_to_4mA_10us integration_3Kavg.dat',
#     '010319_one_tone_3.5GHz_to_5GHz_4.2mA_to_7mA_10us integration_3Kavg.dat',
#     '010319_one_tone_4.5GHz_to_5.3GHz_0mA_to_3mA_10us integration_3Kavg.dat',
#     '010319_one_tone_3GHz_to_4.5GHz_3mA_to_4mA_10us integration_3Kavg.dat',
#     '010319_one_tone_3GHz_to_4.3GHz_4mA_to_7mA_10us integration_3Kavg.dat',
#     '010419_one_tone_4.5GHz_to_5.3GHz_8mA_to_10mA_10us integration_3Kavg.dat',
#     '010519_one_tone_4.07GHz_to_4.12GHz_5.3mA_to_5.5mA_10us integration_3Kavg.dat',
#     '010519_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.4mA_10us integration_30Kavg.dat',
#     '010819_one_tone_5.25GHz_to_5.35GHz_5.0mA_to_5.6mA_10us integration_3Kavg.dat',
#     '010919_one_tone_4.07GHz_to_4.12GHz_5.1mA_to_5.5mA_10us integration_30Kavg.dat',
#     '010919_one_tone_5.2GHz_to_6.3GHz_3mA_to_6mA_10us integration_3Kavg.dat',
#  ]
# OneToneFileList = [
#     '011019_one_tone_4.5GHz_to_5.3GHz_0mA_to_3mA_10us integration_3Kavg.dat',
#     '011019_one_tone_3GHz_to_4.5GHz_3mA_to_4mA_10us integration_3Kavg.dat',
#     '011019_one_tone_3GHz_to_4.3GHz_4mA_to_5.5mA_10us integration_3Kavg.dat',
#     '011119_one_tone_4.07GHz_to_4.12GHz_5.3mA_to_5.5mA_10us integration_3Kavg.dat',
#     # '011119_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.3mA_10us integration_3Kavg.dat',
#     '011119_one_tone_8.2GHz_to_9.4GHz_5.1mA_to_5.4mA_10us integration_3Kavg.dat',
#     # '011119_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.3mA_10us integration_3Kavg_2.dat',
#     '011119_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.4mA_10us integration_3Kavg.dat',
#     '011319_one_tone_8.6GHz_to_8.9GHz_4.9mA_to_5.1mA_10us integration_3Kavg.dat',
#
#
#  ]

# OneToneFileList = [
#     '011519_one_tone_9.54GHz_to_9.66GHz_3.4mA_to_3.5mA_10us integration_3Kavg.dat',
#     '011519_one_tone_3.25GHz_to_3.42GHz_3.4mA_to_3.5mA_10us integration_3Kavg.dat',
#     '011619_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.4mA_10us integration_3Kavg.dat',
# ]
# OneToneFileList = [
#     '012319_one_tone_4.5GHz_to_5.3GHz_0mA_to_3mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3GHz_to_4.5GHz_3mA_to_4mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.3GHz_to_3.4GHz_3.5mA_to_3.6mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.32GHz_to_3.35GHz_3.52mA_to_3.56mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.35GHz_to_3.37GHz_3.52mA_to_3.56mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.32GHz_to_3.357GHz_3.42mA_to_3.46mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.34GHz_to_3.46GHz_3.4mA_to_3.5mA_10us integration_3Kavg.dat',
#     '012319_one_tone_3.32GHz_to_3.35GHz_3.46mA_to_3.52mA_10us integration_3Kavg.dat',
#     '012319_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.4mA_10us integration_3Kavg.dat',
# ]

# OneToneFileList = [
#     'one tone_3.hdf5',
#     'one tone_4_0.hdf5',
#     'one tone_4.hdf5',
#     # '012819_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.4mA_10us integration_3Kavg.dat',
#     '012919_one_tone_8.2GHz_to_9.4GHz_5.2mA_to_5.4mA_10us integration_3Kavg.dat',
#     '020419_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.5mA_10us integration_3Kavg.dat',
#
# ]

# OneToneFileList = [
#     # '021119_one_tone_3.8GHz_to_4.8GHz_2.2mA_to_3.2mA_10us integration_3Kavg.dat',
#     # '021119_one_tone_4.07GHz_to_4.12GHz_5.2mA_to_5.5mA_10us integration_3Kavg.dat',
#     '021119_one_tone_3.8GHz_to_4.8GHz_0mA_to_2.2mA_10us integration_3Kavg.dat',
#     '021119_one_tone_3.8GHz_to_4.8GHz_0.2mA_to_0.5mA_10us integration_3Kavg.dat',
#     '021119_one_tone_4.8GHz_to_5.2GHz_0.0mA_to_0.2mA_10us integration_3Kavg.dat',
#     '021119_one_tone_3.4GHz_to_3.8GHz_0.45mA_to_0.55mA_10us integration_3Kavg.dat',
#     '021119_one_tone_3.8GHz_to_4.2GHz_0.8mA_to_1.1mA_10us integration_3Kavg.dat',
#     # '021119_one_tone_3.99GHz_to_4.05GHz_0.88mA_to_0.92mA_10us integration_6Kavg.dat',
#
# ]

# OneToneFileList = [
#     '021319_one_tone_4.7GHz_to_5.2GHz_0mA_to_0.3mA_10us integration_3Kavg.dat',
#     '021319_one_tone_3.4GHz_to_4.8GHz_0.2mA_to_0.55mA_10us integration_3Kavg.dat',
#     '021319_one_tone_3.8GHz_to_4.2GHz_0.8mA_to_1.1mA_10us integration_3Kavg.dat',
#     '021519_one_tone_8.45GHz_to_8.9GHz_0.8mA_to_1.0mA_10us integration_3Kavg.dat',
#     '021519_one_tone_10GHz_to_12GHz_0.8mA_to_1.0mA_10us integration_3Kavg.dat',
#     '021919_one_tone_5.2GHz_to_5.5GHz_0.8mA_to_0.9mA_10us integration_3Kavg.dat',
#     'two tone_2.hdf5',
#     'two tone_3.hdf5',
# ]
#

#
# OneToneFileList = [
#     # 'one tone_59.hdf5',
#     'one tone_57.hdf5',
#     'one tone_58.hdf5',
#
#     'two tone_28.hdf5',
#     'two tone_29.hdf5',
#     'two tone_30.hdf5',
#     'two tone_31.hdf5',
#     'two tone_32.hdf5',
#     'two tone_33.hdf5',
#     'two tone_34.hdf5',
#     'two tone_35.hdf5',
#     'two tone_36.hdf5',
#     'two tone_37.hdf5',
#     'two tone_40.hdf5',
#     'two tone_64.hdf5',
#     'two tone_66.hdf5',
#     'two tone_69.hdf5',
#     'two tone_70.hdf5',
#     'two tone_71.hdf5',
#     'two tone_91.hdf5',
#     'two tone_92.hdf5',
# ]

OneToneFileList = [
    # 'two tone_113.hdf5',
    'one tone_72.hdf5',
    # 'one tone_75.hdf5',
    # 'two tone_114.hdf5',
    # 'two tone_115.hdf5',
    # 'two tone_116.hdf5',
    # 'two tone_117.hdf5',
    # 'two tone_118.hdf5',
    # 'two tone_119.hdf5',
    # 'two tone_120.hdf5',
    # 'two tone_122.hdf5',
    # 'two tone_123.hdf5',
    # 'two tone_124.hdf5',
    # 'two tone_125.hdf5',
    # 'two tone_126.hdf5',
    # 'two tone_127.hdf5',
    # 'two tone_141.hdf5',
    # 'two tone_142.hdf5',
    # 'two tone_144.hdf5',
    # 'two tone_145.hdf5',
    # 'two tone_146.hdf5',
    # 'two tone_149.hdf5',

]
# OneToneFileList = [
#     'one tone_66.hdf5'
# ]
PlotSpectrum = False
ClickForPoints = False
PickleSave = True
NoCalibrate = False  # for one tone at anti crossing
SelfCalibrate = True  # use background file to calibrate

RawDataIndexForPlot = 0
PickleFile = "SubtractBackgroundPickleDump.dat"
OneTonePower = -15
# EJ = 4.751
# EC = 2.398
# EL = 1.286
# I0 = 0.15
# I_period = 10.1
# EL = 0.5291959920767306
# EC = 1.949345345953156
# EJ = 3.000827450217155
# # I0 = 1.2406555052048085
# # I_period = 8.140664513370277
# I0 = 0.05
EL = 0.30316786766768816
EC = 1.3707449371055807
EJ = 5.081608341619772
# I0 = 1.662483629052521
# I_period = 5.060322002288046
I0 = 0.9
I_period = 3
N = 50
level_num = 7
########################################################################################################################
NumFile = np.size(OneToneFileList)
[OneCurrentUniqList, OneFreqUniqList, OneComplex3List] = [[[] for i in range(NumFile)] for j in range(3)]
for i in range(NumFile):
    if NoCalibrate:
        [OneFreq, OneCurrent, OneComplex3List[i]] = edf.readFISweepLabber(DataPath + OneToneFileList[i])
        OneFreqUniqList[i] = np.unique(OneFreq)
        OneCurrentUniqList[i] = np.unique(OneCurrent)
    elif SelfCalibrate:
        [OneFreqUniqList[i], OneCurrentUniqList[i], OneComplex3List[i]] = sbf.FISweepSelfCalibrate(DataPath,
                                                                                                   OneToneFileList[i])
    else:
        [OneFreqUniqList[i], OneCurrentUniqList[i], OneComplex3List[i]] = sbf.FISweepBackgroundCalibrate(DataPath,
                                                                                                         OneToneFileList[
                                                                                                             i],
                                                                                                         BackgroundFile,
                                                                                                         OneTonePower)
NoBackgroundDataList = [OneCurrentUniqList, OneFreqUniqList, OneComplex3List]
[[Imin, Imax], [freqmin, freqmax]] = sbf.dataListRange(NoBackgroundDataList)

########################################################################################################################
if PlotSpectrum:
    I_array = np.linspace(Imin, Imax, 100)
    with profiler:
        [[freq0x, freq1x, freq2x], [phi0x, phi1x, phi2x], [n0x, n1x, n2x]] = qsf.calcSpectrum(level_num, N, EL, EC, EJ,
                                                                                              I0,
                                                                                              I_period, I_array)
    freqList = [freq0x, freq1x, freq2x]
    print('finish simulating spectrum')
########################################################################################################################
if OneToneFileList[RawDataIndexForPlot].endswith('.dat'):
    [OneFreqRaw, OneCurrentRaw, OneComplexRaw] = edf.readFISweepDat(DataPath + OneToneFileList[RawDataIndexForPlot])
elif OneToneFileList[RawDataIndexForPlot].endswith('.hdf5'):
    if OneToneFileList[RawDataIndexForPlot].startswith('one tone'):
        [OneFreqRaw, OneCurrentRaw, OneComplexRaw] = edf.readFISweepLabber(
            DataPath + OneToneFileList[RawDataIndexForPlot])
    else:
        [OneFreqRaw, OneCurrentRaw, OneComplexRaw] = edf.readFISweepTwoToneLabber(
            DataPath + OneToneFileList[RawDataIndexForPlot])
OneFreqRawUniq = np.unique(OneFreqRaw)
OnePhaseRaw = np.unwrap(np.angle(OneComplexRaw), axis=0)
if len(OnePhaseRaw) > 1000:
    first_ind = 100
else:
    first_ind = 0
PhaseSlope = np.mean((OnePhaseRaw[-1] - OnePhaseRaw[first_ind]) / (OneFreqRawUniq[-1] - OneFreqRawUniq[first_ind]))

if not SelfCalibrate:
    [BackFreqRaw, BackComplexRaw] = edf.readFSweepDat(DataPath + BackgroundFile)
    BackPhaseRaw = np.unwrap(np.angle(BackComplexRaw))

    leg = ('One tone -15dBm', 'Background 5dBm')
# fig, ax = plt.subplots()
# plt.plot(OneFreqRawUniq, np.abs(OneComplexRaw[:, 0]) * 10 ** (15 / 20))
# plt.plot(BackFreqRaw, np.abs(BackComplexRaw) * 10 ** (- 5 / 20), '--')
# plt.xlabel('freq/GHz', fontsize='x-large')
# plt.ylabel('Abs', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')
# plt.legend(leg)
# plt.tight_layout()
#
# fig, ax = plt.subplots()
# plt.plot(OneFreqRawUniq,
#          np.transpose(OnePhaseRaw.transpose() - PhaseSlope * (OneFreqRawUniq - OneFreqRawUniq[100]))[:, 0])
# plt.plot(BackFreqRaw,
#          np.transpose(BackPhaseRaw.transpose() - PhaseSlope * (BackFreqRaw - BackFreqRaw[2000]) + 4 * np.pi))
# plt.xlabel('freq/GHz', fontsize='x-large')
# plt.ylabel('Phase', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')
# plt.legend(leg)
#
# plt.tight_layout()

if PickleSave:
    with open(DataPath + PickleFile, "wb") as fp:  # Pickling
        pickle.dump(NoBackgroundDataList, fp)

# plot amplitude

fig, ax = plt.subplots()
for i in range(NumFile):
    OneAmp = np.abs(OneComplex3List[i])
    OneCurrentUniq = OneCurrentUniqList[i]
    OneFreqUniq = OneFreqUniqList[i]
    pcm = plt.pcolormesh(OneCurrentUniq, OneFreqUniq, OneAmp)
    # pcm = plt.pcolormesh(OneCurrentUniq, OneFreqUniq, OneAmp, vmax=1, vmin=0)

cbar = plt.colorbar(pcm)
cbar.ax.set_ylabel('Abs(r)', fontsize='x-large')
if PlotSpectrum:
    qsf.plotTransitions(I_array, freqList)
plt.ylim(freqmin, freqmax)
plt.xlabel('Current/mA', fontsize='x-large')
plt.ylabel("Freq/GHz", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if ClickForPoints:
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

# plot phase

fig, ax = plt.subplots()

for i in range(NumFile):
    OnePhase = np.angle(OneComplex3List[i])
    if i == 0:
        MaxPhase = OnePhase.max()
        MinPhase = OnePhase.min()
    else:
        MaxPhase = max(OnePhase.max(), MaxPhase)
        MinPhase = min(OnePhase.min(), MinPhase)

for i in range(NumFile):
    OnePhase = np.angle(OneComplex3List[i])
    OneCurrentUniq = OneCurrentUniqList[i]
    OneFreqUniq = OneFreqUniqList[i]
    # pcm = plt.pcolormesh(OneCurrentUniq, OneFreqUniq, OnePhase, vmax=MaxPhase, vmin=MinPhase)
    pcm = plt.pcolormesh(OneCurrentUniq, OneFreqUniq, OnePhase)

# cbar = plt.colorbar(pcm)
# cbar.ax.set_ylabel('Phase', fontsize='x-large')
if PlotSpectrum:
    qsf.plotTransitions(I_array, freqList)
plt.ylim(freqmin, freqmax)
plt.xlabel('Current/mA', fontsize='x-large')
plt.ylabel("Freq/GHz", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if ClickForPoints:
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
#
# fig, ax = plt.subplots()
# for i in range(NumFile):
#     OnePhase = np.angle(OneComplex3List[i])
#     OneCurrentUniq = OneCurrentUniqList[i]
#     OneFreqUniq = OneFreqUniqList[i]
#     plt.plot(OneCurrentUniq, OneFreqUniq[np.argmax(abs(OnePhase), axis=0)], 'b')
# if PlotSpectrum:
#     qsf.plotTransitions(I_array, freqList)
#     # for i in range(level_num - 1):
#     #     ax.plot(I_array, freq0x[:, i], ':r')
#     # for i in range(level_num - 2):
#     #     ax.plot(I_array, freq1x[:, i], ':g')
#     # for i in range(level_num - 3):
#     #     ax.plot(I_array, freq2x[:, i], ':k')
# plt.ylim(freqmin, freqmax)
# plt.xlabel('Current/mA', fontsize='x-large')
# plt.ylabel('freq/GHz', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')
# plt.tight_layout()

plt.show()
