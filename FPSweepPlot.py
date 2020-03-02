import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
import ExtractDataFunc as edf

# DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium022319/'
DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
# BackgroundFile = 'calibration_5.hdf5'
# BackgroundFile = 'power spectroscopy_76.hdf5'
# OneToneFile = 'power spectroscopy_82.hdf5'
BackgroundFile = 'power spectroscopy_138.hdf5'
OneToneFile = 'power spectroscopy_140.hdf5'
# BackgroundFile = 'power spectroscopy_76.hdf5'
# OneToneFile = 'power spectroscopy_77.hdf5'

Calibration = True
UseOneToneRange = False
FitSeparately = False
PlotParamVSPower = False
PlotUnfittedCircle = False
ShiftCircle = False
RotateCircle = False
LineSpec = '.'

StartFreq = 6.538
EndFreq = 6.552
StartPower = -10
EndPower = 5
SelectPower = np.array([])
# SelectPower = np.array([-10, -5, 2.5, 5])

gamma_f_guess = 1.9027e-3
P0_guess = 0.86229
A_guess = 3e3
amp_cor_re_guess = 1
# amp_cor_im_guess = -0.3
P0_im_guess = 0.

# bounds = (
#     [1, 1.9027e-4 - 1e-12, 0.86229 - 1e-9, 1e2, 1 - 1e-9, -1e-9, -1e-9],
#     [20, 1.9027e-2 + 1e-12, 0.86229 + 1e-9, np.inf, 1 + 1e-9, 1e-9,
#      1e-9])  # f0, gamma_f, P0, A, amp_cor_re, amp_cor_im, P0_im
bounds = (
    [1, 1.693e-4 - 1e-12, 0.1, 1e2, 1 - 1e-9, -1e-9, -1e-9],
    [20, 1.693e-2 + 1e-12, 1, np.inf, 1 + 1e-9, 1e-9, 1e-9])  # f0, gamma_f, P0, A, amp_cor_re, amp_cor_im, P0_im
########################################################################################################################
NotUsePowerSpectroscopyCalibrate = BackgroundFile.startswith('one_tone') or BackgroundFile.startswith('calibration')
if OneToneFile.endswith('.dat'):
    FileStrList = OneToneFile.split('_')
    UseTwoTone = FileStrList[1] == 'two'
elif OneToneFile.endswith('.hdf5'):
    UseTwoTone = OneToneFile.startswith('two')
UseCH2 = False

if UseTwoTone:
    Calibration = False
    UseCH2 = FileStrList[3] == 'CH2 drive'
    if UseCH2:
        ReadoutPowerStrList = FileStrList[11].split('z')
        ReadoutPower = float(ReadoutPowerStrList[1][:-3])
    else:
        ReadoutPower = float(FileStrList[10][-6:-3])

if NotUsePowerSpectroscopyCalibrate:
    if BackgroundFile.endswith('.dat'):
        [BackFreq, BackComplex] = edf.readFSweepDat(DataPath + BackgroundFile)
        BackPowerStr = BackgroundFile.split('_')[5][:-3]
        BackPower = float(BackPowerStr)
    elif BackgroundFile.endswith('.hdf5'):
        [BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
        BackPower = edf.readReadoutPowerLabber(DataPath + BackgroundFile)
        BackPowerStr = str(BackPower)
    UseOnePowerCalibrate = True
else:
    if BackgroundFile.endswith('.dat'):
        [BackFreq, BackPower, BackComplex] = edf.readFPSweepDat(DataPath + BackgroundFile)
    elif BackgroundFile.endswith('.hdf5'):
        [BackFreq, BackPower, BackComplex] = edf.readFPSweepLabber(DataPath + BackgroundFile)
    if len(np.unique(BackPower)) == 1:
        UseOnePowerCalibrate = True
        BackPower = BackPower[0, 0]
        BackFreq = BackFreq[:, 0]
        BackComplex = BackComplex[:, 0]
        BackPowerStr = str(BackPower)
    else:
        UseOnePowerCalibrate = False

if OneToneFile.endswith('.dat'):
    [OneFreq, OnePower, OneComplex] = edf.readFPSweepDat(DataPath + OneToneFile)
elif OneToneFile.endswith('.hdf5'):
    [OneFreq, OnePower, OneComplex] = edf.readFPSweepLabber(DataPath + OneToneFile)
# OneFreq -= 50e-3
OnePowerUniq = np.unique(OnePower)
NumPower = np.size(OnePowerUniq)
OneFreqUniq = np.unique(OneFreq)
NumFreq = np.size(OneFreqUniq)
if UseOneToneRange:
    StartFreq = OneFreqUniq.min()
    EndFreq = OneFreqUniq.max()
    StartPower = OnePowerUniq.min()
    EndPower = OnePowerUniq.max()
if Calibration:
    RComplex = sbf.FPSweepBackgroundCalibrate(OneFreq, OnePower, OneComplex, BackFreq, BackComplex, BackPower)

OneComplexNormalized = OneComplex * 10 ** (- OnePower / 20)
if UseTwoTone:
    OneComplexNormalized = OneComplex * 10 ** (- ReadoutPower / 20)

BackComplexNormalized = BackComplex * 10 ** (- BackPower / 20)

FreqInd = (EndFreq >= OneFreqUniq) == (OneFreqUniq >= StartFreq)
if len(SelectPower) > 0:
    for i, p in enumerate(SelectPower):
        if i == 0:
            PowerInd = OnePowerUniq == p
        else:
            PowerInd = PowerInd + (OnePowerUniq == p)
else:
    PowerInd = (EndPower >= OnePowerUniq) == (OnePowerUniq >= StartPower)
OneFreqUniqTrunc = OneFreqUniq[FreqInd]
OnePowerUniqTrunc = OnePowerUniq[PowerInd]
NumFreqTrunc = len(OneFreqUniqTrunc)
NumPowerTrunc = len(OnePowerUniqTrunc)
OneComplexNormalizedTrunc = OneComplexNormalized[FreqInd, :]
OneComplexNormalizedTrunc = OneComplexNormalizedTrunc[:, PowerInd]

if Calibration:
    BackFreqUniq = np.unique(BackFreq)
    BackFreqInd = (EndFreq > BackFreqUniq) == (BackFreqUniq > StartFreq)
    BackFreqUniqTrunc = BackFreqUniq[BackFreqInd]
    if UseOnePowerCalibrate:
        BackComplexNormalizedTrunc = BackComplexNormalized[BackFreqInd]
    else:
        BackComplexNormalizedTrunc = BackComplexNormalized[BackFreqInd, :]
    RComplexTrunc = RComplex[FreqInd, :]
    RComplexTrunc = RComplexTrunc[:, PowerInd]

    if FitSeparately:
        optArray = np.zeros((7, NumPowerTrunc))
        covArray = np.zeros((7, NumPowerTrunc))
        print(OnePowerUniqTrunc)
        LargerFreqRangeList = []
        FittedComplexList = []
        UnfittedComplexList = []
        for i, p in enumerate(OnePowerUniqTrunc):
            f0_ind = np.real(RComplexTrunc[:, i]).argmin()
            f0_guess = OneFreqUniqTrunc[f0_ind]
            amp_cor_re_guess = np.real(RComplexTrunc[0, i] + RComplexTrunc[-1, i]) / 2
            amp_cor_im_guess = np.imag(RComplexTrunc[0, i] + RComplexTrunc[-1, i]) / 2

            amp_cor_re_guess = 1
            amp_cor_im_guess = 0

            guess = [f0_guess, gamma_f_guess, P0_guess, A_guess, amp_cor_re_guess, amp_cor_im_guess, P0_im_guess]
            bounds[0][0] = OneFreqUniqTrunc[0]
            bounds[1][0] = OneFreqUniqTrunc[-1]
            if i == 0:
                print(p)
                print(guess)
                print(bounds)
                print(list(OneFreqUniqTrunc))
                print(list(RComplexTrunc[:, i]))
            try:
                opt, cov, LargerFreqRange, FittedComplex = qsf.fitReflectionCircles(OneFreqUniqTrunc, np.array([p]),
                                                                                    RComplexTrunc[:, i],
                                                                                    guess, bounds)
            except ValueError:
                print("Error - curve_fit failed")
                opt = guess
                cov = np.zeros([len(opt), len(opt)])
                LargerFreqRange = OneFreqUniqTrunc
                FittedComplex = RComplexTrunc[:, i] * 0

            if ShiftCircle:
                amp_cor_re_opt = opt[4]
                amp_cor_im_opt = opt[5]
                RComplexTrunc[:, i] = RComplexTrunc[:, i] / (amp_cor_re_opt + amp_cor_im_opt * 1j)
                FittedComplex = FittedComplex / (amp_cor_re_opt + amp_cor_im_opt * 1j)
                amp_cor_im_guess = 0
                amp_cor_re_guess = 1
            UnfittedComplex = qsf.reflection_for_fit([LargerFreqRange, np.array([1e-3 * 10 ** (p / 10)])],
                                                     f0_guess, gamma_f_guess, P0_guess, A_guess,
                                                     amp_cor_re_guess, amp_cor_im_guess, P0_im_guess)
            split_ind = int(len(UnfittedComplex) / 2)
            UnfittedComplex = UnfittedComplex[:split_ind] + 1j * UnfittedComplex[split_ind:]
            optArray[:, i] = np.array(opt)
            LargerFreqRangeList.append(LargerFreqRange)
            FittedComplexList.append(FittedComplex)
            UnfittedComplexList.append(UnfittedComplex)
        f0_fit, gamma_f_fit, P0_fit, A_fit, amp_cor_re_fit, amp_cor_im_fit, P0_im_fit = optArray.mean(axis=1)
    else:
        f0_ind = np.real(RComplexTrunc[:, 0]).argmin()
        f0_guess = OneFreqUniqTrunc[f0_ind]
        amp_cor_re_guess = np.real(RComplexTrunc[0, 0] + RComplexTrunc[-1, 0]) / 2
        amp_cor_im_guess = np.imag(RComplexTrunc[0, 0] + RComplexTrunc[-1, 0]) / 2

        amp_cor_re_guess = 1
        amp_cor_im_guess = 0

        guess = [f0_guess, gamma_f_guess, P0_guess, A_guess, amp_cor_re_guess, amp_cor_im_guess, P0_im_guess]
        try:
            opt, cov, LargerFreqRange, FittedComplex = qsf.fitReflectionCircles(OneFreqUniqTrunc, OnePowerUniqTrunc,
                                                                                RComplexTrunc,
                                                                                guess, bounds)
            f0_fit, gamma_f_fit, P0_fit, A_fit, amp_cor_re_fit, amp_cor_im_fit, P0_im_fit = opt
        except ValueError:
            f0_fit, gamma_f_fit, P0_fit, A_fit, amp_cor_re_fit, amp_cor_im_fit, P0_im_fit = guess
            LargerFreqRange = OneFreqUniqTrunc
            FittedComplex = RComplexTrunc
        if ShiftCircle:
            RComplexTrunc[:, :] = RComplexTrunc[:, :] / (amp_cor_re_fit + amp_cor_im_fit * 1j)
            FittedComplex = FittedComplex / (amp_cor_re_fit + amp_cor_im_fit * 1j)
            amp_cor_im_fit = 0
            amp_cor_re_fit = 1
        if RotateCircle:
            RComplexTrunc[:, :] = (RComplexTrunc[:, :] - 1) / (P0_fit + P0_im_fit * 1j) * np.abs(
                P0_fit + P0_im_fit * 1j) + 1
            FittedComplex = (FittedComplex - 1) / (P0_fit + P0_im_fit * 1j) * np.abs(P0_fit + P0_im_fit * 1j) + 1
            P0_fit = np.abs(P0_fit + P0_im_fit * 1j)
            P0_im_fit = 0

print('(1+%.5G*1e-3*10**(4.13/10))/2=%.5G' % (A_fit, (1 + A_fit * 1e-3 * 10 ** (4.13 / 10)) / 2))
########################################################################################################################

if Calibration:
    fig, ax = plt.subplots()
    leg = ()
    for i in range(NumPowerTrunc):
        leg += (str(OnePowerUniqTrunc[i]) + 'dBm',)
        # plt.scatter(np.real(RComplexTrunc[:, i]), np.imag(RComplexTrunc[:, i]), s=10)
        plt.plot(OneFreqUniqTrunc, np.real(RComplexTrunc[:, i]), LineSpec)
    for i in range(NumPowerTrunc):
        if FitSeparately:
            plt.plot(LargerFreqRangeList[i][:-1], np.real(FittedComplexList[i][:-1]), 'r')
            if PlotUnfittedCircle:
                plt.plot(LargerFreqRangeList[i][:-1], np.real(UnfittedComplexList[i][:-1]), 'b')
        else:
            plt.plot(LargerFreqRange[:-1], np.real(FittedComplex[:-1, i]), 'r')
    plt.xlim(StartFreq, EndFreq)
    plt.xlabel('freq(GHz)', fontsize='x-large')
    plt.ylabel('Re(r)', fontsize='x-large')
    plt.legend(leg)
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    locs, labels = plt.xticks()
    print(locs)
    print(labels)

if Calibration:
    fig, ax = plt.subplots()
    leg = ()
    for i in range(NumPowerTrunc):
        # plt.scatter(np.real(RComplexTrunc[:, i]), np.imag(RComplexTrunc[:, i]), s=10)
        leg += (str(OnePowerUniqTrunc[i]) + 'dBm',)
        plt.plot(OneFreqUniqTrunc, np.imag(RComplexTrunc[:, i]), LineSpec)
    for i in range(NumPowerTrunc):
        if FitSeparately:
            plt.plot(LargerFreqRangeList[i][:-1], np.imag(FittedComplexList[i][:-1]), 'r')
            if PlotUnfittedCircle:
                plt.plot(LargerFreqRangeList[i][:-1], np.imag(UnfittedComplexList[i][:-1]), 'b')
        else:
            plt.plot(LargerFreqRange[:-1], np.imag(FittedComplex[:-1, i]), 'r')
    plt.xlim(StartFreq, EndFreq)
    plt.xlabel('freq(GHz)', fontsize='x-large')
    plt.ylabel('Im(r)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    plt.legend(leg)

fig, ax = plt.subplots()
leg = ()
for i in range(NumPowerTrunc):
    plt.plot(OneFreqUniqTrunc, np.imag(OneComplexNormalizedTrunc[:, i]), LineSpec)
    # if not UseOnePowerCalibrate:
    #     plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), '--')
    if UseTwoTone:
        leg += ('Two tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
    else:
        leg += ('One tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
if UseOnePowerCalibrate and Calibration:
    plt.plot(BackFreqUniqTrunc, np.imag(BackComplexNormalizedTrunc), '--')
    leg += ('Background ' + BackPowerStr + 'dBm',)
plt.xlabel('freq(GHz)', fontsize='x-large')
plt.ylabel('imag', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.legend(leg)

if not UseOnePowerCalibrate:
    fig, ax = plt.subplots()
    leg = ()
    for i in range(NumPowerTrunc):
        plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), '--')
        leg += ('Background ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
    plt.xlabel('freq(GHz)', fontsize='x-large')
    plt.ylabel('Abs', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    plt.legend(leg)

fig, ax = plt.subplots()
leg = ()
for i in range(NumPowerTrunc):
    phase = np.unwrap(np.angle(OneComplexNormalizedTrunc[:, i]))
    PhaseSlope = (phase[-1] - phase[0]) / (OneFreqUniqTrunc[-1] - OneFreqUniqTrunc[0])
    phase = phase - PhaseSlope * (OneFreqUniqTrunc - OneFreqUniqTrunc[0]) + np.pi
    phase = phase % (2 * np.pi)
    plt.plot(OneFreqUniqTrunc, phase, LineSpec)
    # plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), '--')
    if UseTwoTone:
        leg += ('Two tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
    else:
        leg += ('One tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
if UseOnePowerCalibrate and Calibration:
    phase = np.unwrap(np.angle(BackComplexNormalizedTrunc))
    PhaseSlope = (phase[-1] - phase[0]) / (BackFreqUniqTrunc[-1] - BackFreqUniqTrunc[0])
    phase = phase - PhaseSlope * (BackFreqUniqTrunc - BackFreqUniqTrunc[0]) + np.pi
    phase = phase % (2 * np.pi)
    plt.plot(BackFreqUniqTrunc, phase, '--')
    leg += ('Background ' + BackPowerStr + 'dBm',)
plt.xlabel('freq(GHz)', fontsize='x-large')
plt.ylabel('Phase', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.legend(leg)

#
# fig, ax = plt.subplots()
# leg = ()
# for i in range(NumPowerTrunc):
#     plt.plot(OneFreqUniqTrunc, np.real(OneComplexNormalizedTrunc[:, i]))
#     # plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), '--')
#     leg += ('One tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
# if UseOnePowerCalibrate and Calibration:
#     plt.plot(BackFreqUniqTrunc, np.real(BackComplexNormalizedTrunc), '--')
#     leg += ('Background ' + BackPowerStr + 'dBm',)
# plt.xlabel('freq(GHz)', fontsize='x-large')
# plt.ylabel('Re(r)', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')
# plt.tight_layout()
# plt.legend(leg)

# if Calibration:
#     fig, ax = plt.subplots()
#     leg = ()
#     BackPhaseNormalizedTrunc = np.unwrap(np.angle(BackComplexNormalizedTrunc))
#     PhaseSlope = (BackPhaseNormalizedTrunc[-1] - BackPhaseNormalizedTrunc[0]) / (
#             BackFreqUniqTrunc[-1] - BackFreqUniqTrunc[0])
#     for i in range(NumPowerTrunc):
#         plt.plot(OneFreqUniqTrunc, np.unwrap(np.angle(OneComplexNormalizedTrunc[:, i])) - PhaseSlope * (
#                 OneFreqUniqTrunc - OneFreqUniqTrunc[0]))
#         # plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), '--')
#         leg += ('One tone ' + str(OnePowerUniqTrunc[i]) + 'dBm',)
#     if UseOnePowerCalibrate:
#         plt.plot(BackFreqUniqTrunc, BackPhaseNormalizedTrunc - PhaseSlope * (BackFreqUniqTrunc - BackFreqUniqTrunc[0]),
#                  '--')
#         leg += ('Background ' + BackPowerStr + 'dBm',)
#     plt.xlabel('freq(GHz)', fontsize='x-large')
#     plt.ylabel('Phase', fontsize='x-large')
#     plt.tick_params(axis='both', which='major', labelsize='x-large')
#     plt.tight_layout()
#     plt.legend(leg)

if Calibration:
    if not UseOnePowerCalibrate:
        fig, ax = plt.subplots()
        for i in range(NumPowerTrunc):
            plt.plot(BackFreqUniqTrunc, np.abs(BackComplexNormalizedTrunc[:, i]), LineSpec)
        plt.xlabel('freq(GHz)', fontsize='x-large')
        plt.ylabel('Abs', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

    fig, ax = plt.subplots()
    limit = 1.7
    leg = ()
    for i in range(NumPowerTrunc):
        leg += (str(OnePowerUniqTrunc[i]) + 'dBm',)
        plt.plot(np.real(RComplexTrunc[:, i]), np.imag(RComplexTrunc[:, i]), LineSpec)
    for i in range(NumPowerTrunc):
        if FitSeparately:
            plt.plot(np.real(FittedComplexList[i]), np.imag(FittedComplexList[i]), 'r')
            if PlotUnfittedCircle:
                plt.plot(np.real(UnfittedComplexList[i]), np.imag(UnfittedComplexList[i]), 'b')
                if i == 0:
                    print(list(UnfittedComplexList[i]))
        else:
            plt.plot(np.real(FittedComplex[:, i]), np.imag(FittedComplex[:, i]), 'r')
    plt.plot([-2, 2], [0, 0], '--')
    plt.plot([1], [0], 'ro')
    plt.xlim(-limit, limit)
    plt.ylim(-limit, limit)
    plt.title('f0=%.6GGHz, gamma_f=%.5GMHz, P0=%.5G, A=%.5G,\n amp_cor_re=%0.2G, amp_cor_im=%0.2G, P0_im=%0.2G' % (
        f0_fit, gamma_f_fit * 1e3, P0_fit, A_fit, amp_cor_re_fit, amp_cor_im_fit, P0_im_fit))
    plt.xlabel('Re(r)', fontsize='x-large')
    plt.ylabel('Im(r)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    ax.set_aspect('equal')
    plt.legend(leg)
    print('Optimal power =', 10 * np.log10(1 / A_fit * 1000), 'dBm')
param_name = ['f0', 'gamma_f', 'P0', 'A', 'amp_cor_re', 'amp_cor_im', 'P0_im']
if FitSeparately and PlotParamVSPower and Calibration:
    for i, p in enumerate(param_name):
        fig, ax = plt.subplots()
        plt.plot(OnePowerUniqTrunc, optArray[i, :], 'o')
        plt.xlabel('Power/dBm', fontsize='x-large')
        plt.ylabel(p, fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

plt.show()
