import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, rabi_curve, FitTransientTime, AutoRotate
import ExtractDataFunc as edf
import h5py

DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
BackgroundFile = 'power spectroscopy_116.hdf5'
BackgroundFile = 'power spectroscopy_76.hdf5'
RabiFileList = [
# 'transient_10.hdf5',
'transient_9.hdf5',
]

IQModFreq = 0.05
# FitForGamma = True
Gamma_r = 1.807 * np.pi * 2
FitCorrectedR = False
LogScale = False
Calibration = True
RotateComplex = False

PhaseSlope = 326.7041108065019
PhaseRefrenceFreq = 4.105

NumFile = len(RabiFileList)
DrivePowerArray = np.zeros([NumFile, ])
# analyze background file

[BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
BackPower = edf.readReadoutPowerLabber(DataPath + BackgroundFile)
BackPowerStr = str(BackPower)

for i, RabiFile in enumerate(RabiFileList):
    RabiFileStrList = RabiFile[:-5].split('_')
    MeasurementType = RabiFileStrList[0]

    ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
    ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
    QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)
    # DrivePower = edf.readPumpPowerLabber(DataPath + RabiFile)

    if MeasurementType in ('rabi', 'transient'):
        [Time, DrivePower, ComplexRabi] = edf.readRabiPowerSweepLabber(DataPath + RabiFile)
    elif MeasurementType == 't1':
        [Time, DrivePower, ComplexRabi] = edf.readT1PowerSweepLabber(DataPath + RabiFile)

    if RabiFile == 'transient_28.hdf5':
        ind = DrivePower > -45
        DrivePower = DrivePower[ind]
        ComplexRabi = ComplexRabi[:, ind]

    ComplexRabiNormalized = ComplexRabi * 10 ** (- ReadoutPower / 20)
    if Calibration:
        RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexRabi, BackFreq,
                                                  BackComplex, BackPower)
    else:
        RComplex = ComplexRabiNormalized

    for j, power in enumerate(DrivePower):
        if RotateComplex:
            RComplex[:, j] = AutoRotate(RComplex[:, j])
        y_data = np.array(RComplex[:, j].real, dtype='float64')
        x_data = np.array(Time, dtype='float64')
        if MeasurementType in ('t1', 'transient'):
            B_guess = y_data[-1]
            A_guess = y_data[0].real - B_guess
            T1_guess = x_data[-1] / 2
            bounds = (
                (-2, 1, -1),
                (2, 1e6, 1)
            )
            opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], bounds=bounds, maxfev=30000)
            A_fit, T1_fit, B_fit = opt
            FitR = T1_curve(Time, A_fit, T1_fit, B_fit)
            ParamList = ['A', 'Decay time/ns', 'B']
        elif MeasurementType in ('rabi'):
            B_guess = y_data.mean()
            A_guess = y_data[0] - B_guess
            T1_guess = x_data[-1]
            MaxInd = y_data.argmax()
            MinInd = y_data.argmin()
            Tpi_guess = np.abs(x_data[MaxInd] - x_data[MinInd])
            phi0_guess = 0
            guess = ([A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess])
            bounds = (
                (-2, 1, -1, 1, - np.pi / 2),
                (2, np.inf, 1, np.inf, np.pi / 2)
            )
            opt, cov = curve_fit(rabi_curve, x_data, y_data, p0=guess, bounds=bounds)
            A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit = opt
            FitR = rabi_curve(Time, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
            ParamList = ['A', 'Decay time/ns', 'B', 'Tpi', 'phi0']
        if i == 0 and j == 0:
            TimeList = []
            RComplexList = []
            FitRList = []
            DrivePowerArray = np.array([power])
            OptMatrix = np.reshape(opt, (len(opt), 1))
            ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
        else:
            OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
            ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)
            DrivePowerArray = np.concatenate((DrivePowerArray, np.array([power])))
        TimeList.append(Time)
        RComplexList.append(RComplex)
        # print(FitR)
        FitRList.append(FitR)
        for i_p, par in enumerate(opt):
            if np.abs(ErrMatrix[i_p, -1]) > np.abs(par) or opt[1] > 5 * Time[-1]:
                OptMatrix[i_p, -1] = np.nan
                # OptMatrix[i_p, -1] = par
            else:
                OptMatrix[i_p, -1] = par

if MeasurementType == 'transient':
    print(OptMatrix)
    power_for_plot = np.linspace(np.min(DrivePowerArray) - 20, np.max(DrivePowerArray) + 10, 101)
    opt, cov, FitTime = FitTransientTime(DrivePowerArray, Gamma_r, OptMatrix, power_for_plot=power_for_plot)
    Gamma_in_fit, Gamma_out_fit, pow_ratio_fit = opt
    TransientOpt = opt
    TransientErr = np.sqrt(cov.diagonal())
    TransientParamList = ['Gamma_in', 'Gamma_out', 'pow_ratio']
    print(opt)

limit = 1.7

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for i in range(NumFile):
    plt.plot(np.real(RComplexList[i]), np.imag(RComplexList[i]))
if Calibration:
    plt.plot([-2, 2], [0, 0], '--')
    plt.plot([1], [0], 'ro')
plt.xlabel('Re', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if Calibration:
    plt.xlim(-limit, limit)
    plt.ylim(-limit, limit)
ax.set_aspect('equal')

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for i in range(len(TimeList)):
    plt.plot(TimeList[i], np.real(RComplexList[i]), 'o')
    plt.plot(TimeList[i], FitRList[i])
plt.xlabel('Time/ns', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
ax.grid(linestyle='--')
plotInd = 1
# plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
ax.errorbar(DrivePowerArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
if MeasurementType == 'transient':
    plt.plot(power_for_plot, FitTime)
plt.xlabel('Power/dBm', fontsize='x-large')
plt.ylabel('Decay time/us', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
if MeasurementType == 'transient':
    plt.title('Gamma_in=%.3G$\pm$%.3Gus$^{-1}$, Gamma_out=%.3G$\pm$%.3Gus$^{-1}$,\n pow_ratio=%.3G$\pm$%.3G' % (
        Gamma_in_fit, TransientErr[0], Gamma_out_fit, TransientErr[1], pow_ratio_fit, TransientErr[2]))
plt.tight_layout()
if LogScale:
    ax.set_yscale('log')

fig, ax = plt.subplots()
ax.grid(linestyle='--')
plotInd = 1
DrivePowerArray_MHz = np.sqrt(1e-3 * 10 ** (DrivePowerArray / 10) * pow_ratio_fit) / 2 / np.pi
power_for_plot_MHz = np.sqrt(1e-3 * 10 ** (power_for_plot / 10) * pow_ratio_fit) / 2 / np.pi
branching_ratio = Gamma_r / Gamma_out_fit
ax.errorbar(DrivePowerArray_MHz, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
if MeasurementType == 'transient':
    plt.plot(power_for_plot_MHz, FitTime)
plt.xlabel('$\Omega/2\pi$(MHz)', fontsize='x-large')
plt.ylabel('Decay tim(us)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
if MeasurementType == 'transient':
    plt.title('Gamma_in=(%.3G$\pm$%.3Gus)$^{-1}$, Gamma_out=(%.3G$\pm$%.3Gus)$^{-1}$,\n branching_ratio=%.3G' % (
        1 / Gamma_in_fit, TransientErr[0] / Gamma_in_fit ** 2, 1 / Gamma_out_fit, TransientErr[1] / Gamma_out_fit ** 2,
        branching_ratio))
plt.tight_layout()
ax.set_xscale('log')
if LogScale:
    ax.set_yscale('log')

if MeasurementType == 'rabi':
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plotInd = 2
    # plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
    ax.errorbar(DrivePowerArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (dBm)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    # plt.ylim(0, 20000)

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plotInd = 3
    # plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
    ax.errorbar(DrivePowerArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (dBm)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    ax.set_yscale('log')

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plotInd = 3
    DrivePowerWattArray = 10 ** (DrivePowerArray / 10) / 1000
    ax.errorbar(DrivePowerWattArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (W)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    ax.set_yscale('log')

plt.show()
