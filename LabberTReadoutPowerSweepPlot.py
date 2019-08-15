import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, rabi_curve, FitTransientTime, AutoRotate
import ExtractDataFunc as edf
import h5py

DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium032619/'
BackgroundFile = 'calibration_5.hdf5'


LabberFileList = [
    't1_2019-08-14-10.hdf5',
]

IQModFreq = 0.05
FitCorrectedR = True
LogScale = False
Calibration = False
RotateComplex = True

PhaseSlope = 326.7041108065019
PhaseRefrenceFreq = 4.105

NumFile = len(LabberFileList)
# analyze background file

[BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
BackPower = edf.readQubitPowerLabber(DataPath + BackgroundFile)
BackPowerStr = str(BackPower)

for i, RabiFile in enumerate(LabberFileList):
    RabiFileStrList = RabiFile[:-5].split('_')
    MeasurementType = RabiFileStrList[0]

    ReadoutFreq = edf.readQubitFreqLabber(DataPath + RabiFile)
    QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)
    # DrivePower = edf.readPumpPowerLabber(DataPath + RabiFile)

    if MeasurementType in ('rabi', 'transient'):
        [Time, DrivePower, ComplexRabi] = edf.readRabiPowerSweepLabber(DataPath + RabiFile)
    elif MeasurementType == 't1':
        [Time, ReadoutPower, ComplexRabi] = edf.readT1ReadoutPowerSweepLabber(DataPath + RabiFile)
    TimeFit = np.linspace(Time.min(), Time.max(), 200)
    ComplexRabiNormalized = ComplexRabi * 10 ** (- ReadoutPower / 20)
    if Calibration:
        RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexRabi, BackFreq,
                                                  BackComplex, BackPower)
    else:
        RComplex = ComplexRabiNormalized

    for j, power in enumerate(ReadoutPower):
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
            opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
            A_fit, T1_fit, B_fit = opt
            FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
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
            FitR = rabi_curve(TimeFit, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
            ParamList = ['A', 'Decay time/ns', 'B', 'Tpi', 'phi0']
        if i == 0 and j == 0:
            TimeList = [Time]
            TimeFitList = [TimeFit]
            RComplexList = [RComplex]
            FitRList = [FitR]
            ReadoutPowerArray = ReadoutPower
            OptMatrix = np.reshape(opt, (len(opt), 1))
            ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
        else:
            OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
            ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)
            ReadoutPowerArray = np.concatenate((ReadoutPowerArray, ReadoutPowerArray))
            TimeList.append(Time)
            RComplexList.append(RComplex)
            FitRList.append(FitR)
        for i_p, par in enumerate(opt):
            if np.abs(ErrMatrix[i_p, -1]) > np.abs(par) or opt[1] > 5 * Time[-1]:
                OptMatrix[i_p, -1] = np.nan
            else:
                OptMatrix[i_p, -1] = par


limit = 1.7

fig, ax = plt.subplots()
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
for i in range(NumFile):
    plt.plot(TimeList[i], np.real(RComplexList[i]), 'o')
    plt.plot(TimeFitList[i], FitRList[i])
plt.xlabel('Time(ns)', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
plotInd = 1
# plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
ax.errorbar(ReadoutPowerArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
plt.xlabel('Power(dBm)', fontsize='x-large')
plt.ylabel('Decay time(us)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if LogScale:
    ax.set_yscale('log')

if MeasurementType == 'rabi':
    fig, ax = plt.subplots()
    plotInd = 2
    # plt.plot(ReadoutPowerArray, OptMatrix[plotInd, :]/1000, 'o')
    ax.errorbar(ReadoutPowerArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (dBm)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    # plt.ylim(0, 20000)

    fig, ax = plt.subplots()
    plotInd = 3
    # plt.plot(ReadoutPowerArray, OptMatrix[plotInd, :]/1000, 'o')
    ax.errorbar(ReadoutPowerArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (dBm)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    ax.set_yscale('log')

    fig, ax = plt.subplots()
    plotInd = 3
    ReadoutPowerWattArray = 10 ** (ReadoutPowerArray / 10) / 1000
    ax.errorbar(ReadoutPowerWattArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Power (W)', fontsize='x-large')
    plt.ylabel('T_pi (ns)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    ax.set_yscale('log')

plt.show()
