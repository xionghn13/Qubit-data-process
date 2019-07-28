import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, rabi_curve, FitTransientTime, AutoRotate
import ExtractDataFunc as edf
import h5py

FixedFolder = 'E:/Projects\Fluxonium\data_process/Fluxonium032619/'
FixedFolder = ''
T1FileList = [
    't1_2019-07-23-18-15-23.hdf5',
    't1_2019-07-23-16-56-18.hdf5',
    't1_2019-07-23-15-40-42.hdf5',
    't1_2019-07-23-11-25-21.hdf5',
    't1_2019-07-23-12-23-53.hdf5',
    't1_2019-07-23-13-26-36.hdf5',
    't1_2019-07-23-14-29-44.hdf5',
    't1_2019-07-23-10-30-49.hdf5',
    't1_2019-07-23-09-40-32.hdf5',
    't1_2019-07-23-08-53-59.hdf5',
]

LogScale = False
RotateComplex = True

NumFile = len(T1FileList)
PulseLengthArray = np.zeros([NumFile, ])


for i, T1File in enumerate(T1FileList):
    if FixedFolder == '':
        DataPath = edf.getFolder(T1File)
    else:
        DataPath = FixedFolder
    T1FileStrList = T1File[:-5].split('_')
    MeasurementType = T1FileStrList[0]

    ReadoutFreq = edf.readQubitFreqLabber(DataPath + T1File)
    ReadoutPower = edf.readQubitPowerLabber(DataPath + T1File)
    QubitFreq = edf.readPumpFreqLabber(DataPath + T1File)
    # DrivePower = edf.readPumpPowerLabber(DataPath + T1File)
    DrivingPulseLen = edf.readDrivingPulseLenLabber(DataPath + T1File) * 1e6

    [Time, Complex] = edf.readT1Labber(DataPath + T1File)

    if RotateComplex:
        Complex = AutoRotate(Complex)
    y_data = np.array(Complex.real, dtype='float64')
    x_data = np.array(Time, dtype='float64')
    B_guess = y_data[-1]
    A_guess = y_data[0].real - B_guess
    T1_guess = x_data[-1] / 2
    bounds = (
        (-2, 1, -1),
        (2, 1e6, 1)
    )
    opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
    A_fit, T1_fit, B_fit = opt
    FitY = T1_curve(Time, A_fit, T1_fit, B_fit)
    ParamList = ['A', 'Decay time/ns', 'B']

    if i == 0:
        TimeList = []
        ComplexList = []
        FitYList = []
        PulseLengthArray = np.array([DrivingPulseLen])
        OptMatrix = np.reshape(opt, (len(opt), 1))
        ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
    else:
        OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
        ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)
        PulseLengthArray = np.concatenate((PulseLengthArray, np.array([DrivingPulseLen])))
    TimeList.append(Time)
    ComplexList.append(Complex)
    FitYList.append(FitY)
    for i_p, par in enumerate(opt):
        if np.abs(ErrMatrix[i_p, -1]) > np.abs(par) or opt[1] > 5 * Time[-1]:
            OptMatrix[i_p, -1] = np.nan
        else:
            OptMatrix[i_p, -1] = par

limit = 1.7

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(np.real(ComplexList[i]), np.imag(ComplexList[i]))
plt.xlabel('Re', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_aspect('equal')

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(TimeList[i], np.real(ComplexList[i]), 'o')
    plt.plot(TimeList[i], FitYList[i])
plt.xlabel('Time/ns', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
plotInd = 1
# plt.plot(PulseLengthArray, OptMatrix[plotInd, :]/1000, 'o')
ax.errorbar(PulseLengthArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
plt.xlabel('Driving Pulse Length/us', fontsize='x-large')
plt.ylabel('Decay time/us', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if LogScale:
    ax.set_yscale('log')
plt.show()
