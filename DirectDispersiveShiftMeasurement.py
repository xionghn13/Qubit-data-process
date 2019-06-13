from scipy import *
import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
import ExtractDataFunc as edf


def lorenztian(f, f0, kappa, A, B):
    t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B
    return t


# DataPath = 'E:/Projects\Fluxonium\data_process/cavity/8.5GHz/'
DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium032619/'
OneToneFile = 'rabi_175.hdf5'

TruncateFreq = False

StartFreq = 7.78
EndFreq = 7.8

[Time, Freq, Complex] = edf.readRabiFreqSweepLabber(DataPath + OneToneFile)
# [Freq, Complex] = edf.readVNAS21(DataPath + OneToneFile)
# Complex = np.sqrt(Complex)
Complex = Complex ** 2
NumTime = len(Time)

if TruncateFreq:
    FreqInd = (EndFreq >= Freq) == (Freq >= StartFreq)
    FreqTrunc = Freq[FreqInd]
    ComplexTrunc = Complex[:, FreqInd]
else:
    FreqTrunc = Freq
    ComplexTrunc = Complex

AbsComplex = np.abs(ComplexTrunc)
MaxAbs = np.max(AbsComplex)
ComplexTrunc /= MaxAbs
AbsComplex = np.abs(ComplexTrunc)
for i in range(NumTime):
    MaxAbs = np.max(AbsComplex[i, :])
    MinAbs = np.min(AbsComplex[i, :])
    MaxInd = AbsComplex[i, :].argmax()
    f0_guess = FreqTrunc[MaxInd]
    kappa_guess = (FreqTrunc[-1] - FreqTrunc[0]) / 4
    B_guess = MinAbs
    A_guess = (MaxAbs - MinAbs) * (kappa_guess / 2) ** 2

    guess = ([f0_guess, kappa_guess, A_guess, B_guess])
    bounds = (
        (Freq[0], 0, 0, 0),
        (Freq[-1], kappa_guess * 4, MaxAbs * 10, MaxAbs)
    )

    qopt, qcov = curve_fit(lorenztian, FreqTrunc, AbsComplex[i, :], guess, bounds=bounds)
    f0_fit, kappa_fit, A_fit, B_fit = qopt
    # kappa_fit = np.abs(kappa_fit)
    # AbsGuess = lorenztian(FreqTrunc, f0_guess, kappa_guess, A_guess, B_guess)
    AbsFit = lorenztian(FreqTrunc, f0_fit, kappa_fit, A_fit, B_fit)
    if i == 0:
        OptMatrix = np.reshape(qopt, (len(qopt), 1))
        ErrMatrix = np.reshape(np.sqrt(qcov.diagonal()), (len(qopt), 1))
        AbsFitList = []
    else:
        OptMatrix = np.concatenate((OptMatrix, np.reshape(qopt, (len(qopt), 1))), axis=1)
        ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(qcov.diagonal()), (len(qopt), 1))), axis=1)
        AbsFitList.append(AbsFit)

f0Array = OptMatrix[0, :]
MaxInd = f0Array.argmax()
MinInd = f0Array.argmin()
Maxf0 = np.max(f0Array)
Minf0 = np.min(f0Array)
print('Min f0=%.5GGHz, max f0=%.5GGHz, difference=%.5GMHz' % (Minf0, Maxf0, (Maxf0 - Minf0) * 1e3))

fig, ax = plt.subplots()
leg = []
for i in [MinInd, MaxInd]:
    plt.plot(FreqTrunc, AbsComplex[i, :], '.')
    plt.plot(FreqTrunc, AbsFitList[i], 'r')
    leg += ['driving time=%.3Gns' % Time[i]]
    # plt.plot(FreqTrunc, AbsGuess, 'y')
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Abs', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
plotInd = 0
leg = ()
ax.errorbar(Time, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
plt.xlabel('Time(ns)', fontsize='x-large')
plt.ylabel('f0(GHz)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
