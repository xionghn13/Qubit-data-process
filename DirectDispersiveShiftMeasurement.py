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


def two_lorenztian_same_kappa(f, kappa, f0, f1, A0, A1, B):
    t = lorenztian(f, f0, kappa, A0, B) + lorenztian(f, f1, kappa, A1, B)
    return t


def two_lorenztian_same_kappa_two_curves(f, kappa, f0, f1, A00, A01, A10, A11, B):
    t1 = lorenztian(f, f0, kappa, A00, B) + lorenztian(f, f1, kappa, A01, B)
    t2 = lorenztian(f, f0, kappa, A10, B) + lorenztian(f, f1, kappa, A11, B)
    t = np.concatenate((t1, t2), axis=0)
    return t

# DataPath = 'E:/Projects\Fluxonium\data_process/cavity/8.5GHz/'
DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process\Fluxonium032619/'
OneToneFile = 'rabi_3.hdf5'

TruncateFreq = False

StartFreq = 7.96
EndFreq = 7.97

fit_single_lorenztian = False
fit_all_together = True

MaxFreqInd = 8
MinFreqInd = 0

[Time, Freq, Complex] = edf.readRabiFreqSweepLabber(DataPath + OneToneFile)
# [Freq, Complex] = edf.readVNAS21(DataPath + OneToneFile)
# Complex = np.sqrt(Complex)
Complex = Complex ** 2
NumTime = len(Time)
# print(Freq)

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
if fit_all_together:
    MaxAbs = np.max(AbsComplex[0, :])
    MinAbs = np.min(AbsComplex[0, :])
    MaxInd = AbsComplex[0, :].argmax()
    f0_guess = FreqTrunc[MaxInd]
    kappa_guess = (FreqTrunc[-1] - FreqTrunc[0]) / 4
    B_guess = MinAbs
    A_guess = (MaxAbs - MinAbs) * (kappa_guess / 2) ** 2
    f1_guess = (FreqTrunc[-1] + FreqTrunc[0]) - f0_guess

    guess = ([kappa_guess, f0_guess, f1_guess, A_guess, A_guess, A_guess, A_guess, B_guess])

    bounds = (
        (0, Freq[0], Freq[0], 0, 0, 0, 0, 0),
        (kappa_guess * 4, Freq[-1], Freq[-1], MaxAbs * 10, MaxAbs * 10, MaxAbs * 10, MaxAbs * 10, MaxAbs)
    )
    AbsComplexToFit = AbsComplex[[MinFreqInd, MaxFreqInd], :].ravel()
    qopt, qcov = curve_fit(two_lorenztian_same_kappa_two_curves, FreqTrunc, AbsComplexToFit, guess, bounds=bounds)
    AbsFit = two_lorenztian_same_kappa_two_curves(FreqTrunc, *qopt)

    OptMatrix = np.tensordot(qopt, np.ones(NumTime), axes=0)
    qerr = np.sqrt(qcov.diagonal())
    ErrMatrix = np.tensordot(qerr, np.ones(NumTime), axes=0)
    AbsFitList = list(AbsFit.reshape((2, len(FreqTrunc))))
else:
    for i in range(NumTime):
        MaxAbs = np.max(AbsComplex[i, :])
        MinAbs = np.min(AbsComplex[i, :])
        MaxInd = AbsComplex[i, :].argmax()
        f0_guess = FreqTrunc[MaxInd]
        kappa_guess = (FreqTrunc[-1] - FreqTrunc[0]) / 4
        B_guess = MinAbs
        A_guess = (MaxAbs - MinAbs) * (kappa_guess / 2) ** 2

        if fit_single_lorenztian:
            guess = ([f0_guess, kappa_guess, A_guess, B_guess])
            bounds = (
                (Freq[0], 0, 0, 0),
                (Freq[-1], kappa_guess * 4, MaxAbs * 10, MaxAbs)
            )
            qopt, qcov = curve_fit(lorenztian, FreqTrunc, AbsComplex[i, :], guess, bounds=bounds)
            f0_fit, kappa_fit, A_fit, B_fit = qopt
            AbsFit = lorenztian(FreqTrunc, f0_fit, kappa_fit, A_fit, B_fit)
        else:
            f1_guess = (FreqTrunc[-1] + FreqTrunc[0]) - f0_guess
            guess = ([kappa_guess, f0_guess, f1_guess, A_guess, A_guess, B_guess])
            bounds = (
                (0, Freq[0], Freq[0], 0, 0, 0),
                (kappa_guess * 4, Freq[-1], Freq[-1], MaxAbs * 10, MaxAbs * 10, MaxAbs)
            )
            qopt, qcov = curve_fit(two_lorenztian_same_kappa, FreqTrunc, AbsComplex[i, :], guess, bounds=bounds)
            AbsFit = two_lorenztian_same_kappa(FreqTrunc, *qopt)


        if i == 0:
            OptMatrix = np.reshape(qopt, (len(qopt), 1))
            ErrMatrix = np.reshape(np.sqrt(qcov.diagonal()), (len(qopt), 1))
            AbsFitList = [AbsFit]
        else:
            OptMatrix = np.concatenate((OptMatrix, np.reshape(qopt, (len(qopt), 1))), axis=1)
            ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(qcov.diagonal()), (len(qopt), 1))), axis=1)
            AbsFitList.append(AbsFit)

if fit_single_lorenztian:
    f0Array = OptMatrix[0, :]
    kappa_array = OptMatrix[1, :]
else:
    f0Array = OptMatrix[1, :]
    kappa_array = OptMatrix[0, :]

# MaxFreqInd = f0Array.argmax()
# MinFreqInd = f0Array.argmin()

Maxf0 = f0Array[MaxFreqInd]
Minf0 = f0Array[MinFreqInd]
print('Min f0=%.5GGHz, max f0=%.5GGHz, difference=%.5GMHz' % (Minf0, Maxf0, (Maxf0 - Minf0) * 1e3))
print('Min kappa=%.5GMHz, max kappa=%.5GMHz' % (kappa_array[MinFreqInd] * 1e3, kappa_array[MaxFreqInd] * 1e3))
# print(MaxInd)
# print(MinInd)
# print(AbsFitList.__len__())
# print(f0Array.__len__())
if fit_all_together:
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    for i in [MinFreqInd, MaxFreqInd]:
        plt.plot(FreqTrunc, AbsComplex[i, :], '-')
        leg += ['driving time=%.3Gns' % Time[i]]
        # plt.plot(FreqTrunc, AbsGuess, 'y')
    print(leg)
    leg = ['no driving', 'after $\pi$ pulse']
    plt.legend(leg)
    for i in range(2):
        plt.plot(FreqTrunc, AbsFitList[i], '--')
    plt.xlabel('freq/GHz', fontsize='x-large')
    plt.ylabel('Abs', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

else:
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    for i in [MinFreqInd, MaxFreqInd]:
        plt.plot(FreqTrunc, AbsComplex[i, :], '-')
        leg += ['driving time=%.3Gns' % Time[i]]
        # plt.plot(FreqTrunc, AbsGuess, 'y')
    print(leg)
    leg = ['no driving', 'after $\pi$ pulse']
    plt.legend(leg)
    for i in [MinFreqInd, MaxFreqInd]:
        plt.plot(FreqTrunc, AbsFitList[i], '--')
    plt.xlabel('freq/GHz', fontsize='x-large')
    plt.ylabel('Abs', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plotInd = 0
    leg = ()
    ax.errorbar(Time, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
    plt.xlabel('Time(ns)', fontsize='x-large')
    plt.ylabel('f0(GHz)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()


    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plotInd = 1
    leg = ()
    ax.errorbar(Time, OptMatrix[plotInd, :] * 1e3, yerr=ErrMatrix[plotInd, :] * 1e3, fmt='o')
    plt.xlabel('Time(ns)', fontsize='x-large')
    plt.ylabel('kappa(MHz)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

plt.show()
