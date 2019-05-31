import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import ExtractDataFunc as edf
import MeasurementControlFunc as mcf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ExperimentName = 'wg6 in 7.5GHz cavity'
# CoolDownDate = 'test'
Current = 5.95e-3
Anchor1 = [6.0419e-3, 283.9e6]
Anchor2 = [6.0619e-3, 195.9e6]

ConfigName = 'two tone sweep for cavity.hdf5'
MeasLabel = 'two tone'

slope = (Anchor1[1] - Anchor2[1]) / (Anchor1[0] - Anchor2[0])
PredictFreq = slope * (Current - Anchor1[0]) + Anchor1[1]
Span = 20e6

ItemDict = {
    'Pump - Frequency': [[PredictFreq - Span / 2, 'START'], [PredictFreq + Span / 2, 'STOP']],
    'Yoko - Current': Current,
}

[OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict)
[Freq, Complex] = edf.readFSweepTwoToneLabber(OutPath + OutFile)

y_data = Complex.imag()
f0_guess = PredictFreq
kappa_guess = Span / 4
B_guess = y_data[0] / 2 + y_data[-1] / 2
Ind = np.argmax(np.abs(y_data - B_guess))
A_guess = (y_data[Ind] - B_guess) * (kappa_guess / 2) ** 2

def lorenztian(f, f0, kappa, A, B):
    t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B
    return t

guess = ([f0_guess, kappa_guess, A_guess, B_guess])
bounds = (
    (Freq[0], 0, - A_guess * 10, - B_guess * 10),
    (Freq[-1], kappa_guess * 4, A_guess * 10, B_guess * 10)
)

qopt, qcov = curve_fit(lorenztian, Freq, y_data, guess, bounds=bounds)
f0_fit, kappa_fit, A_fit, B_fit = qopt

CurveGuess = lorenztian(Freq, f0_guess, kappa_guess, A_guess, B_guess)
CurveFit = lorenztian(Freq, f0_fit, kappa_fit, A_fit, B_fit)

fig, ax = plt.subplots()
plt.plot(Freq, y_data)
plt.plot(Freq, CurveFit, 'r.')
plt.plot(Freq, CurveGuess, 'y.')
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.title('f0=%.3GGHz, kappa/2pi=%.3GMHz, A=%.3G, B=%.3G' % (f0_fit, kappa_fit * 1e3, A_fit, B_fit))
plt.tight_layout()
plt.show()