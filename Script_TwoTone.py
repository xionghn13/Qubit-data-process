import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import ExtractDataFunc as edf
import MeasurementControlFunc as mcf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




def FindQubitFreqTwoTone(Current, Anchor1, Anchor2, Power=0, Span=10e6, SaveFig=True):
    ConfigName = 'two tone sweep for cavity.hdf5'
    MeasLabel = 'two tone'

    slope = (Anchor1[1] - Anchor2[1]) / (Anchor1[0] - Anchor2[0])
    PredictFreq = slope * (Current - Anchor1[0]) + Anchor1[1]

    ItemDict = {
        'Pump - Frequency': [[PredictFreq - Span / 2, 'START'], [PredictFreq + Span / 2, 'STOP']],
        'Pump - Power': Power,
        'Yoko - Current': Current,
    }

    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict)
    # OutFile = 'two tone_2019-06-04-22-34-04.hdf5'
    # OutPath = 'C:/Users/admin\\Labber\\Data/2019/06\\Data_0604/'
    [Freq, Complex] = edf.readFSweepTwoToneLabber(OutPath + OutFile)

    y_data = np.angle(Complex)
    f0_guess = PredictFreq / 1e9
    kappa_guess = Span / 4 / 1e9
    B_guess = y_data[0] / 2 + y_data[-1] / 2
    Ind = np.argmax(np.abs(y_data - B_guess))
    A_guess = (y_data[Ind] - B_guess) * (kappa_guess / 2) ** 2

    def lorenztian(f, f0, kappa, A, B):
        t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B
        return t

    guess = ([f0_guess, kappa_guess, A_guess, B_guess])
    bounds = (
        (Freq[0], 0, - A_guess.__abs__() * 10, - B_guess.__abs__() * 10),
        (Freq[-1], kappa_guess * 4, A_guess.__abs__() * 10, B_guess.__abs__() * 10)
    )
    # print(guess)
    # print(bounds)
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
    plt.title('f0=%.4GGHz, kappa/2pi=%.3GMHz, A=%.3G, B=%.3G' % (f0_fit, kappa_fit * 1e3, A_fit, B_fit))
    plt.tight_layout()
    if SaveFig:
        FigPath = OutPath + 'figures/'
        if not os.path.exists(FigPath):
            os.makedirs(FigPath)
        FigName = OutFile.split('.')[0] + '.PNG'
        plt.savefig(FigPath + FigName)
    else:
        plt.show()

    return f0_fit

if __name__ == '__main__':
    # ExperimentName = 'wg6 in 7.5GHz cavity'
    # CoolDownDate = 'test'
    Current = 6.2e-3
    Anchor1 = [6.22e-3, 581e6]
    Anchor2 = [6.23e-3, 572.2e6]
    SaveFig = False
    Span = 10e6
    Power = 0
    f0 = FindQubitFreqTwoTone(Current, Anchor1, Anchor2, Power=Power, Span=Span, SaveFig=SaveFig)
    print(f0)