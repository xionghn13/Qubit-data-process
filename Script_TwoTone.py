import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import ExtractDataFunc as edf
import MeasurementControlFunc as mcf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def FindQubitFreqTwoTone(Current, Anchor1, Anchor2, Avg=500e3, Power=0, ReadoutFreq=7.292e9, ReadoutPower=0, Span=10e6,
                         SeqLen=1e4, SaveFig=True,
                         DataFolderName='10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'):
    ConfigName = 'two tone sweep for cavity.hdf5'
    MeasLabel = 'two tone'

    if (Anchor1[0] - Anchor2[0]) == 0:
        slope = 0
    else:
        slope = (Anchor1[1] - Anchor2[1]) / (Anchor1[0] - Anchor2[0])
    PredictFreq = slope * (Current - Anchor1[0]) + Anchor1[1]

    ItemDict = {
        'Drive2 - Frequency': [[PredictFreq - Span / 2, 'START'], [PredictFreq + Span / 2, 'STOP']],
        'Drive2 - Power': Power,
        'Readout - Power': ReadoutPower,
        'Readout - Frequency': ReadoutFreq,
        'Alazar - Number of records': Avg,
        'Yoko - Current': Current,
        'Counter - Number of points': 0,
        'Pulse Generator - Number of points': SeqLen,
        'Alazar - Channel B - Range': 4  # 200mV
    }
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict, DataFolderName=DataFolderName)
    # OutFile = 'two tone_2019-06-04-22-34-04.hdf5'
    # OutPath = 'C:/Users/admin\\Labber\\Data/2019/06\\Data_0604/'
    [Freq, Complex] = edf.readFSweepTwoToneLabber(OutPath + OutFile)

    y_data = np.angle(Complex)
    y_data = np.abs(Complex)
    N = 4
    y_data_conv = np.convolve(y_data, np.ones((N,)) / N, mode='valid')
    Freq_conv = np.convolve(Freq, np.ones((N,)) / N, mode='valid')
    kappa_guess = Span / 4 / 1e9
    C_guess = (y_data_conv[-1] - y_data_conv[0]) / (Freq_conv[-1] - Freq_conv[0])
    Ind = np.argmax(np.abs(y_data_conv - C_guess * (Freq_conv - Freq_conv[0]) - y_data_conv[0]))
    f0_guess = Freq_conv[Ind]
    B_guess = y_data_conv[0] + C_guess * (f0_guess - Freq_conv[0])
    A_guess = (y_data_conv[Ind] - B_guess) * (kappa_guess / 2) ** 2

    def lorenztian(f, f0, kappa, A, B, C):
        t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B + C * (f - f0)
        return t

    guess = ([f0_guess, kappa_guess, A_guess, B_guess, C_guess])
    bounds = (
        (Freq[0], 0, - A_guess.__abs__() * 10, - B_guess.__abs__() * 10, - np.abs(C_guess) * 10),
        (Freq[-1], kappa_guess * 4, A_guess.__abs__() * 10, B_guess.__abs__() * 10, np.abs(C_guess) * 10)
    )
    # print(guess)
    # print(bounds)
    qopt, qcov = curve_fit(lorenztian, Freq, y_data, guess, bounds=bounds)
    f0_fit, kappa_fit, A_fit, B_fit, C_fit = qopt
    f0_std, kappa_std, A_std, B_std, C_std = np.sqrt(qcov.diagonal())

    CurveGuess = lorenztian(Freq, f0_guess, kappa_guess, A_guess, B_guess, C_guess)
    CurveFit = lorenztian(Freq, f0_fit, kappa_fit, A_fit, B_fit, C_fit)

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plt.plot(Freq, y_data)
    plt.plot(Freq, CurveFit, 'r.')
    plt.plot(Freq, CurveGuess, 'y.')
    plt.xlabel('freq/GHz', fontsize='x-large')
    # plt.ylabel('Arb.', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.title('f0=%.4GGHz, kappa/2pi=%.3GMHz, A=%.3G, B=%.3G' % (f0_fit, kappa_fit * 1e3, A_fit, B_fit))
    plt.tight_layout()
    if SaveFig:
        FigPath = OutPath + 'figures/'
        if not os.path.exists(FigPath):
            os.makedirs(FigPath)
        FigName = OutFile.split('.')[0] + '.PNG'
        plt.savefig(FigPath + FigName)
        plt.close('all')
    else:
        plt.show()

    return [f0_fit, f0_std]


def FindReadoutFreqOneTone(Current, Anchor1, Anchor2, Avg=10e3, ReadoutPower=5, Span=10e6,
                           SeqLen=2e4, SaveFig=True,
                           DataFolderName='10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'):
    ConfigName = 'one tone sweep for script.hdf5'
    MeasLabel = 'one tone'

    if (Anchor1[0] - Anchor2[0]) == 0:
        slope = 0
    else:
        slope = (Anchor1[1] - Anchor2[1]) / (Anchor1[0] - Anchor2[0])
    PredictFreq = slope * (Current - Anchor1[0]) + Anchor1[1]

    ItemDict = {
        'Readout - Power': ReadoutPower,
        'Readout - Frequency': [[PredictFreq - Span / 2, 'START'], [PredictFreq + Span / 2, 'STOP']],
        'Alazar - Number of records': Avg,
        'Yoko - Current': Current,
        # 'Counter - Number of points': 0,
        'Pulse Generator - Number of points': SeqLen,
        'Alazar - Channel B - Range': 4  # 200mV
    }
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict,
                                            DataFolderName=DataFolderName)
    [Freq, Complex] = edf.readFSweepLabber(OutPath + OutFile)

    # y_data = np.angle(Complex)
    y_data = np.abs(Complex)
    kappa_guess = Span / 4 / 1e9
    C_guess = (y_data[-1] - y_data[0]) / (Freq[-1] - Freq[0])
    Ind = np.argmax(np.abs(y_data - C_guess * (Freq - Freq[0]) - y_data[0]))
    f0_guess = Freq[Ind]
    B_guess = y_data[0] + C_guess * (f0_guess - Freq[0])
    A_guess = (y_data[Ind] - B_guess) * (kappa_guess / 2) ** 2

    def lorenztian(f, f0, kappa, A, B, C):
        t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B + C * (f - f0)
        return t

    guess = ([f0_guess, kappa_guess, A_guess, B_guess, C_guess])
    bounds = (
        (Freq[0], 0, - A_guess.__abs__() * 10, - B_guess.__abs__() * 10, - np.abs(C_guess) * 10),
        (Freq[-1], kappa_guess * 4, A_guess.__abs__() * 10, B_guess.__abs__() * 10, np.abs(C_guess) * 10)
    )
    # print(guess)
    # print(bounds)
    qopt, qcov = curve_fit(lorenztian, Freq, y_data, guess, bounds=bounds)
    f0_fit, kappa_fit, A_fit, B_fit, C_fit = qopt

    CurveGuess = lorenztian(Freq, f0_guess, kappa_guess, A_guess, B_guess, C_guess)
    CurveFit = lorenztian(Freq, f0_fit, kappa_fit, A_fit, B_fit, C_fit)

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    plt.plot(Freq, y_data)
    plt.plot(Freq, CurveFit, 'r.')
    plt.plot(Freq, CurveGuess, 'y.')
    plt.xlabel('freq/GHz', fontsize='x-large')
    # plt.ylabel('Arb.', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.title('f0=%.4GGHz, kappa/2pi=%.3GMHz, A=%.3G, B=%.3G' % (f0_fit, kappa_fit * 1e3, A_fit, B_fit))
    plt.tight_layout()
    if SaveFig:
        FigPath = OutPath + 'figures/'
        if not os.path.exists(FigPath):
            os.makedirs(FigPath)
        FigName = OutFile.split('.')[0] + '.PNG'
        plt.savefig(FigPath + FigName)
        plt.close('all')
    else:
        plt.show()

    return [f0_fit, kappa_fit]


if __name__ == '__main__':
    # ExperimentName = 'wg6 in 7.5GHz cavity'
    # CoolDownDate = 'test'
    DataFolderName = '10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'
    Current = 1.9e-3
    AnchorCurrents = [2.54e-3, 2.542e-3]
    ExpectFreq = [2.128e9, 7.978e9]
    Anchor1 = [AnchorCurrents[0], ExpectFreq[0]]
    Anchor2 = [AnchorCurrents[1], ExpectFreq[0]]
    AnchorReadout1 = [AnchorCurrents[0], ExpectFreq[1]]
    AnchorReadout2 = [AnchorCurrents[1], ExpectFreq[1]]
    SaveFig = False
    Span = 200e6
    Power = -15
    Avg = 10e3
    SeqLen = 10e3
    ReadoutPower = -10
    [fc, width] = FindReadoutFreqOneTone(Current, AnchorReadout1, AnchorReadout2, SaveFig=SaveFig,
                                         ReadoutPower=ReadoutPower, DataFolderName=DataFolderName)
    fc *= 1e9
    width *= 1e9
    ReadoutFreq = fc - width / 2
    # ReadoutFreq = 7.967e9
    f0 = FindQubitFreqTwoTone(Current, Anchor1, Anchor2, Avg=Avg, Power=Power, Span=Span, SeqLen=SeqLen,
                              ReadoutFreq=ReadoutFreq,
                              SaveFig=SaveFig, ReadoutPower=ReadoutPower, DataFolderName=DataFolderName)
    print(f0)
