import numpy as np
import scipy.interpolate as itp
import ExtractDataFunc as edf


def FISweepSelfCalibrate(DataPath, OneToneFile):
    if OneToneFile.endswith('.dat'):
        [OneFreq, OneCurrent, OneComplex] = edf.readFISweepDat(DataPath + OneToneFile)
    elif OneToneFile.endswith('.hdf5'):
        if OneToneFile.startswith('one tone'):
            [OneFreq, OneCurrent, OneComplex] = edf.readFISweepLabber(DataPath + OneToneFile)
        else:
            [OneFreq, OneCurrent, OneComplex] = edf.readFISweepTwoToneLabber(DataPath + OneToneFile)
    OneFreqUniq = np.unique(OneFreq)
    OneCurrentUniq = np.unique(OneCurrent)
    if OneToneFile.startswith('two tone'):
        OneComplex = (OneComplex / OneComplex.mean(axis=0))
    else:
        OneComplex = np.transpose(OneComplex.transpose() / OneComplex.mean(axis=1))
    return [OneFreqUniq, OneCurrentUniq, OneComplex]


def FISweepBackgroundCalibrate(DataPath, OneToneFile, BackgroundFile, OneTonePower):
    [OneFreq, OneCurrent, OneComplex] = edf.readFISweepDat(DataPath + OneToneFile)
    [BackFreq, BackComplex] = edf.readFSweepDat(DataPath + BackgroundFile)
    BackPowerStr = BackgroundFile.split('_')[5][:-3]
    BackPower = float(BackPowerStr)
    BackComplexITP = itp.interp1d(BackFreq, BackComplex)
    RComplex = OneComplex / BackComplexITP(OneFreq) * 10 ** (BackPower / 20 - OneTonePower / 20)
    OneFreqUniq = np.unique(OneFreq)
    OneCurrentUniq = np.unique(OneCurrent)
    return [OneFreqUniq, OneCurrentUniq, RComplex]



def FPSweepBackgroundCalibrate(OneFreq, OnePower, OneComplex, BackFreq, BackComplex, BackPower):
    if isinstance(BackPower, float):
        if len(BackFreq) > 1:
            BackComplexITP = itp.interp1d(BackFreq, BackComplex)
            RComplex = OneComplex / BackComplexITP(OneFreq) * 10 ** (BackPower / 20 - OnePower / 20)
        else:
            RComplex = OneComplex / BackComplex * 10 ** (BackPower / 20 - OnePower / 20)
    else:
        BackFreqUniq = np.unique(BackFreq)
        OneFreqUniq = np.unique(OneFreq)
        BackPowerUniq = np.unique(BackPower)
        OnePowerUniq = np.unique(OnePower)
        if np.all(BackPowerUniq != OnePowerUniq):
            raise ValueError('The power of the background does not match the one tone power.')
        RComplex = OneComplex * 0
        for i in range(len(BackPowerUniq)):
            BackComplexITP = itp.interp1d(BackFreqUniq, BackComplex[:, i])
            RComplex[:, i] = OneComplex[:, i] / BackComplexITP(OneFreqUniq)
    return RComplex

def FSweepPhaseCalibrate(Freq, RComplex):

    return

def dataListRange(NoBackgroundDataList):
    [OneCurrentUniqList, OneFreqUniqList, OneComplex3List] = NoBackgroundDataList
    NumFile = len(OneCurrentUniqList)
    for i in range(NumFile):
        if i == 0:
            Imin = OneCurrentUniqList[i].min()
            Imax = OneCurrentUniqList[i].max()
            freqmin = OneFreqUniqList[i].min()
            freqmax = OneFreqUniqList[i].max()
        else:
            Imin = min(Imin, OneCurrentUniqList[i].min())
            Imax = max(Imax, OneCurrentUniqList[i].max())
            freqmin = min(freqmin, OneFreqUniqList[i].min())
            freqmax = max(freqmax, OneFreqUniqList[i].max())
    return [[Imin, Imax], [freqmin, freqmax]]
