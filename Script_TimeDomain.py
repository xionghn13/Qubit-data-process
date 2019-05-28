import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import MeasurementControlFunc as mcf
from ReferencedTSweepPlot import plotReferencedTSweep
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot


def timeDomainMeasurement(Current, QubitFreq, DrivingPower, PiPulseLength, Detuning=0.5e6, MeasurementType='rabi',
                          FitAndPlot=True, ShowFig=True):
    # ExperimentName = 'wg6 in 7.5GHz cavity'
    # CoolDownDate = 'test'

    if MeasurementType in ('t2_ramsey'):
        DrivingFreq = QubitFreq - Detuning
    else:
        DrivingFreq = QubitFreq

    if MeasurementType == 'rabi':
        ItemDict = {
            'Yoko - Current': Current,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': [[300e-9, 'STOP']]
        }
    elif MeasurementType == 't2_ramsey':
        ItemDict = {
            'Yoko - Current': Current,
            'Pump - Frequency': [[DrivingFreq, 'START'], [QubitFreq, 'STOP']],
            'Pump - Power': DrivingPower,
            'Pulse Generator - Sequence duration': [[6e-6, 'STOP']]
        }
    elif MeasurementType == 't1':
        ItemDict = {
            'Yoko - Current': Current,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            # 'Pulse Generator - Readout delay': [[450e-6, 'STOP']],
            'Alazar - Number of records': 300e3,
            'Pulse Generator - Pulse type': 1,  # 0 Gaussian
            # 'Counter - Number of points': [[10, 'STOP']]
            'Counter - Number of points': 0
        }
    else:
        ItemDict = {
            'Yoko - Current': Current,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength
        }

    MeasLabel = MeasurementType
    ConfigName = MeasLabel + '.hdf5'
    ConfigName = ConfigName.replace('_', ' ')
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict)
    if FitAndPlot:
        if MeasurementType == 't1_t2_interleaved' or 'Counter - Number of points' in ItemDict and ItemDict['Counter - Number of points'] is list:
            plotLabberRepeatedTSweepPlot(OutPath, OutFile)
        else:
            FitDict = plotReferencedTSweep(OutPath, OutFile, ShowFig=ShowFig)
            return FitDict


if __name__ == '__main__':
    Current = 50.787e-3
    QubitFreq = 109.6e6
    DrivingPower = 0
    PiPulseLength = 126e-9
    # PiPulseLength = 500e-6
    # MeasTypeList = ['rabi']
    MeasTypeList = ['t2_ramsey', 't1', 't2_echo']
    # MeasTypeList = ['t1']
    # MeasTypeList = ['t2echo']
    # MeasTypeList = ['rabi', 't1', 't2_ramsey', 't2_echo']
    # MeasTypeList = ['t1_t2_interleaved']
    # MeasType = 'rabi'
    # MeasType = 't1'
    # MeasType = 't2_echo'
    # MeasType = 't2_ramsey'
    # MeasType = 't1_t2_interleaved'
    # for PiPulseLength in np.linspace(300e-6, 350e-6, 2):
    for MeasType in MeasTypeList:
        FitDict = timeDomainMeasurement(Current, QubitFreq, DrivingPower, PiPulseLength, MeasurementType=MeasType, FitAndPlot=True, ShowFig=False)
        if MeasType == 'rabi':
            PiPulseLength = round(FitDict['opt'][3]) * 1e-9
            print('PiPulseLength now is %.3Gs' % PiPulseLength)
