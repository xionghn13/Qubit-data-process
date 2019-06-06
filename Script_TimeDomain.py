import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import MeasurementControlFunc as mcf
from ReferencedTSweepPlot import plotReferencedTSweep
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot


def timeDomainMeasurement(Current, ReadoutFreq, QubitFreq, DrivingPower, PiPulseLength, Detuning=0.5e6,
                          MeasurementType='rabi', T1MaxDelay=150e-6,
                          PulseType=0, FitAndPlot=True, ShowFig=True):
    # ExperimentName = 'wg6 in 7.5GHz cavity'
    # CoolDownDate = 'test'

    if MeasurementType in ('t2_ramsey'):
        DrivingFreq = QubitFreq - Detuning
    else:
        DrivingFreq = QubitFreq

    if MeasurementType == 'rabi':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': [[1500e-9, 'STOP']],
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
        }
    elif MeasurementType == 't2_ramsey':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Pump - Frequency': [[DrivingFreq, 'START'], [QubitFreq, 'STOP']],
            'Pump - Power': DrivingPower,
            'Pulse Generator - Sequence duration': [[60e-6, 'STOP']],
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
        }
    elif MeasurementType == 't1':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Readout delay': [[T1MaxDelay, 'STOP']],
            'Alazar - Number of records': 300e3,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            # 'Counter - Number of points': [[10, 'STOP']]
            'Counter - Number of points': 0
        }
    elif MeasurementType == 't2_echo':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Sequence duration': [[100e-6, 'STOP']],
            'Alazar - Number of records': 300e3,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            # 'Counter - Number of points': [[10, 'STOP']]
            'Counter - Number of points': 0
        }
    elif MeasurementType == 't1_t2_interleaved':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Sequence duration': [[240e-6, 'STOP']],
            'Alazar - Number of records': 300e3,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            'Counter - Number of points': [[6, 'STOP']]
            # 'Counter - Number of points': 0
        }
    # else:
    #     ItemDict = {
    #         'Yoko - Current': Current,
    #         'Qubit - Frequency': ReadoutFreq,
    #         'Pump - Frequency': DrivingFreq,
    #         'Pump - Power': DrivingPower,
    #         'Pulse Generator - Width #1': PiPulseLength,
    #         'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
    #     }

    MeasLabel = MeasurementType
    ConfigName = MeasLabel + '.hdf5'
    ConfigName = ConfigName.replace('_', ' ')
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict)
    if FitAndPlot:
        if MeasurementType == 't1_t2_interleaved' or (
                'Counter - Number of points' in ItemDict and ItemDict['Counter - Number of points'] is list):

            plotLabberRepeatedTSweepPlot(OutPath, [OutFile])
        else:
            FitDict = plotReferencedTSweep(OutPath, OutFile, ShowFig=ShowFig)
            return FitDict


if __name__ == '__main__':
    Current = 6.27e-3
    QubitFreq = 541.6e6
    print('Current = %.3GmA\nQubitFreq = %.4GGHz' % (Current * 1e3, QubitFreq / 1e9))
    # Current = 6.35e-3
    # QubitFreq = 514.9e6
    DrivingPower = 13
    PiPulseLength = 657e-9
    ReadoutFreq = 8.4205e9
    T1MaxDelay = 500e-6
    PulseType = 0  # 0 Gaussian, 1 Square
    # PiPulseLength = 500e-6
    # MeasTypeList = ['rabi']
    # MeasTypeList = ['t2_ramsey', 't1', 't2_echo']
    # MeasTypeList = ['t1']
    # MeasTypeList = ['t2_echo']
    MeasTypeList = ['rabi', 't1', 't2_ramsey', 't2_echo']
    # MeasTypeList = ['t1_t2_interleaved']
    # MeasTypeList = ['rabi', 't2_ramsey', 't1_t2_interleaved']
    # MeasType = 't1'
    # for PiPulseLength in np.linspace(600e-6, 1000e-6, 5):
    for MeasType in MeasTypeList:
        FitDict = timeDomainMeasurement(Current, ReadoutFreq, QubitFreq, DrivingPower, PiPulseLength, T1MaxDelay=T1MaxDelay,
                                        PulseType=PulseType, MeasurementType=MeasType, FitAndPlot=True, ShowFig=False)
        if MeasType == 'rabi':
            PiPulseLength = round(FitDict['opt'][3]) * 1e-9
            print('PiPulseLength now is %.3Gs' % PiPulseLength)
        # elif MeasType == 't1':
        #     T1MaxDelay = round(FitDict['opt'][1] * 5 / 1e3) * 1e-6
        #     print('T1MaxDelay now is %.3Gs' % T1MaxDelay)

