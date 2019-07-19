import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import MeasurementControlFunc as mcf
from ReferencedTSweepPlot import plotReferencedTSweep
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot


def timeDomainMeasurement(Current, ReadoutFreq, QubitFreq, DrivingPower, PiPulseLength, ReadoutPower=-5, Detuning=0.5e6,
                          MeasurementType='rabi', T1MaxDelay=150e-6, Avg=300e3, DutyCyclePoints=400e3,
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
            'Qubit - Power': ReadoutPower,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': [[1200e-9, 'STOP'], [31, 'N_PTS']],
            'Pulse Generator - Number of points': DutyCyclePoints,
            'Alazar - Number of records': Avg,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            # 'Pulse Generator - Number of points': 20e3,
        }
    elif MeasurementType == 't2_ramsey':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Qubit - Power': ReadoutPower,
            'Pump - Frequency': [[DrivingFreq, 'START'], [QubitFreq, 'STOP']],
            'Pump - Power': DrivingPower,
            'Pulse Generator - Sequence duration': [[10e-6, 'STOP']],
            'Pulse Generator - Number of points': DutyCyclePoints,
            'Alazar - Number of records': Avg,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
        }
    elif MeasurementType == 't1':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Qubit - Power': ReadoutPower,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Readout delay': [[T1MaxDelay, 'STOP']],
            # 'Pulse Generator - Number of points': DutyCyclePoints,
            'Alazar - Number of records': Avg,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            # 'Counter - Number of points': [[40, 'STOP']]
            'Counter - Number of points': 0
        }
    elif MeasurementType == 't2_echo':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Qubit - Power': ReadoutPower,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Sequence duration': [[105e-6, 'STOP']],
            'Pulse Generator - Number of points': DutyCyclePoints,
            'Alazar - Number of records': Avg * 4,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            # 'Counter - Number of points': [[10, 'STOP']]
            'Counter - Number of points': 0
        }
    elif MeasurementType == 't1_t2_interleaved':
        ItemDict = {
            'Yoko - Current': Current,
            'Qubit - Frequency': ReadoutFreq,
            'Qubit - Power': ReadoutPower,
            'Pump - Frequency': DrivingFreq,
            'Pump - Power': DrivingPower,
            'Pulse Generator - Width #1': PiPulseLength,
            'Pulse Generator - Sequence duration': [[240e-6, 'STOP']],
            'Pulse Generator - Number of points': DutyCyclePoints,
            'Alazar - Number of records': Avg,
            'Pulse Generator - Pulse type': PulseType,  # 0 Gaussian
            'Counter - Number of points': [[10, 'STOP']]
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
                'Counter - Number of points' in ItemDict and isinstance(ItemDict['Counter - Number of points'], list)):

            plotLabberRepeatedTSweepPlot(OutPath, [OutFile])
        else:
            FitDict = plotReferencedTSweep(OutPath, OutFile, ShowFig=ShowFig)
            return FitDict


if __name__ == '__main__':
    Current = 7.596e-3
    QubitFreq = 0.515e9
    ReadoutPower = -15
    print('Current = %.4GmA\nQubitFreq = %.5GGHz' % (Current * 1e3, QubitFreq / 1e9))
    # Current = 6.35e-3
    # QubitFreq = 514.9e6
    DrivingPower = 10
    PiPulseLength = 8000e-9
    ReadoutFreq = 7.8305e9
    T1MaxDelay = 60e-6
    PulseType = 0  # 0 Gaussian, 1 Square
    Avg = 400e3
    T2RamseyDetuning = 0.5e6
    CyclePoints = 400e3
    # PiPulseLength = 500e-6
    # MeasTypeList = ['rabi']
    # MeasTypeList = ['t1', 't2_echo']
    # MeasTypeList = ['t1']
    # MeasTypeList = ['t2_echo']
    # MeasTypeList = ['rabi', 't1', 't2_ramsey', 't2_echo']
    # MeasTypeList = ['t1_t2_interleaved']
    # MeasTypeList = ['rabi', 't2_ramsey', 't1_t2_interleaved']
    MeasTypeList = ['rabi', 't1']
    # MeasType = 't1'
    # for PiPulseLength in np.linspace(400e-6, 1000e-6, 7):
    for MeasType in MeasTypeList:
    # PiPulseLength = 300e-6
    # for DrivingPower in [5, 10, 15]:
        FitDict = timeDomainMeasurement(Current, ReadoutFreq, QubitFreq, DrivingPower, PiPulseLength, ReadoutPower=ReadoutPower,
                                        T1MaxDelay=T1MaxDelay, Detuning=T2RamseyDetuning, DutyCyclePoints=CyclePoints,
                                        Avg=Avg, PulseType=PulseType, MeasurementType=MeasType, FitAndPlot=True,
                                        ShowFig=False)
        if MeasType == 'rabi':
            PiPulseLength = round(FitDict['opt'][3]) * 1e-9
            print('PiPulseLength now is %.3Gs' % PiPulseLength)
        # elif MeasType == 't1':
        #     T1MaxDelay = round(FitDict['opt'][1] * 5 / 1e3) * 1e-6
        #     print('T1MaxDelay now is %.3Gs' % T1MaxDelay)
