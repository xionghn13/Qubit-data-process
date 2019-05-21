import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import MeasurementControlFunc as mcf
from ReferencedTSweepPlot import plotReferencedTSweep
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot


def timeDomainMeasurement(Current, QubitFreq, DrivingPower, PiPulseLength, Detuning=0.5e6, MeasurementType='rabi'):
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
            'Pump - Power': DrivingPower
        }
    elif MeasurementType == 't2_ramsey':
        ItemDict = {
            'Yoko - Current': Current,
            'Pump - Frequency': [[DrivingFreq, 'START'], [QubitFreq, 'STOP']],
            'Pump - Power': DrivingPower
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
    if MeasurementType == 't1_t2_interleaved':
        plotLabberRepeatedTSweepPlot(OutPath, OutFile)
    else:
        plotReferencedTSweep(OutPath, OutFile)


if __name__ == '__main__':
    Current = 6.095e-3
    QubitFreq = 109.6e6
    DrivingPower = 0
    PiPulseLength = 126e-9
    # MeasType = 'rabi'
    # MeasType = 't1'
    # MeasType = 't2_echo'
    # MeasType = 't2_ramsey'
    MeasType = 't1_t2_interleaved'
    timeDomainMeasurement(Current, QubitFreq, DrivingPower, PiPulseLength, MeasurementType=MeasType)
