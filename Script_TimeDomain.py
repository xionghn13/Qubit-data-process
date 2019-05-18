import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import MeasurementControlFunc as mcf


def TimeDomainMeasurement(Current, QubitFreq, DrivingPower):
    # ExperimentName = 'wg6 in 7.5GHz cavity'
    # CoolDownDate = 'test'

    DrivingFreq = QubitFreq
    ItemDict = {
        'Yoko - Current': Current,
        'Pump - Frequency': DrivingFreq,
        'Pump - Power': DrivingPower
    }

    ConfigName = 'rabi.hdf5'
    MeasLabel = 'rabi'
    mcf.RunMeasurement(ConfigName, MeasLabel)


if __name__ == '__main__':
    Current = 6.04e-3
    QubitFreq = 293.27e6
    DrivingPower = 0
    TimeDomainMeasurement(Current, QubitFreq, DrivingPower)
