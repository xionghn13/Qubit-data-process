import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools
import ExtractDataFunc as edf
import MeasurementControlFunc as mcf
import matplotlib.pyplot as plt


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
fig, ax = plt.subplots()
plt.plot(Freq, Complex.imag)
# plt.xlabel('Power (dBm)', fontsize='x-large')
# plt.ylabel('T_pi (ns)', fontsize='x-large')
# plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()