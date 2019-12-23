from QubitDataProcessPackages import *

# figure 2

DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process/paper_figures\Fluorescence/'
BackgroundFile = 'power spectroscopy_76.hdf5'
OneToneFile = 'power spectroscopy_77.hdf5'

Calibration = True
UseOneToneRange = True

StartFreq = 6.542
EndFreq = 6.554
StartPower = -10
EndPower = 10
SelectPower = np.array([])
# SelectPower = np.array([-25, -20, -15, -10])
########################################################################################################################
NotUsePowerSpectroscopyCalibrate = BackgroundFile.startswith('one_tone') or BackgroundFile.startswith('calibration')

[BackFreq, BackPower, BackComplex] = edf.readFPSweepLabber(DataPath + BackgroundFile)
BackPower = BackPower[0, 0]
BackFreq = BackFreq[:, 0]
BackComplex = BackComplex[:, 0]
BackPowerStr = str(BackPower)

[OneFreq, OnePower, OneComplex] = edf.readFPSweepLabber(DataPath + OneToneFile)
OnePowerUniq = np.unique(OnePower)
NumPower = np.size(OnePowerUniq)
OneFreqUniq = np.unique(OneFreq)
NumFreq = np.size(OneFreqUniq)
if UseOneToneRange:
    StartFreq = OneFreqUniq.min()
    EndFreq = OneFreqUniq.max()
    StartPower = OnePowerUniq.min()
    EndPower = OnePowerUniq.max()
if Calibration:
    RComplex = sbf.FPSweepBackgroundCalibrate(OneFreq, OnePower, OneComplex, BackFreq, BackComplex, BackPower)

OneComplexNormalized = OneComplex * 10 ** (- OnePower / 20)
BackComplexNormalized = BackComplex * 10 ** (- BackPower / 20)
FreqInd = (EndFreq >= OneFreqUniq) == (OneFreqUniq >= StartFreq)
if len(SelectPower) > 0:
    for i, p in enumerate(SelectPower):
        if i == 0:
            PowerInd = OnePowerUniq == p
        else:
            PowerInd = PowerInd + (OnePowerUniq == p)
else:
    PowerInd = (EndPower >= OnePowerUniq) == (OnePowerUniq >= StartPower)
OneFreqUniqTrunc = OneFreqUniq[FreqInd]
OnePowerUniqTrunc = OnePowerUniq[PowerInd]
NumFreqTrunc = len(OneFreqUniqTrunc)
NumPowerTrunc = len(OnePowerUniqTrunc)
OneComplexNormalizedTrunc = OneComplexNormalized[FreqInd, :]
OneComplexNormalizedTrunc = OneComplexNormalizedTrunc[:, PowerInd]

BackFreqUniq = np.unique(BackFreq)
BackFreqInd = (EndFreq > BackFreqUniq) == (BackFreqUniq > StartFreq)
BackFreqUniqTrunc = BackFreqUniq[BackFreqInd]
BackComplexNormalizedTrunc = BackComplexNormalized[BackFreqInd]
RComplexTrunc = RComplex[FreqInd, :]
RComplexTrunc = RComplexTrunc[:, PowerInd]

f = h5py.File(DataPath + 'circles_data.hdf5', 'w')
fset = f.create_dataset('freq', data=OneFreqUniqTrunc)
pset = f.create_dataset('power', data=OnePowerUniqTrunc)
r_re_set = f.create_dataset('R_real', data=RComplexTrunc.real)
r_im_set = f.create_dataset('R_imag', data=RComplexTrunc.imag)
f.close()

# figure 3

BackgroundFile = 'RefRabiCal_4.hdf5'
RabiFile = 'RefRabi_12.hdf5'
PopulationConversionConst = [1, 1.]
# use rabi to calibrate rabi. Plus or minus are merged in one array.
[BackFreq, BackPower, BackComplex] = edf.readRefRabiCalLabber(DataPath + BackgroundFile)
Plus50MHzBackPower = BackPower
Minus50MHzBackPower = BackPower
Plus50MHzBackFreq = BackFreq
Plus50MHzBackComplex = BackComplex
Minus50MHzBackFreq = Plus50MHzBackFreq
Minus50MHzBackComplex = Plus50MHzBackComplex
# read data
ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
MeasurementType = 'rabi'
[Time, ReadoutLowerFreq, ReadoutHigherFreq, ComplexLowerFreq, ComplexHigherFreq] = edf.readRefRabiLabber(
    DataPath + RabiFile)
RComplexLowerFreq = sbf.FPSweepBackgroundCalibrate(ReadoutLowerFreq, ReadoutPower, ComplexLowerFreq,
                                                   Plus50MHzBackFreq,
                                                   Plus50MHzBackComplex, Plus50MHzBackPower)
RComplexHigherFreq = sbf.FPSweepBackgroundCalibrate(ReadoutHigherFreq, ReadoutPower, ComplexHigherFreq,
                                                    Minus50MHzBackFreq,
                                                    Minus50MHzBackComplex, Minus50MHzBackPower)
ComplexLowerFreqNormalized = ComplexLowerFreq * 10 ** (- ReadoutPower / 20)
ComplexHigherFreqNormalized = ComplexHigherFreq * 10 ** (- ReadoutPower / 20)

y_data = np.array(RComplexLowerFreq.real, dtype='float64')
population = (PopulationConversionConst[0] - y_data) * PopulationConversionConst[1]
f = h5py.File(DataPath + 'rabi_data.hdf5', 'w')
yset = f.create_dataset('y', data=population)
xset = f.create_dataset('x', data=Time)
f.close()

# t1 t2
RabiFileList = [
    'ref_t1_ramsey_echo_interleaved_1.hdf5',
]

NumFile = len(RabiFileList)
CounterArray = np.zeros([NumFile, ])
# analyze background file

RabiFile = 'ref_t1_ramsey_echo_interleaved_1.hdf5'
RabiFileStrList = RabiFile[:-5].split('_')
MeasurementType = RabiFileStrList[0]
MeasurementType = 'ref_t1_ramsey_echo_interleaved'
ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)

[Time, Counter, ReadoutFreq, ComplexVoltageT1, ComplexVoltageRamsey,
 ComplexVoltageEcho] = edf.readReferencedRepeatedT1RamseyEchoInterleavedSweepLabber(
    DataPath + RabiFile)

ComplexVoltageT1Normalized = ComplexVoltageT1 * 10 ** (- ReadoutPower / 20)
ComplexVoltageRamseyNormalized = ComplexVoltageRamsey * 10 ** (- ReadoutPower / 20)
ComplexVoltageEchoNormalized = ComplexVoltageEcho * 10 ** (- ReadoutPower / 20)
RComplexT1 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT1, BackFreq,
                                            BackComplex, BackPower)
RComplexRamsey = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageRamsey,
                                                BackFreq, BackComplex, BackPower)
RComplexEcho = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageEcho,
                                              BackFreq, BackComplex, BackPower)

y_dataList = [np.array(RComplexT1[:, 0, 0].real, dtype='float64'),
              np.array(RComplexRamsey[:, 0, 0].real, dtype='float64'),
              np.array(RComplexEcho[:, 0, 0].real, dtype='float64')]
for i in range(3):
    y_dataList[i] = (PopulationConversionConst[0] - y_dataList[i]) * PopulationConversionConst[1]
x_data = np.array(Time, dtype='float64')

f = h5py.File(DataPath + 't1_ramsey_echo_data.hdf5', 'w')
yset = f.create_dataset('y', data=y_dataList)
xset = f.create_dataset('x', data=x_data)
f.close()

# figure 4

BackgroundFile = 'power spectroscopy_101.hdf5'
RabiFileList = [
    'transient_30.hdf5',
    'transient_29.hdf5',
    'transient_27.hdf5',
]

NumFile = len(RabiFileList)

[BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
BackPower = edf.readReadoutPowerLabber(DataPath + BackgroundFile)
BackPowerStr = str(BackPower)

for i, RabiFile in enumerate(RabiFileList):
    RabiFileStrList = RabiFile[:-5].split('_')
    MeasurementType = RabiFileStrList[0]

    ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
    ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
    QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)
    [Time, DrivePower, ComplexRabi] = edf.readRabiPowerSweepLabber(DataPath + RabiFile)

    ComplexRabiNormalized = ComplexRabi * 10 ** (- ReadoutPower / 20)
    RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexRabi, BackFreq, BackComplex, BackPower)

    for j, power in enumerate(DrivePower):
        if i == 0 and j == 0:
            TimeList = []
            RComplexList = []
            FitRList = []
            DrivePowerArray = np.array([power])
        else:
            DrivePowerArray = np.concatenate((DrivePowerArray, np.array([power])))
        TimeList.append(Time)
        RComplexList.append(RComplex[:, j].real)

f = h5py.File(DataPath + 'transient.hdf5', 'w')
timeset = f.create_dataset('time', data=TimeList)
Rset = f.create_dataset('R', data=RComplexList)
Pset = f.create_dataset('power', data=DrivePowerArray)
f.close()
#
#
#
#
# Sup:
# 240uA
# BackgroundFile = 'power spectroscopy_76.hdf5'
# OneToneFile = 'power spectroscopy_82.hdf5'
#
# 01 rabi
# BackgroundFile = 'power spectroscopy_76.hdf5'
# RabiFile = 'rabi_14.hdf5'
#
# p0/p1=3.094827586206897
# T = 49.2mK
# 02 rabi
# BackgroundFile = 'power spectroscopy_76.hdf5'
# RabiFile = 'rabi_12.hdf5'
#
# p0/p2=25.068965517241356
# T = 57.9mK
# p0=0.734
#
# transient
#
