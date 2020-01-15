from QubitDataProcessPackages import *

# figure 2

# DataPath = 'C:\SC Lab\GitHubRepositories\Qubit-data-process\PaperDataProcess\Fluorescence shelving of a superconducting circuit\Fluorescence/'
DataPath = 'D:\GitHubRepository\Qubit-data-process\PaperDataProcess\Fluorescence shelving of a superconducting circuit\Fluorescence/'
BackgroundFile = 'power spectroscopy_105.hdf5'
OneToneFile = 'power spectroscopy_108.hdf5'

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

BackgroundFile = 'power spectroscopy_105.hdf5'
RabiFileList = [
    # 'transient_9.hdf5',
    'transient_36.hdf5',
]

PopulationConversionConst = [1, 0.973600492844838]
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

    if RabiFile in ('transient_28.hdf5'):
        ind = DrivePower > -45
        DrivePower = DrivePower[ind]
        ComplexRabi = ComplexRabi[:, ind]
    # elif RabiFile in ('transient_32.hdf5'):
    #     ind = DrivePower > -40
    #     DrivePower = DrivePower[ind]
    #     ComplexRabi = ComplexRabi[:, ind]

    ComplexRabiNormalized = ComplexRabi * 10 ** (- ReadoutPower / 20)
    RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexRabi, BackFreq, BackComplex, BackPower)

    for j, power in enumerate(DrivePower):
        if i == 0 and j == 0:
            TimeList = []
            y_dataList = []
            FitRList = []
            DrivePowerArray = np.array([power])
        else:
            DrivePowerArray = np.concatenate((DrivePowerArray, np.array([power])))
        TimeList.append(list(Time))
        y_dataList.append(list((PopulationConversionConst[0] - RComplex[:, j].real) * PopulationConversionConst[1]))

f = h5py.File(DataPath + 'transient_data.hdf5', 'w')
# timeset0 = f.create_dataset('time0', data=TimeList[:1])
# Yset0 = f.create_dataset('y0', data=y_dataList[:1])
timeset = f.create_dataset('time', data=TimeList)
Yset = f.create_dataset('y', data=y_dataList)
Pset = f.create_dataset('power', data=DrivePowerArray)
f.close()

# population measurement
BackgroundFile = 'power spectroscopy_105.hdf5'
RabiFileList = [
    'transient_P2_P1_24.hdf5',
    'transient_P2_P1_26.hdf5',
    'transient_P2_P1_31.hdf5',
    't1_P2_P1_23.hdf5',
]
OutFileList = [
    'optimal_power_transient_population.hdf5',
    'high_power_transient_population.hdf5',
    'low_power_transient_population.hdf5',
    '02_pi_pulse_decay.hdf5'

]
PopulationConversionConstList = [
    [1., 1 / 1.0789138211341804],
    [1., 1 / 1.0282720690126486],
    [1., 1. / 1.0621686624421236],
    [1, 1. / 0.9761871987220584],
]
[BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
BackPower = edf.readReadoutPowerLabber(DataPath + BackgroundFile)

for i in range(len(RabiFileList)):
    PopulationConversionConst = PopulationConversionConstList[i]
    RabiFile = RabiFileList[i]
    ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
    ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
    TransientPower = edf.readDrive1PowerLabber(DataPath + RabiFile)
    if RabiFile.startswith('transient'):
        [Time, Complex] = edf.readMultiRabiLabber(DataPath + RabiFile)
    else:
        [Time, Complex] = edf.readMultiT1Labber(DataPath + RabiFile)
    ComplexNormalized = Complex * 10 ** (- ReadoutPower / 20)
    RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, BackFreq, BackComplex, BackPower)
    num_curve = RComplex.shape[1]
    Population = (PopulationConversionConst[0] - RComplex.real) * PopulationConversionConst[1]
    CorrectP2 = True
    if CorrectP2:
        P2PiPulse = 108
        P2RabiT1 = 604
        k = np.exp(- P2PiPulse / P2RabiT1)
        P0 = np.mean(Population[:, [0, 2]], axis=1)
        P2 = Population[:, 3]
        Population[:, 3] = 2 / (k + 1) * (P2 + (k - 1) / 2 * P0)

    P0 = np.mean(Population[:, [0, 2]], axis=1)
    P1 = Population[:, 1]
    P2 = Population[:, 3]
    PList = [P0, P1, P2]

    f = h5py.File(DataPath + OutFileList[i], 'w')
    timeset = f.create_dataset('time', data=Time)
    Yset = f.create_dataset('y', data=PList)
    Pset = f.create_dataset('power', data=TransientPower)
    f.close()

#
#
#
# Sup:
# 240uA
BackgroundFile = 'power spectroscopy_76.hdf5'
OneToneFile = 'power spectroscopy_82.hdf5'

Calibration = True
UseOneToneRange = True

StartFreq = 6.542
EndFreq = 6.554
StartPower = -5
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

f = h5py.File(DataPath + 'sup_circles_data.hdf5', 'w')
fset = f.create_dataset('freq', data=OneFreqUniqTrunc)
pset = f.create_dataset('power', data=OnePowerUniqTrunc)
r_re_set = f.create_dataset('R_real', data=RComplexTrunc.real)
r_im_set = f.create_dataset('R_imag', data=RComplexTrunc.imag)
f.close()

# 01 rabi
BackgroundFile = 'power spectroscopy_76.hdf5'
RabiFile = 'rabi_14.hdf5'
# PopulationConversionConst = [1, 1.]
[BackFreq, BackPower, BackComplex] = edf.readFPSweepLabber(DataPath + BackgroundFile)
BackPower = BackPower[0, 0]
BackFreq = BackFreq[:, 0]
BackComplex = BackComplex[:, 0]

# read data
ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
[Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, BackFreq, BackComplex, BackPower)

y_data = np.array(RComplex.real, dtype='float64')
# population = (PopulationConversionConst[0] - y_data) * PopulationConversionConst[1]
f = h5py.File(DataPath + 'sup_01_rabi_data.hdf5', 'w')
yset = f.create_dataset('y', data=y_data)
xset = f.create_dataset('x', data=Time)
f.close()
# 02 rabi

BackgroundFile = 'power spectroscopy_76.hdf5'
RabiFile = 'rabi_12.hdf5'
# PopulationConversionConst = [1, 1.]
[BackFreq, BackPower, BackComplex] = edf.readFPSweepLabber(DataPath + BackgroundFile)
BackPower = BackPower[0, 0]
BackFreq = BackFreq[:, 0]
BackComplex = BackComplex[:, 0]

# read data
ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
[Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, BackFreq, BackComplex, BackPower)

y_data = np.array(RComplex.real, dtype='float64')
# population = (PopulationConversionConst[0] - y_data) * PopulationConversionConst[1]
f = h5py.File(DataPath + 'sup_02_rabi_data.hdf5', 'w')
yset = f.create_dataset('y', data=y_data)
xset = f.create_dataset('x', data=Time)
f.close()
