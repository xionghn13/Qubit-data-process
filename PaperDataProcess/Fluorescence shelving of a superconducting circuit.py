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
if Calibration:
    if BackgroundFile == []:
        [Plus50MHzBackFreq, Plus50MHzBackComplex] = edf.readFSweepDat(DataPath + Plus50MHzBackgroundFile)
        Plus50MHzBackPowerStr = Plus50MHzBackgroundFile.split('_')[5][:-3]
        Plus50MHzBackPower = float(Plus50MHzBackPowerStr)

        [Minus50MHzBackFreq, Minus50MHzBackComplex] = edf.readFSweepDat(DataPath + Minus50MHzBackgroundFile)
        Minus50MHzBackPowerStr = Minus50MHzBackgroundFile.split('_')[5][:-3]
        Minus50MHzBackPower = float(Minus50MHzBackPowerStr)
    elif BackgroundFile.startswith('one_tone'):
        [Plus50MHzBackFreq, Plus50MHzBackComplex] = edf.readFSweepDat(BackgroundFolder + BackgroundFile)
        Plus50MHzBackPowerStr = BackgroundFile.split('_')[5][:-3]
        Plus50MHzBackPower = float(Plus50MHzBackPowerStr)

        [Minus50MHzBackFreq, Minus50MHzBackComplex] = [Plus50MHzBackFreq, Plus50MHzBackComplex]
        Minus50MHzBackPowerStr = Plus50MHzBackPowerStr
        Minus50MHzBackPower = Plus50MHzBackPower
    elif BackgroundFile.endswith('dat'):
        BackgroundFileStrList = BackgroundFile.split('_')
        MeasurementType = BackgroundFileStrList[1]
        if MeasurementType == 'rabi' and BackgroundFileStrList[3] != 'no pump':
            MeasurementType = 'pump rabi'
        BackFreqDict = {'t1': 3, 'rabi': 5, 'pump rabi': 6}
        BackFreqInd = BackFreqDict[MeasurementType]
        BackPowerDict = {'t1': 5, 'rabi': 7, 'pump rabi': 8}
        BackPowerInd = BackPowerDict[MeasurementType]
        BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
        BackLowerFreq = round(float(BackgroundFileStrList[BackFreqInd][:-3]) - IQModFreq, 4)
        BackHigherFreq = round(BackLowerFreq + IQModFreq * 2, 15)
        BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
        Plus50MHzBackPower = BackPower
        Minus50MHzBackPower = BackPower
        if MeasurementType in ('rabi', 'transient'):
            [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(BackgroundFolder + BackgroundFile)
        elif MeasurementType == 't1':
            [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(BackgroundFolder + BackgroundFile)
        Plus50MHzBackFreq = np.array([BackLowerFreq, BackHigherFreq])
        Plus50MHzBackComplex = np.array([ComplexLowerFreq.mean(), ComplexHigherFreq.mean()])
        Minus50MHzBackFreq = Plus50MHzBackFreq
        Minus50MHzBackComplex = Plus50MHzBackComplex
    elif BackgroundFile.endswith('hdf5'):
        if BackgroundFile.startswith('transient'):
            [Time, ComplexLowerFreq] = edf.readRabiLabber(BackgroundFolder + BackgroundFile)
            BackPower = edf.readReadoutPowerLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            BackFreq = edf.readReadoutFreqLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackFreq = np.array([BackFreq])
            Plus50MHzBackComplex = np.array([ComplexLowerFreq.mean()])
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex
        elif BackgroundFile.startswith('calibration') or BackgroundFile.startswith('power'):
            [BackFreq, BackComplex] = edf.readFSweepLabber(BackgroundFolder + BackgroundFile)
            BackPower = edf.readReadoutPowerLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            Plus50MHzBackFreq = BackFreq
            Plus50MHzBackComplex = BackComplex
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex
        elif BackgroundFile.startswith('RefRabiCal'):
            # use rabi to calibrate rabi. Plus or minus are merged in one array.
            [BackFreq, BackPower, BackComplex] = edf.readRefRabiCalLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            Plus50MHzBackFreq = BackFreq
            Plus50MHzBackComplex = BackComplex
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex

# read data file
if RabiFile.endswith('dat'):
    RabiFileStrList = RabiFile.split('_')
    MeasurementType = RabiFileStrList[1]
    if MeasurementType == 'rabi' and RabiFileStrList[3] != 'no pump':
        MeasurementType = 'pump rabi'
    if MeasurementType == 'Ch1':
        MeasurementType = 'Ch1 rabi'
    if MeasurementType == 'Ch2 rabi 2Vpp':
        MeasurementType = 'Ch1 pump rabi'

    ReadoutFreqDict = {'t1': 3, 'rabi': 5, 'Ch1 rabi': 6, 'pump rabi': 6, 'Ch1 pump rabi': 4, 't2': 3, 't2echo': 3}
    ReadoutFreqInd = ReadoutFreqDict[MeasurementType]
    ReadoutPowerDict = {'t1': 5, 'rabi': 7, 'Ch1 rabi': 8, 'pump rabi': 8, 'Ch1 pump rabi': 6, 't2': 5, 't2echo': 5}
    ReadoutPowerInd = ReadoutPowerDict[MeasurementType]
    ReadoutFreq = round(float(RabiFileStrList[ReadoutFreqInd][:-3]), 4)
    ReadoutLowerFreq = round(float(RabiFileStrList[ReadoutFreqInd][:-3]) - IQModFreq, 4)
    ReadoutHigherFreq = round(ReadoutLowerFreq + IQModFreq * 2, 15)
    ReadoutPower = float(RabiFileStrList[ReadoutPowerInd][:-3])
    QubitFreqDict = {'rabi': 8, 'Ch1 rabi': 9, 't1': 6, 'Ch1 pump rabi': 7, 't2': 6, 't2echo': 7}
    QubitFreqInd = QubitFreqDict[MeasurementType]
    QubitFreq = float(RabiFileStrList[QubitFreqInd][5:-3])
    if QubitFreq in (ReadoutLowerFreq, ReadoutHigherFreq) and MeasurementType == 'rabi':
        MeasurementType = 'transient'

    if MeasurementType in ('rabi', 'Ch1 rabi', 'transient', 'Ch1 pump rabi'):
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(DataPath + RabiFile)
        if LimitTimeRange:
            TimeInd = (EndTime >= Time) == (Time >= StartTime)
            Time = Time[TimeInd]
            ComplexLowerFreq = ComplexLowerFreq[TimeInd]
            ComplexHigherFreq = ComplexHigherFreq[TimeInd]
    elif MeasurementType == 't1':
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(DataPath + RabiFile)
        if LimitTimeRange:
            TimeInd = (EndTime >= Time) == (Time >= StartTime)
            Time = Time[TimeInd]
            ComplexLowerFreq = ComplexLowerFreq[TimeInd]
            ComplexHigherFreq = ComplexHigherFreq[TimeInd]
    elif MeasurementType in ('t2', 't2echo'):
        [Time, Complex] = edf.readT2H5(DataPath + RabiFile)
        if LimitTimeRange:
            TimeInd = (EndTime >= Time) == (Time >= StartTime)
            Time = Time[TimeInd]
            Complex = Complex[TimeInd]
elif RabiFile.endswith('hdf5'):
    ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
    ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
    if RabiFile.startswith('transient'):
        MeasurementType = 'transient no ref'
        [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
    elif RabiFile.startswith('rabi'):
        MeasurementType = 'rabi no ref'
        if RabiFile.startswith('rabi CH1 drive'):
            [Time, Complex] = edf.readRabiCH1DriveLabber(DataPath + RabiFile)
            # [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
        elif RabiFile.startswith('rabi CH1 pumped'):
            [Time, Complex] = edf.readRabiCH1PumpedLabber(DataPath + RabiFile)
        else:
            [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
    elif RabiFile.startswith('int_t1'):
        MeasurementType = 't1 no ref'
        [Time, Complex] = edf.readIntegratedT1Labber(DataPath + RabiFile)
    elif RabiFile.startswith('t1'):
        MeasurementType = 't1 no ref'
        [Time, Complex] = edf.readT1Labber(DataPath + RabiFile)
    elif RabiFile.startswith('t2_ramsey'):
        MeasurementType = 't2'
        [Time, Complex] = edf.readT2Labber(DataPath + RabiFile)
    elif RabiFile.startswith('t2_echo'):
        MeasurementType = 't2echo'
        [Time, Complex] = edf.readT2Labber(DataPath + RabiFile)
    elif RabiFile.startswith('RefRabi'):
        MeasurementType = 'rabi'
        [Time, ReadoutLowerFreq, ReadoutHigherFreq, ComplexLowerFreq, ComplexHigherFreq] = edf.readRefRabiLabber(
            DataPath + RabiFile)

    if LimitTimeRange:
        TimeInd = (EndTime >= Time) == (Time >= StartTime)
        Time = Time[TimeInd]
        if RabiFile.startswith('RefRabi'):
            ComplexLowerFreq = ComplexLowerFreq[TimeInd]
            ComplexHigherFreq = ComplexHigherFreq[TimeInd]
        else:
            Complex = Complex[TimeInd]

if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
    ComplexNormalized = Complex * 10 ** (- ReadoutPower / 20)
    if Calibration:
        RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, Plus50MHzBackFreq,
                                                  Plus50MHzBackComplex, Plus50MHzBackPower)
    else:
        RComplex = ComplexNormalized

    if RotateComplex:
        RComplex = AutoRotate(RComplex)
    y_data = np.array(RComplex.real, dtype='float64')
    # y_data = np.array(RComplex.imag, dtype='float64')

else:
    RComplexLowerFreq = sbf.FPSweepBackgroundCalibrate(ReadoutLowerFreq, ReadoutPower, ComplexLowerFreq,
                                                       Plus50MHzBackFreq,
                                                       Plus50MHzBackComplex, Plus50MHzBackPower)
    RComplexHigherFreq = sbf.FPSweepBackgroundCalibrate(ReadoutHigherFreq, ReadoutPower, ComplexHigherFreq,
                                                        Minus50MHzBackFreq,
                                                        Minus50MHzBackComplex, Minus50MHzBackPower)
    ComplexLowerFreqNormalized = ComplexLowerFreq * 10 ** (- ReadoutPower / 20)
    ComplexHigherFreqNormalized = ComplexHigherFreq * 10 ** (- ReadoutPower / 20)

    [amp_cor_re_fit, amp_cor_im_fit, P0_fit, P0_im_fit] = CorrectionParam
    if CircleCorrection:
        RComplexLowerFreq /= (amp_cor_re_fit + amp_cor_im_fit * 1j)
        RComplexLowerFreq = (RComplexLowerFreq - 1) / (P0_fit + P0_im_fit * 1j) * np.abs(
            P0_fit + P0_im_fit * 1j) + 1

    # y_data = np.array(RComplexLowerFreq.real, dtype='float64')
    if FitCorrectedR:
        y_data = np.array(np.real(RComplexLowerFreq / RComplexHigherFreq), dtype='float64')
    else:
        y_data = np.array(RComplexLowerFreq.real, dtype='float64')

# figure 4

# DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
# BackgroundFile = 'power spectroscopy_101.hdf5'
# RabiFileList = [
#     'transient_28.hdf5',
#     'transient_27.hdf5',
# ]
#
# IQModFreq = 0.05
# # FitForGamma = True
# Gamma_r = 2.6 * np.pi * 2
# FitCorrectedR = False
# LogScale = False
# Calibration = True
# RotateComplex = False
#
# PhaseSlope = 326.7041108065019
# PhaseRefrenceFreq = 4.105
#
# NumFile = len(RabiFileList)
# DrivePowerArray = np.zeros([NumFile, ])
# # analyze background file
#
# [BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
# BackPower = edf.readReadoutPowerLabber(DataPath + BackgroundFile)
# BackPowerStr = str(BackPower)
#
# for i, RabiFile in enumerate(RabiFileList):
#     RabiFileStrList = RabiFile[:-5].split('_')
#     MeasurementType = RabiFileStrList[0]
#
#     ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
#     ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
#     QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)
#     # DrivePower = edf.readPumpPowerLabber(DataPath + RabiFile)
#
#     if MeasurementType in ('rabi', 'transient'):
#         [Time, DrivePower, ComplexRabi] = edf.readRabiPowerSweepLabber(DataPath + RabiFile)
#     elif MeasurementType == 't1':
#         [Time, DrivePower, ComplexRabi] = edf.readT1PowerSweepLabber(DataPath + RabiFile)
#
#     if i == 0:
#         ind = DrivePower > -45
#         DrivePower = DrivePower[ind]
#         ComplexRabi = ComplexRabi[:, ind]
#
#     ComplexRabiNormalized = ComplexRabi * 10 ** (- ReadoutPower / 20)
#     if Calibration:
#         RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexRabi, BackFreq,
#                                                   BackComplex, BackPower)
#     else:
#         RComplex = ComplexRabiNormalized
