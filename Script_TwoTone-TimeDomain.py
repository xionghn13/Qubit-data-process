import Script_TwoTone
import Script_TimeDomain
import MeasurementControlFunc as mcf
from QubitDataProcessPackages import *

WaitHour = 0  # hr
DataFolderName = '10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'
# current list
CurrentList = np.array([-0.102e-3])
# CurrentList = np.concatenate((np.linspace(2.516e-3, 2.56e-3, 23), np.linspace(2.562e-3, 2.602e-3, 11)))
# CurrentList = np.concatenate((np.linspace(2.606e-3, 2.682e-3, 20), np.linspace(2.684e-3, 2.76e-3, 39)))

# CurrentList = np.linspace(2.734e-3, 2.76e-3, 14)

# CurrentList = np.insert(CurrentList, 0, 6.18e-3)
print(CurrentList)

# anchor points
AnchorCurrents = [2.73e-3, 2.732e-3]
Anchor1 = [AnchorCurrents[0], 6.829e9]
Anchor2 = [AnchorCurrents[1], 6.829e9]
AnchorReadout1 = [AnchorCurrents[0], 7.967e9]
AnchorReadout2 = [AnchorCurrents[1], 7.967e9]

# common setting
FixReadoutFreq = True
SaveFig = True
DrivingPower = 12
PiPulseLength = 0.233e-6
ReadoutFreq = 7.959e9
ReadoutPower = 5
T1MaxDelay = 60e-6
PulseType = 0  # 0 Gaussian, 1 Square
Avg = 10e3
T2RamseyDetuning = 0.01e6
CyclePoints = 50e3

# two tone setting
OneToneSpan = 20e6
TwoToneSpan = 20e6
TwoTonePower = -10
TwoToneAvg = 50e3
TwoToneSeqLen = 10e3

# pumped rabi setting
PumpFreq = 6.70635e9
PumpPower = -5
PumpLength = 61e-9
DriveDelayList = [1e-9]

# time domain measurements
MeasTypeList = ['rabi']
# MeasTypeList = ['rabi CH1 drive']
# MeasTypeList = ['rabi CH1 pumped']
# MeasTypeList = ['t2_ramsey', 't1', 't2_echo']
# MeasTypeList = ['t1']
# MeasTypeList = ['rabi', 't2_ramsey']
# MeasTypeList = ['t2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't1', 't2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't1', 't2_echo']
# MeasTypeList = ['rabi', 't1']
# MeasTypeList = ['rabi', 't1_short']
# MeasTypeList = ['rabi', 't2_echo']
# MeasTypeList = []
# MeasTypeList = ['rabi', 't1_t2_interleaved']
# MeasTypeList = ['rabi', 't1_ramsey_echo_interleaved']
# MeasTypeList = ['rabi', 't2_ramsey', 't1_t2_interleaved']

mcf.Wait(WaitHour * 3600)
if 'rabi CH1 pumped' not in MeasTypeList:
    DriveDelayList = [1e-9]
for i, cur in enumerate(CurrentList):
    if not FixReadoutFreq:
        [fc, width] = Script_TwoTone.FindReadoutFreqOneTone(cur, AnchorReadout1, AnchorReadout2, Span=OneToneSpan,
                                                            ReadoutPower=ReadoutPower, SaveFig=SaveFig,
                                                            DataFolderName=DataFolderName)
        fc *= 1e9
        width *= 1e9
        AnchorReadout2 = AnchorReadout1
        AnchorReadout1 = [cur, fc]
        ReadoutFreq = fc - width / 2
        ReadoutFreq = fc
    [f0, f0_err] = Script_TwoTone.FindQubitFreqTwoTone(cur, Anchor1, Anchor2, Power=TwoTonePower,
                                                       ReadoutFreq=ReadoutFreq,
                                                       ReadoutPower=ReadoutPower, Avg=TwoToneAvg, Span=TwoToneSpan,
                                                       SeqLen=TwoToneSeqLen, SaveFig=SaveFig,
                                                       DataFolderName=DataFolderName)
    # f0 = 0.5286e9
    # f0_err = 0
    f0 *= 1e9
    f0_err *= 1e9

    print(u'Current = %.4GmA\nReadoutFreq = %.5GGHz\nQubitFreq = %.5G\u00B1%.3GGHz' % (
        cur * 1e3, ReadoutFreq / 1e9, f0 / 1e9, f0_err / 1e9))
    if f0_err > f0 * 0.5e-4:
        print('Fail to find qubit. Go to next flux point.')
    else:
        Anchor2 = Anchor1
        Anchor1 = [cur, f0]
        for MeasType in MeasTypeList:
            for DriveDelay in DriveDelayList:
                FitDict = Script_TimeDomain.timeDomainMeasurement(cur, ReadoutFreq, f0, DrivingPower, PiPulseLength,
                                                                  T1MaxDelay=T1MaxDelay, Avg=Avg,
                                                                  Detuning=T2RamseyDetuning,
                                                                  DutyCyclePoints=CyclePoints,
                                                                  ReadoutPower=ReadoutPower,
                                                                  PulseType=PulseType, MeasurementType=MeasType,
                                                                  PumpFreq=PumpFreq, PumpPower=PumpPower,
                                                                  PumpLength=PumpLength,
                                                                  DriveDelay=DriveDelay, FitAndPlot=True, ShowFig=False,
                                                                  DataFolderName=DataFolderName)
            if MeasType == 'rabi':
                PiPulseLength = round(FitDict['opt'][3]) * 1e-9
                print('PiPulseLength now is %.3Gs' % PiPulseLength)
                # elif MeasType == 't1':
                #     T1MaxDelay = round(FitDict['opt'][1] * 5 / 1e3) * 1e-6
                #     print('T1MaxDelay now is %.3Gs' % T1MaxDelay)
