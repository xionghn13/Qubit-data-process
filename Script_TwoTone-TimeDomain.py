import Script_TwoTone
import Script_TimeDomain
from QubitDataProcessPackages import *

# current list
CurrentList = np.array([2.63e-3])
# CurrentList = np.concatenate((np.linspace(2.456e-3, 2.5e-3, 23), np.linspace(2.502e-3, 2.542e-3, 11)))
# CurrentList = np.concatenate((np.linspace(2.546e-3, 2.622e-3, 20), np.linspace(2.624e-3, 2.7e-3, 39)))

# CurrentList = np.linspace(0.1e-3, 0e-3, 2)
#
# CurrentList = np.insert(CurrentList, 0, 6.18e-3)
print(CurrentList)

# anchor points
AnchorCurrents = [0.3e-3, 0.2e-3]
Anchor1 = [AnchorCurrents[0], 0.5246e9]
Anchor2 = [AnchorCurrents[1], 0.5246e9]
AnchorReadout1 = [AnchorCurrents[0], 7.968e9]
AnchorReadout2 = [AnchorCurrents[1], 7.965e9]

# common setting
FixReadoutFreq = True
SaveFig = True
DrivingPower = 5
PiPulseLength = 25e-9
ReadoutFreq = 7.967e9
ReadoutPower = -10
T1MaxDelay = 60e-6
PulseType = 0  # 0 Gaussian, 1 Square
Avg = 50e3
T2RamseyDetuning = 0.3e6

# two tone setting
OneToneSpan = 40e6
TwoToneSpan = 20e6
CyclePoints = 400e3
TwoTonePower = -15
TwoToneAvg = 20e3
TwoToneSeqLen = 10e3

# pumped rabi setting
PumpFreq = 6.70635e9
PumpPower = -5
PumpLength = 61e-9
DriveDelayList = [1e-9]

# time domain measurements
# MeasTypeList = ['rabi']
# MeasTypeList = ['rabi CH1 drive']
# MeasTypeList = ['rabi CH1 pumped']
# MeasTypeList = ['t2_ramsey', 't1', 't2_echo']
# MeasTypeList = ['t1']
# MeasTypeList = ['t2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't1', 't2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't2_ramsey', 't2_echo']
# MeasTypeList = ['rabi', 't1', 't2_echo']
MeasTypeList = ['rabi', 't1']
# MeasTypeList = ['rabi', 't1_short']
# MeasTypeList = []
# MeasTypeList = ['rabi', 't1_t2_interleaved']
# MeasTypeList = ['rabi', 't2_ramsey', 't1_t2_interleaved']

if 'rabi CH1 pumped' not in MeasTypeList:
    DriveDelayList = [1e-9]
for i, cur in enumerate(CurrentList):
    if not FixReadoutFreq:
        [fc, width] = Script_TwoTone.FindReadoutFreqOneTone(cur, AnchorReadout1, AnchorReadout2, Span=OneToneSpan,
                                                            ReadoutPower=ReadoutPower, SaveFig=SaveFig)
        fc *= 1e9
        width *= 1e9
        AnchorReadout2 = AnchorReadout1
        AnchorReadout1 = [cur, fc]
        ReadoutFreq = fc - width / 2
    f0 = 1e9 * Script_TwoTone.FindQubitFreqTwoTone(cur, Anchor1, Anchor2, Power=TwoTonePower, ReadoutFreq=ReadoutFreq,
                                                   ReadoutPower=ReadoutPower, Avg=TwoToneAvg, Span=TwoToneSpan,
                                                   SeqLen=TwoToneSeqLen, SaveFig=SaveFig)
    # f0 = 0.5246e9
    Anchor2 = Anchor1
    Anchor1 = [cur, f0]
    print('Current = %.4GmA\nReadoutFreq = %.4GGHz\nQubitFreq = %.4GGHz' % (cur * 1e3, ReadoutFreq / 1e9, f0 / 1e9))
    for MeasType in MeasTypeList:
        for DriveDelay in DriveDelayList:
            FitDict = Script_TimeDomain.timeDomainMeasurement(cur, ReadoutFreq, f0, DrivingPower, PiPulseLength,
                                                              T1MaxDelay=T1MaxDelay, Avg=Avg, Detuning=T2RamseyDetuning,
                                                              DutyCyclePoints=CyclePoints, ReadoutPower=ReadoutPower,
                                                              PulseType=PulseType, MeasurementType=MeasType,
                                                              PumpFreq=PumpFreq, PumpPower=PumpPower,
                                                              PumpLength=PumpLength,
                                                              DriveDelay=DriveDelay, FitAndPlot=True, ShowFig=False)
        if MeasType == 'rabi':
            PiPulseLength = round(FitDict['opt'][3]) * 1e-9
            print('PiPulseLength now is %.3Gs' % PiPulseLength)
            # elif MeasType == 't1':
            #     T1MaxDelay = round(FitDict['opt'][1] * 5 / 1e3) * 1e-6
            #     print('T1MaxDelay now is %.3Gs' % T1MaxDelay)
