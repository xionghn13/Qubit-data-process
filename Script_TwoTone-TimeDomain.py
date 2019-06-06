import Script_TwoTone
import Script_TimeDomain
from QubitDataProcessPackages import *

CurrentList = np.linspace(6.22e-3, 6.2e-3, 3)
Anchor1 = [6.23e-3, 572.2e6]
Anchor2 = [6.24e-3, 563.5e6]
SaveFig = True
Span = 10e6
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
for i, cur in enumerate(CurrentList):
    f0 = 1e9 * Script_TwoTone.FindQubitFreqTwoTone(cur, Anchor1, Anchor2, Span=Span, SaveFig=SaveFig)
    Anchor2 = Anchor1
    Anchor1 = [cur, f0]
    print('Current = %.3GmA\nQubitFreq = %.4GGHz' % (cur * 1e3, f0 / 1e9))
    for MeasType in MeasTypeList:
        FitDict = Script_TimeDomain.timeDomainMeasurement(cur, ReadoutFreq, f0, DrivingPower, PiPulseLength,
                                                          T1MaxDelay=T1MaxDelay,
                                                          PulseType=PulseType, MeasurementType=MeasType,
                                                          FitAndPlot=True, ShowFig=False)
        if MeasType == 'rabi':
            PiPulseLength = round(FitDict['opt'][3]) * 1e-9
            print('PiPulseLength now is %.3Gs' % PiPulseLength)
        elif MeasType == 't1':
            T1MaxDelay = round(FitDict['opt'][1] * 5 / 1e3) * 1e-6
            print('T1MaxDelay now is %.3Gs' % T1MaxDelay)
