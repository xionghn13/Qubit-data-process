from QubitDataProcessPackages import *
from datetime import timedelta
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot

DataPath = 'E:\Projects\Fluxonium\data_process\Fluxonium032619\\'
FileList = [
    't1_t2_interleaved_2019-08-12-11-55-08.hdf5',
]

FitDoubleExp = True

TimeStartStop = [
    {'start': [11, 55, 8], 'stop': [56, 17, 43]},
]

time_origin = [11, 0, 0]
time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()

TimeLimit = (0, 46)

CountList = []
OptList = []
ErrList = []
# fig = plt.figure()
# ax = plt.subplot(211)
fig, ax = plt.subplots()
for i, file in enumerate(FileList):
    [Count, Opt, Err] = plotLabberRepeatedTSweepPlot(DataPath, file, Calibration=False, FitCorrectedR=True, RotateComplex=True,
                                 LogScale=False, FitDoubleExponential=FitDoubleExp, PlotNumber=10, MinPlotInd=0, MaxPlotInd=50,
                                 PlotIndex=[0, 3], T2MaxTime=3e5, ShowFig=False)
    # print(Err)
    TimeStamp = TimeStartStop[i]
    start_sec = timedelta(hours=TimeStamp['start'][0], minutes=TimeStamp['start'][1],
                          seconds=TimeStamp['start'][2]).total_seconds()
    stop_sec = timedelta(hours=TimeStamp['stop'][0], minutes=TimeStamp['stop'][1],
                         seconds=TimeStamp['stop'][2]).total_seconds()
    start_sec -= time_origin_sec
    stop_sec -= time_origin_sec
    time = np.linspace(start_sec, stop_sec, len(Count)) / 3600
    fmt_list = ['bo', 'ro']
    for ind in range(2):
        # print(Opt[ind][1, :] / 1000)
        # print(Err[ind][1, :] / 1000)
        ax.errorbar(time, Opt[ind][1, :] / 1000, yerr=Err[ind][1, :] / 1000, fmt=fmt_list[ind])
        if FitDoubleExp and ind == 0:
            ax.errorbar(time, Opt[ind][3, :] / 1000, yerr=Err[ind][3, :] / 1000, fmt=fmt_list[ind])
    if i == 0:
        if FitDoubleExp:
            plt.legend(['TR', 'Tqp', 'T2echo'])
        else:
            plt.legend(['T1', 'T2echo'])
    fmt_list = ['b--', 'r--']
    for ind in range(2):
        plt.plot(time, Opt[ind][1, :] / 1000, fmt_list[ind])
    CountList += [Count]
    OptList += [Opt]
    ErrList += [Err]

plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel("Decay time(us)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.xlim(TimeLimit)
if FitDoubleExp:
    ax.set_yscale('log')
plt.show()