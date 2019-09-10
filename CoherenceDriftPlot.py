from QubitDataProcessPackages import *
from datetime import timedelta
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot

DataPath = 'E:\Projects\Fluxonium\data_process\Fluxonium032619\\'
FileList = [
    't1_t2_interleaved_2019-08-12-11-55-08.hdf5',
]
Setup2FileList = [
    'T1_T2E_41.491_mA_2.hdf5'
]

FitDoubleExp = False
PlotSetup2Data = True

TimeStartStop = [
    {'start': [11, 55, 8], 'stop': [56, 17, 43]},
]
Setup2TimeStartStop = [
    {'start': [14, 41, 0], 'stop': [38, 3, 15]},
]
time_origin = [11, 55, 0]
TimeLimit = (2.5, 27.5)

time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()

TimeList = []
OptList = []
ErrList = []
# fig = plt.figure()
# ax = plt.subplot(211)
fig, ax = plt.subplots()
ax.grid(linestyle='--')
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
    if file.startswith('t1_t2_interleaved'):
        for ind in range(2):
            # print(Opt[ind][1, :] / 1000)
            # print(Err[ind][1, :] / 1000)
            ax.errorbar(time, Opt[ind][1, :] / 1000, yerr=Err[ind][1, :] / 1000, fmt=fmt_list[ind])
            if FitDoubleExp and ind == 0:
                ax.errorbar(time, Opt[ind][3, :] / 1000, yerr=Err[ind][3, :] / 1000, fmt='go')
    else:
        fmt_list = ['bo', 'ro']
        ax.errorbar(time, Opt[1, :] / 1000, yerr=Err[1, :] / 1000, fmt=fmt_list[0])
        if FitDoubleExp:
            ax.errorbar(time, Opt[3, :] / 1000, yerr=Err[3, :] / 1000, fmt='go')
    if i == 0:
        if FitDoubleExp:
            if file.startswith('t1_t2_interleaved'):
                plt.legend(['TR', 'Tqp', 'T2echo'])
            else:
                plt.legend(['TR', 'Tqp'])
        else:
            if file.startswith('t1_t2_interleaved'):
                plt.legend(['T1', 'T2echo'])
            else:
                plt.legend(['T1'])
    fmt_list = ['royalblue', 'violet']
    if file.startswith('t1_t2_interleaved'):
        for ind in range(2):
            plt.plot(time, Opt[ind][1, :] / 1000, '--', color=fmt_list[ind])
            if FitDoubleExp and ind == 0:
                plt.plot(time, Opt[ind][3, :] / 1000, 'g--')
    else:
        plt.plot(time, Opt[1, :] / 1000, '--', color=fmt_list[0])
        if FitDoubleExp:
            plt.plot(time, Opt[3, :] / 1000, 'g--')
    TimeList += [time]
    OptList += [Opt]
    ErrList += [Err]

plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel("Decay time(us)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.xlim(TimeLimit)
if FitDoubleExp:
    ax.set_yscale('log')
    plt.ylim([10, 1e3])

if PlotSetup2Data:
    Setup2TimeList = []
    Setup2OptList = []
    Setup2ErrList = []
    # fig = plt.figure()
    # ax = plt.subplot(211)
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    for i, file in enumerate(Setup2FileList):
        [Count, Opt, Err] = plotLabberRepeatedTSweepPlot(DataPath, file, Calibration=False, FitCorrectedR=True,
                                                         RotateComplex=True,
                                                         LogScale=False, FitDoubleExponential=FitDoubleExp,
                                                         PlotNumber=10, MinPlotInd=0, MaxPlotInd=50,
                                                         PlotIndex=[0, 3], T2MaxTime=3e5, ShowFig=False)
        # print(Err)
        TimeStamp = Setup2TimeStartStop[i]
        start_sec = timedelta(hours=TimeStamp['start'][0], minutes=TimeStamp['start'][1],
                              seconds=TimeStamp['start'][2]).total_seconds()
        stop_sec = timedelta(hours=TimeStamp['stop'][0], minutes=TimeStamp['stop'][1],
                             seconds=TimeStamp['stop'][2]).total_seconds()
        start_sec -= time_origin_sec
        stop_sec -= time_origin_sec
        time = np.linspace(start_sec, stop_sec, len(Count)) / 3600
        fmt_list = ['bo', 'ro']
        if file.startswith('T1_T2E'):
            for ind in range(2):
                # print(Opt[ind][1, :] / 1000)
                # print(Err[ind][1, :] / 1000)
                ax.errorbar(time, Opt[ind][1, :] / 1000, yerr=Err[ind][1, :] / 1000, fmt=fmt_list[ind])
                if FitDoubleExp and ind == 0:
                    ax.errorbar(time, Opt[ind][3, :] / 1000, yerr=Err[ind][3, :] / 1000, fmt='go')
        else:
            fmt_list = ['bo', 'ro']
            ax.errorbar(time, Opt[1, :] / 1000, yerr=Err[1, :] / 1000, fmt=fmt_list[0])
            if FitDoubleExp:
                ax.errorbar(time, Opt[3, :] / 1000, yerr=Err[3, :] / 1000, fmt='go')
        if i == 0:
            if FitDoubleExp:
                if file.startswith('T1_T2E'):
                    plt.legend(['TR', 'Tqp', 'T2echo'])
                else:
                    plt.legend(['TR', 'Tqp'])
            else:
                if file.startswith('T1_T2E'):
                    plt.legend(['T1', 'T2echo'])
                else:
                    plt.legend(['T1'])
        fmt_list = ['royalblue', 'violet']
        if file.startswith('T1_T2E'):
            for ind in range(2):
                plt.plot(time, Opt[ind][1, :] / 1000, '--', color=fmt_list[ind])
                if FitDoubleExp and ind == 0:
                    plt.plot(time, Opt[ind][3, :] / 1000, 'g--')
        else:
            plt.plot(time, Opt[1, :] / 1000, '--', color=fmt_list[0])
            if FitDoubleExp:
                plt.plot(time, Opt[3, :] / 1000, 'g--')
        Setup2TimeList += [time]
        Setup2OptList += [Opt]
        Setup2ErrList += [Err]

    plt.xlabel('Time(hr)', fontsize='x-large')
    plt.ylabel("Decay time(us)", fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    plt.xlim(TimeLimit)
    if FitDoubleExp:
        ax.set_yscale('log')
        plt.ylim([10, 1e3])

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    ax.errorbar(TimeList[0], OptList[0][0][1, :] / 1000, yerr=ErrList[0][0][1, :] / 1000, fmt='bo')
    ax.errorbar(Setup2TimeList[0], Setup2OptList[0][0][1, :] / 1000, yerr=Setup2ErrList[0][0][1, :] / 1000, fmt='ro')
    plt.legend(('Sample 1', 'Sample 2'))
    plt.plot(TimeList[0], OptList[0][0][1, :] / 1000, '--', color=fmt_list[0])
    plt.plot(Setup2TimeList[0], Setup2OptList[0][0][1, :] / 1000, '--', color=fmt_list[1])
    plt.xlabel('Time(hr)', fontsize='x-large')
    plt.ylabel("T1(us)", fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    plt.xlim(TimeLimit)
plt.show()