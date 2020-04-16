from QubitDataProcessPackages import *
from datetime import timedelta
from LabberRepeatedTSweepPlot import plotLabberRepeatedTSweepPlot
import matplotlib.dates as mdates

# DataFolderName = '10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'
# DataFolderName = 'Data'
# DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2019/10\Data_1029\\'
# DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process\Fluxonium032619\\'
DataPath = 'C:\SC Lab\Projects\Fluxonium\data_process\ziggy4\\'
FileList = [
    # 't1_ramsey_echo_interleaved_1.hdf5',
    # 't1_ramsey_echo_interleaved_2.hdf5',
    # 't1_ramsey_echo_interleaved_3.hdf5',
    # 't1_ramsey_echo_interleaved_4.hdf5',
    # 't1_ramsey_echo_interleaved_5.hdf5',
    # 't1_ramsey_echo_interleaved_6.hdf5',
    # 't1_ramsey_echo_interleaved_7.hdf5',
    # 't1_ramsey_echo_interleaved_8.hdf5',
    # 't1_ramsey_echo_interleaved_9.hdf5',
    # 't1_ramsey_echo_interleaved_10.hdf5',
    # 't1_ramsey_echo_interleaved_11.hdf5',

    # 't1_t2_interleaved_2019-10-11-23-41-22_2.hdf5',
    # 't1_t2_interleaved_2019-10-12-14-11-02.hdf5',
    # 't1_t2_interleaved_2019-10-13-03-08-16.hdf5',
    # 't1_ramsey_echo_interleaved_2019-10-14-14-53-39_7.hdf5',
    # 't1_ramsey_echo_interleaved_2019-10-17.hdf5',
    # 't1_ramsey_echo_interleaved_2019-10-17_7.hdf5',
    # 't1_ramsey_echo_interleaved_2019-10-18_3.hdf5',

    't1_ramsey_echo_interleaved_15.hdf5',

]
Setup2DataPath = 'Z:\Projects\Transmon_Palmer\\2019\\10\Data_1017\\'
Setup2FileList = [
    'T1_4.hdf5',
    'T1_5.hdf5'
]

FitDoubleExp = False
PlotSetup2Data = False
NoErrorBar = False
SetTimeLimit = False
TimeLimit = (-1, 30)
MarkerSize = 1
LineSpec = '-'
# TimeStartStop = [
#     {'start': [21, 39, 30], 'stop': [24 + 15, 56, 1]},
#     # {'start': [17, 59, 21], 'stop': [17 + 22, 59 + 40, 21 + 3]},
#     # {'start': [19, 15, 25], 'stop': [19 + 15, 15 + 41, 25 + 40]},
#     # {'start': [11, 31, 23], 'stop': [48 + 12, 43, 0]},
#     # {'start': [12, 48, 17], 'stop': [12 + 24, 48 + 30, 17 + 10]},
#     # {'start': [13, 42, 50], 'stop': [13 + 22, 42 + 38, 50 + 20]},
#
# ]

# time_origin = [17, 59, 21]
# time_origin = TimeStartStop[0]['start']


# Setup2TimeStartStop = [
#     {'start': [14, 41, 0], 'stop': [38, 3, 15]},
# ]

# time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()
# time_origin_sec = calendar.timegm(time.strptime('2019-'))

TimeList = []
OptList = []
ErrList = []
# fig = plt.figure()
# ax = plt.subplot(211)
fig, ax = plt.subplots()
ax.grid(linestyle='--')
for i, file in enumerate(FileList):
    print(file)
    [Count, Opt, Err] = plotLabberRepeatedTSweepPlot(DataPath, file, Calibration=False, FitCorrectedR=True,
                                                     RotateComplex=True,
                                                     LogScale=False, FitDoubleExponential=FitDoubleExp, PlotNumber=10,
                                                     MinPlotInd=0, MaxPlotInd=50,
                                                     PlotIndex=[0, 3], T2MaxTime=3e6, ShowFig=False)
    if i == 0:
        [start_sec, stop_sec] = edf.readStartStopTime(DataPath + file)
        time_origin_sec = start_sec
        print('start time is ' + str(start_sec),
              '(s), which is ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_sec)))
    t_monitor = edf.getTimeStampArray(DataPath + file, len(Count)) - time_origin_sec
    t_monitor /= 3600
    t_monitor = edf.SecArrayToDateTimeList(edf.getTimeStampArray(DataPath + file, len(Count)))
    # print(t_monitor)
    # TimeStamp = TimeStartStop[i]
    # start_sec = timedelta(hours=TimeStamp['start'][0], minutes=TimeStamp['start'][1],
    #                       seconds=TimeStamp['start'][2]).total_seconds()
    # stop_sec = timedelta(hours=TimeStamp['stop'][0], minutes=TimeStamp['stop'][1],
    #                      seconds=TimeStamp['stop'][2]).total_seconds()
    # start_sec -= time_origin_sec
    # stop_sec -= time_origin_sec
    # t_monitor = np.linspace(start_sec, stop_sec, len(Count)) / 3600
    color_map_dic = {
        'T1': 0,
        'T2ramsey': 1,
        'T2echo': 2,
    }
    if file.startswith('t1_t2_interleaved'):
        meas_names = ['T1', 'T2echo']
    elif file.startswith('t1_ramsey_echo_interleaved'):
        meas_names = ['T1', 'T2ramsey', 'T2echo']
    fmt_list = ['bo', 'ro', 'ko', 'go']
    if file.startswith('t1_t2_interleaved') or file.startswith('t1_ramsey_echo_interleaved'):
        num_meas = len(Opt)
        for ind in range(num_meas):
            if ind != 1:
                if NoErrorBar:
                    plt.plot(t_monitor, Opt[ind][1, :] / 1000, fmt_list[ind], ms=MarkerSize)
                else:
                    ax.errorbar(t_monitor, Opt[ind][1, :] / 1000, yerr=Err[ind][1, :] / 1000,
                                fmt=fmt_list[color_map_dic[meas_names[ind]]], label=meas_names[ind])
                if FitDoubleExp and ind == 0:
                    ax.errorbar(t_monitor, Opt[ind][3, :] / 1000, yerr=Err[ind][3, :] / 1000, fmt='go')
    else:
        fmt_list = ['bo', 'ro']
        ax.errorbar(t_monitor, Opt[1, :] / 1000, yerr=Err[1, :] / 1000, fmt=fmt_list[0])
        if FitDoubleExp:
            ax.errorbar(t_monitor, Opt[3, :] / 1000, yerr=Err[3, :] / 1000, fmt='go')
    # if i == 0:
    #     if FitDoubleExp:
    #         if file.startswith('t1_t2_interleaved'):
    #             plt.legend(['TR', 'Tqp', 'T2echo'])
    #         else:
    #             plt.legend(['TR', 'Tqp'])
    #     else:
    #         if file.startswith('t1_t2_interleaved'):
    #             plt.legend(['T1', 'T2echo'])
    #         elif file.startswith('t1_ramsey_echo_interleaved'):
    #             plt.legend(['T1', 'T2ramsey', 'T2echo', 'T_phi'])
    #         else:
    #             plt.legend(['T1'])


    # fmt_list = ['royalblue', 'violet', 'grey']
    # if file.startswith('t1_t2_interleaved') or file.startswith('t1_ramsey_echo_interleaved'):
    #     for ind in range(num_meas):
    #         if i == len(FileList) - 1:
    #             plt.plot(t_monitor, Opt[ind][1, :] / 1000, LineSpec, label=meas_names[ind],
    #                      color=fmt_list[color_map_dic[meas_names[ind]]])
    #         else:
    #             plt.plot(t_monitor, Opt[ind][1, :] / 1000, LineSpec, color=fmt_list[color_map_dic[meas_names[ind]]])
    #         if FitDoubleExp and ind == 0:
    #             plt.plot(t_monitor, Opt[ind][3, :] / 1000, 'g' + LineSpec)
    # else:
    #     plt.plot(t_monitor, Opt[1, :] / 1000, LineSpec, color=fmt_list[0])
    #     if FitDoubleExp:
    #         plt.plot(t_monitor, Opt[3, :] / 1000, 'g' + LineSpec)

    TimeList += [t_monitor]
    OptList += [Opt]
    ErrList += [Err]
plt.legend()
plt.xlabel('Time', fontsize='x-large')
plt.ylabel("Decay time(us)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
if SetTimeLimit:
    plt.xlim(TimeLimit)
if FitDoubleExp:
    ax.set_yscale('log')
    plt.ylim([10, 1e3])
# ax.format_xdata = mdates.DateFormatter('%H-%M')
x_fmt = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(x_fmt)
plt.gcf().autofmt_xdate()
plt.gcf().set_size_inches(18, 4, forward=True)
plt.tight_layout()

if file.startswith('t1_t2_interleaved') or file.startswith('t1_ramsey_echo_interleaved'):
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    for i in range(len(FileList)):
        t_monitor = TimeList[i]
        Opt = OptList[i]
        fmt_list = ['bo', 'ro', 'ko', 'go']
        T_phi = 1 / (1 / Opt[-1][1, :] - 0.5 / Opt[0][1, :]) / 1000
        plt.plot(t_monitor, Opt[0][1, :] / 1000, fmt_list[0], ms=MarkerSize)
        plt.plot(t_monitor, T_phi, fmt_list[1], ms=MarkerSize)
        # plt.legend(['T1', 'T_phi'])
        fmt_list = ['royalblue', 'violet', 'grey']
        if i == len(FileList) - 1:
            plt.plot(t_monitor, Opt[0][1, :] / 1000, LineSpec, label='T1', color=fmt_list[0])
            plt.plot(t_monitor, T_phi, LineSpec, label='T_phi', color=fmt_list[1])
        else:
            plt.plot(t_monitor, Opt[0][1, :] / 1000, LineSpec, color=fmt_list[0])
            plt.plot(t_monitor, T_phi, LineSpec, color=fmt_list[1])
    plt.legend()
    plt.xlabel('Time', fontsize='x-large')
    plt.ylabel("Decay time(us)", fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    if SetTimeLimit:
        plt.xlim(TimeLimit)
    ax.set_yscale('log')
    plt.ylim([10, 1e3])
    plt.gcf().autofmt_xdate()
    plt.gcf().set_size_inches(15, 5, forward=True)
    plt.tight_layout()

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    T1 = np.array([])
    T2 = np.array([])
    for i in range(len(FileList)):
        t_monitor = TimeList[i]
        Opt = OptList[i]
        T1 = np.concatenate((T1, Opt[0][1, :] / 1e3))
        T2 = np.concatenate((T2, Opt[-1][1, :] / 1e3))
    print('cov mat = ')
    print(np.corrcoef(T1, T2))
    pcm = plt.hist2d(T1, T2, bins=20, cmap='PuBu', cmin=0.5)
    min_T = np.min(np.concatenate((T1, T2)))
    max_T = np.max(np.concatenate((T1, T2)))
    plt.xlabel('T1(us)', fontsize='x-large')
    plt.ylabel("T2echo(us)", fontsize='x-large')
    plt.xlim((min_T, max_T))
    plt.ylim((min_T, max_T))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Count', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

if PlotSetup2Data:
    Setup2TimeList = []
    Setup2OptList = []
    Setup2ErrList = []
    # fig = plt.figure()
    # ax = plt.subplot(211)
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    for i, file in enumerate(Setup2FileList):
        [Count, Opt, Err] = plotLabberRepeatedTSweepPlot(Setup2DataPath, file, Calibration=False, FitCorrectedR=True,
                                                         RotateComplex=True,
                                                         LogScale=False, FitDoubleExponential=FitDoubleExp,
                                                         PlotNumber=10, MinPlotInd=0, MaxPlotInd=50,
                                                         PlotIndex=[0, 3], T2MaxTime=3e5, ShowFig=False)
        t_monitor = edf.getTimeStampArray(Setup2DataPath + file, len(Count)) - time_origin_sec
        t_monitor /= 3600
        fmt_list = ['bo', 'ro']
        if file.startswith('T1_T2E'):
            for ind in range(2):
                if NoErrorBar:
                    plt.plot(t_monitor, Opt[ind][1, :] / 1000, fmt_list[ind], ms=MarkerSize)
                else:
                    ax.errorbar(t_monitor, Opt[ind][1, :] / 1000, yerr=Err[ind][1, :] / 1000, fmt=fmt_list[ind])
                if FitDoubleExp and ind == 0:
                    ax.errorbar(t_monitor, Opt[ind][3, :] / 1000, yerr=Err[ind][3, :] / 1000, fmt='go')
        else:
            fmt_list = ['bo', 'ro']
            if NoErrorBar:
                plt.plot(t_monitor, Opt[1, :] / 1000, fmt_list[0], ms=MarkerSize)
            else:
                ax.errorbar(t_monitor, Opt[1, :] / 1000, yerr=Err[1, :] / 1000, fmt=fmt_list[0])
            if FitDoubleExp:
                ax.errorbar(t_monitor, Opt[3, :] / 1000, yerr=Err[3, :] / 1000, fmt='go')
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
                plt.plot(t_monitor, Opt[ind][1, :] / 1000, LineSpec, color=fmt_list[ind])
                if FitDoubleExp and ind == 0:
                    plt.plot(t_monitor, Opt[ind][3, :] / 1000, 'g' + LineSpec)
        else:
            plt.plot(t_monitor, Opt[1, :] / 1000, LineSpec, color=fmt_list[0])
            if FitDoubleExp:
                plt.plot(t_monitor, Opt[3, :] / 1000, 'g' + LineSpec)
        Setup2TimeList += [t_monitor]
        Setup2OptList += [Opt]
        Setup2ErrList += [Err]

    plt.xlabel('Time', fontsize='x-large')
    plt.ylabel("Decay time(us)", fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    if SetTimeLimit:
        plt.xlim(TimeLimit)
    if FitDoubleExp:
        ax.set_yscale('log')
        plt.ylim([10, 1e3])

    # fig, ax = plt.subplots()
    # ax.grid(linestyle='--')
    # ax.errorbar(TimeList[0], OptList[0][0][1, :] / 1000, yerr=ErrList[0][0][1, :] / 1000, fmt='bo')
    # ax.errorbar(Setup2TimeList[0], Setup2OptList[0][0][1, :] / 1000, yerr=Setup2ErrList[0][0][1, :] / 1000, fmt='ro')
    # plt.legend(('Sample 1', 'Sample 2'))
    # plt.plot(TimeList[0], OptList[0][0][1, :] / 1000, LineSpec, color=fmt_list[0])
    # plt.plot(Setup2TimeList[0], Setup2OptList[0][0][1, :] / 1000, LineSpec, color=fmt_list[1])
    # plt.xlabel('Time', fontsize='x-large')
    # plt.ylabel("T1(us)", fontsize='x-large')
    # plt.tick_params(axis='both', which='major', labelsize='x-large')
    # plt.tight_layout()
    # if SetTimeLimit:
    #     plt.xlim(TimeLimit)

plt.show()
