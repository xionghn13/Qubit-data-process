from QubitDataProcessPackages import *
from datetime import timedelta


DataPath = 'E:\Projects\Fluxonium\data_process\Fluxonium032619\\'
FileList = [
    't1_t2_interleaved_2019-08-12-11-55-08.hdf5',
]

TimeStartStop = [
    {'start': [11, 55, 8], 'stop': [5, 18, 44]},
]

time_origin = [11, 55, 8]
time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()

TimeLimit = (0, 8.7)

CountList = []
OptList = []
ErrList = []
# fig = plt.figure()
# ax = plt.subplot(211)
fig, ax = plt.subplots()
for i, file in enumerate(FileList):
    [Count, Opt, Err] = edf.readRepeatedFSweepTwoToneLabber(DataPath + file)
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
        ax.errorbar(CounterArray, OptMatrixList[ind][1, :] / 1000, yerr=ErrMatrixList[ind][1, :] / 1000, fmt=fmt_list[ind])
    if i == 0:
        plt.legend(['T1', 'T2echo'])
    CountList += [Count]
    OptList += [Opt]
    ErrList += [Err]

plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel("Decay time(us)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.xlim(TimeLimit)