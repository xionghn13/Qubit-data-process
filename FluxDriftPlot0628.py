from QubitDataProcessPackages import *
from datetime import timedelta

# Haonan's data
DataPath = 'E:\Projects\Fluxonium\data_process\Fluxonium032619\\'
FileList = [
    'two tone_656.hdf5',
    'two tone_657.hdf5',
    'two tone_658.hdf5',
    'two tone_655.hdf5',
]

TimeStartStop = [
    {'start': [4, 35, 13], 'stop': [5, 18, 44]},
    {'start': [5, 35, 52], 'stop': [6, 12, 42]},
    {'start': [6, 15, 2], 'stop': [10, 26, 35]},
    {'start': [3, 5, 9], 'stop': [4, 32, 2]}
]

time_origin = [3, 5, 9]
time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()

Anchor1 = [7.7e-3, 658.51e-3]
Anchor2 = [7.78e-3, 544.1e-3]
slope = (Anchor1[1] - Anchor2[1]) / (Anchor1[0] - Anchor2[0])
hI = 2.57e-3  # A
Area = 459.7e-12 # m^2

FreqList = []
TimeList = []
ComplexList = []
fig, ax = plt.subplots()
for i, file in enumerate(FileList):
    [Freq, Count, Complex] = edf.readRepeatedFSweepTwoToneLabber(DataPath + file)
    TimeStamp = TimeStartStop[i]
    start_sec = timedelta(hours=TimeStamp['start'][0], minutes=TimeStamp['start'][1],
                          seconds=TimeStamp['start'][2]).total_seconds()
    stop_sec = timedelta(hours=TimeStamp['stop'][0], minutes=TimeStamp['stop'][1],
                         seconds=TimeStamp['stop'][2]).total_seconds()
    start_sec -= time_origin_sec
    stop_sec -= time_origin_sec
    FreqUniq = np.unique(Freq)
    if Freq[0, 0] > Freq[-1, -1]:
        FreqUniq = np.flip(FreqUniq, axis=0)
    CountUniq = np.unique(Count)
    if Count[0, 0] > Count[-1, -1]:
        CountUniq = np.flip(CountUniq, axis=0)
    Complex = (Complex / Complex.mean(axis=0))
    time = np.linspace(start_sec, stop_sec, len(CountUniq)) / 3600
    Current = (FreqUniq - Anchor1[1]) / slope
    Flux = Current / hI / 2 * Phi0
    B_field = Flux / Area * 1e8
    # B_field = Current / 50e-3 * 0.91e-4
    pcm = plt.pcolormesh(time, B_field, np.abs(Complex))
    FreqList += [FreqUniq]
    TimeList += [time]
    ComplexList += [Complex]

plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel("B($10^{-8} $ T)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()


# Long's data
FileList = [
    'Two_tone_twpa_on_38.hdf5',
]

TimeStartStop = [
    {'start': [6, 10, 0], 'stop': [11, 45, 0]},
]

fig, ax1 = plt.subplots(211)
for i, file in enumerate(FileList):
    [Freq, Count, Complex] = edf.readRepeatedFSweepTwoToneLabber(DataPath + file)
    TimeStamp = TimeStartStop[i]
    start_sec = timedelta(hours=TimeStamp['start'][0], minutes=TimeStamp['start'][1],
                          seconds=TimeStamp['start'][2]).total_seconds()
    stop_sec = timedelta(hours=TimeStamp['stop'][0], minutes=TimeStamp['stop'][1],
                         seconds=TimeStamp['stop'][2]).total_seconds()
    start_sec -= time_origin_sec
    stop_sec -= time_origin_sec
    FreqUniq = np.unique(Freq)
    if Freq[0, 0] > Freq[-1, -1]:
        FreqUniq = np.flip(FreqUniq, axis=0)
    CountUniq = np.unique(Count)
    if Count[0, 0] > Count[-1, -1]:
        CountUniq = np.flip(CountUniq, axis=0)
    Complex = (Complex / Complex.mean(axis=0))
    time = np.linspace(start_sec, stop_sec, len(CountUniq)) / 3600
    Current = (FreqUniq - Anchor1[1]) / slope
    Flux = Current / hI / 2 * Phi0
    B_field = Flux / Area * 1e8
    # B_field = Current / 50e-3 * 0.91e-4
    pcm = plt.pcolormesh(time, B_field, np.abs(Complex))
    FreqList += [FreqUniq]
    TimeList += [time]
    ComplexList += [Complex]

plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel("B($10^{-8} $ T)", fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
