from QubitDataProcessPackages import*
from datetime import timedelta

FolderName = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619\\Shared log\\'
FileName = 'magnet_20190813.txt'
DataTable = np.genfromtxt(FolderName + FileName, skip_header=23, delimiter=',')
Time = DataTable[:, 0]
MagCur = DataTable[:, 2]
lines = []
File = open(FolderName + FileName)
for i in range(22):
    lines += [File.readline()]
DateLine = lines[10]
TimeLine = lines[11]
DateStr = DateLine[5:]
TimeStr = TimeLine[5:]

time_origin = [11, 0, 0]
time_origin_sec = timedelta(hours=time_origin[0], minutes=time_origin[1], seconds=time_origin[2]).total_seconds()

TimeLimit = (0, 46)

fig, ax = plt.subplots()
plt.plot(Time / 3600, MagCur)
# plt.title()
plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel('Current(A)', fontsize='x-large')
plt.xlim(TimeLimit)
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()