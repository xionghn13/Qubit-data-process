from QubitDataProcessPackages import*
from datetime import timedelta

FolderName = 'C:\SC Lab\Projects\Fluxonium\data_process\\flux\\'
FileName = 'loaner_ips_20201021_4pm.txt'

time_origin = [21, 16, 0, 0]

TimeLimit = (-0.5, 24 * 1)


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
date_str_list = DateStr.split('/')
day = float(date_str_list[2])
time_str_list = TimeStr.split(':')
hour = float(time_str_list[0])
minute = float(time_str_list[1])
second = float(time_str_list[2])
StartTime = [day, hour, minute, second]
start_time_sec = timedelta(days=StartTime[0], hours=StartTime[1], minutes=StartTime[2], seconds=StartTime[3]).total_seconds()
time_origin_sec = timedelta(days=time_origin[0], hours=time_origin[1], minutes=time_origin[2], seconds=time_origin[3]).total_seconds()
Time = Time + start_time_sec - time_origin_sec
TruncInd = (TimeLimit[0] < (Time / 3600)) * ((Time / 3600) < TimeLimit[1])
fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(Time[TruncInd] / 3600, MagCur[TruncInd])
# print(MagCur[TruncInd])
# print(Time[620000:620500] / 3600)
# print(Time[TruncInd]/3600)

# plt.title()
plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel('Current(A)', fontsize='x-large')
plt.xlim(TimeLimit)
plt.ylim()
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()