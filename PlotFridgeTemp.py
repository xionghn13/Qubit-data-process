from QubitDataProcessPackages import *
from datetime import timedelta

DateList = ['19-08-30', '19-08-31', '19-09-01']
lines = []
for date in DateList:
    FolderName = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619\\BF temp\\' + date + '\\'
    FileName = 'CH7 T ' + date + '.log'
    File = open(FolderName + FileName)
    lines += File.readlines()

time_origin = [8, 30, 19, 0, 0]
TimeLimit = (0, 36)
Num = len(lines)
Time = np.zeros((Num,))
Temp = np.zeros((Num,))
for i, line in enumerate(lines):
    line = line[1:-1]
    line_str = line.split(',')
    date_str = line_str[0]
    time_str = line_str[1]
    temp = float(line_str[2])
    if temp > 0.7e2 or temp < 0.00001:
        Temp[i] = np.nan
    else:
        Temp[i] = temp
    date_str_list = date_str.split('-')
    if i == 0:
        month = float(date_str_list[1])
        day_offset = 0
    else:
        if month < float(date_str_list[1]):
            day_offset = day
            month = float(date_str_list[1])
    day = float(date_str_list[0]) + day_offset
    time_str_list = time_str.split(':')
    hour = float(time_str_list[0])
    minute = float(time_str_list[1])
    second = float(time_str_list[2])
    Time[i] = timedelta(days=day, hours=hour, minutes=minute, seconds=second).total_seconds()
time_origin_sec = timedelta(days=time_origin[1], hours=time_origin[2], minutes=time_origin[3],
                            seconds=time_origin[4]).total_seconds()
Time = Time - time_origin_sec
TruncInd = (TimeLimit[0] < Time / 3600) * (Time / 3600 < TimeLimit[1])
fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(Time[TruncInd] / 3600, Temp[TruncInd] * 1e3)
# plt.title()
plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel('Temperature(mK)', fontsize='x-large')
plt.xlim(TimeLimit)
plt.ylim()
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()
