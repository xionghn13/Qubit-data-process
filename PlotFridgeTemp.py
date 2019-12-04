from QubitDataProcessPackages import *
from datetime import timedelta


def getTimeFromStrings(day, time_str):
    time_str_list = time_str.split(':')
    hour = float(time_str_list[0])
    minute = float(time_str_list[1])
    second = float(time_str_list[2])
    time_sec = timedelta(days=day, hours=hour, minutes=minute, seconds=second).total_seconds()
    return time_sec


def readFridgeTemperatureLog(flie_name_list):
    lines = []
    for file_name in flie_name_list:
        try:
            File = open(file_name)
            lines += File.readlines()
        except FileNotFoundError:
            print(file_name + ' not found')
    Num = len(lines)
    Time = np.zeros((Num,))
    Temp = np.zeros((Num,))
    for i, line in enumerate(lines):
        line = line[1:-1]
        line_str = line.split(',')
        date_str = line_str[0]
        time_str = line_str[1]
        temp = float(line_str[2])
        CorrectValue = True
        log_name = file_name.split('\\')[-1]
        if log_name.startswith('CH7'):
            if temp > 0.7e2 or temp < 7.11e-3:
                CorrectValue = False
        elif log_name.startswith('CH5'):
            if temp < 0.5:
                CorrectValue = False
        elif log_name.startswith('CH2'):
            if temp < 1.5:
                CorrectValue = False
        elif log_name.startswith('CH1'):
            if temp > 60:
                CorrectValue = False
        CorrectValue = True
        if not CorrectValue:
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
        Time[i] = getTimeFromStrings(day, time_str)
    return [Time, Temp]


def readFridgePressureLog(flie_name_list):
    lines = []
    for file_name in flie_name_list:
        File = open(file_name)
        lines += File.readlines()
    Num = len(lines)
    Time = np.zeros((Num,))
    Pressure = np.zeros((Num, 6))
    for i, line in enumerate(lines):
        line = line[0:-1]
        line_str = line.split(',')
        date_str = line_str[0]
        time_str = line_str[1]
        pres_ind_list = [5, 11, 17, 23, 29, 35]
        for j, ind in enumerate(pres_ind_list):
            pres = float(line_str[ind])
            if pres < 0.0000:
                Pressure[i, j] = np.nan
            else:
                Pressure[i, j] = pres
        date_str_list = date_str.split('-')
        if i == 0:
            month = float(date_str_list[1])
            day_offset = 0
        else:
            if month < float(date_str_list[1]):
                day_offset = day
                month = float(date_str_list[1])
        day = float(date_str_list[0]) + day_offset
        print(line)
        print(time_str)
        Time[i] = getTimeFromStrings(day, time_str)
    return [Time, Pressure]


def readFridgeFlowLog(flie_name_list):
    lines = []
    for file_name in flie_name_list:
        File = open(file_name)
        lines += File.readlines()
    Num = len(lines)
    Time = np.zeros((Num,))
    Flow = np.zeros((Num,))
    for i, line in enumerate(lines):
        line = line[0:-1]
        line_str = line.split(',')
        date_str = line_str[0]
        time_str = line_str[1]
        flow = float(line_str[2])
        Flow[i] = flow
        date_str_list = date_str.split('-')
        if i == 0:
            month = float(date_str_list[1])
            day_offset = 0
        else:
            if month < float(date_str_list[1]):
                day_offset = day
                month = float(date_str_list[1])
        day = float(date_str_list[0]) + day_offset
        Time[i] = getTimeFromStrings(day, time_str)
    return [Time, Flow]


DateList = ['19-10-30', '19-10-31', '19-11-01', '19-11-02', '19-11-03', '19-11-04']
TempChannelList = [1, 2, 5, 7]
AllChannelTempFileList = []
AllChannelPressureFileList = []
FlowFileList = []
for ch in TempChannelList:
    TempFileList = []
    for date in DateList:
        FolderName = 'Z:\Projects\BlueFors1\BlueFors Log\\' + date + '\\'
        TempFileList += [FolderName + 'CH' + str(ch) + ' T ' + date + '.log']
    AllChannelTempFileList += [TempFileList]

for date in DateList:
    FolderName = 'Z:\Projects\BlueFors1\BlueFors Log\\' + date + '\\'
    AllChannelPressureFileList += [FolderName + 'maxigauge ' + date + '.log']

for date in DateList:
    FolderName = 'Z:\Projects\BlueFors1\BlueFors Log\\' + date + '\\'
    FlowFileList += [FolderName + 'Flowmeter ' + date + '.log']

print(AllChannelTempFileList)

time_origin = [10, 30]
time_origin += [16, 42, 11]
TimeLimit = (-6.2, 130.3)

time_origin_sec = timedelta(days=time_origin[1], hours=time_origin[2], minutes=time_origin[3],
                                seconds=time_origin[4]).total_seconds()

[Time_P, Pressure] = readFridgePressureLog(AllChannelPressureFileList)
[Time_Flow, Flow] = readFridgeFlowLog(FlowFileList)
Time_P -= time_origin_sec
TruncInd_P = (TimeLimit[0] < Time_P / 3600) * (Time_P / 3600 < TimeLimit[1])
Time_Flow -= time_origin_sec
TruncInd_Flow = (TimeLimit[0] < Time_Flow / 3600) * (Time_Flow / 3600 < TimeLimit[1])
for i, ch in enumerate(TempChannelList):
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    [Time, Temp] = readFridgeTemperatureLog(AllChannelTempFileList[i])
    # Temp -= np.mean(Temp)
    Time = Time - time_origin_sec
    TruncInd = (TimeLimit[0] < Time / 3600) * (Time / 3600 < TimeLimit[1])
    plt.plot(Time[TruncInd] / 3600, Temp[TruncInd] * 1e3)
    leg += ['CH' + str(ch)]
    # plt.title()
    plt.xlabel('Time(hr)', fontsize='x-large')
    plt.ylabel('Temperature(mK)', fontsize='x-large')
    plt.legend(leg)
    plt.xlim(TimeLimit)
    plt.ylim()
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    # ax.set_yscale('log')
    plt.tight_layout()

for i in range(6):
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    plt.plot(Time_P[TruncInd_P] / 3600, Pressure[TruncInd_P, i])
    leg += ['CH' + str(i + 1)]
    # plt.title()
    plt.xlabel('Time(hr)', fontsize='x-large')
    plt.ylabel('Pressure(mbar)', fontsize='x-large')
    plt.legend(leg)
    plt.xlim(TimeLimit)
    plt.ylim()
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    # ax.set_yscale('log')
    plt.tight_layout()

fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(Time_Flow[TruncInd_Flow] / 3600, Flow[TruncInd_Flow])
plt.xlabel('Time(hr)', fontsize='x-large')
plt.ylabel('Flow(mmol/sec)', fontsize='x-large')
plt.xlim(TimeLimit)
plt.ylim()
plt.tick_params(axis='both', which='major', labelsize='x-large')
# ax.set_yscale('log')
plt.tight_layout()
plt.show()
