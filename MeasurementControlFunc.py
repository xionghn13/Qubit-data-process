import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools


def RunMeasurement(ConfigName, MeasLabel, ItemDict={}):

    # set path to executable
    ScriptTools.setExePath('C:\Program Files (x86)\Labber\Program')

    TimeNow = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    TimeStr = str(TimeNow)
    Year = TimeStr[:4]
    Month = TimeStr[5:7]
    Day = TimeStr[8:10]
    # define measurement objects
    ConfigPath = 'E:\Data/fluxonium waveguide\labber config/'
    DataPath = 'C:/Users/admin\Labber\Data/' + Year + '/' + Month + '/' + 'Data_' + Month + Day + '/'


    if not os.path.exists(DataPath):
        os.makedirs(DataPath)


    ConfigFile = ConfigPath + ConfigName
    OutputFile = DataPath + MeasLabel + '_' + TimeStr
    print(MeasLabel + TimeStr)
    MeasObj = ScriptTools.MeasurementObject(ConfigFile, OutputFile)

    for item, value in ItemDict.items():
        if type(value) is list:
            MeasObj.updateValue(item, value[0], itemType=value[1])
        else:
            MeasObj.updateValue(item, value)

    # measure resonator
    MeasObj.performMeasurement(return_data=False)

