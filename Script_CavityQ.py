import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools

# ExperimentName = 'wg5 in 8.5GHz cavity'
# CoolDownDate = 'test'

ConfigName = 'cavity Q.hdf5'

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
OutputFile = DataPath + ConfigName + TimeStr
MeasObj = ScriptTools.MeasurementObject(ConfigFile, OutputFile)

MeasObj.updateValue('Flux bias', value_1)
# measure resonator
MeasObj.performMeasurement(return_data=False)

