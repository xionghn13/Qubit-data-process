import os
import numpy as np
from datetime import datetime
from Labber import ScriptTools

ExperimentName = 'wg6 in 7.5GHz cavity'
CoolDownDate = 'test'

ConfigName = 'two tone sweep for cavity.hdf5'

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
OutputFile = DataPath + 'two tone' + '_' + TimeStr
MeasObj = ScriptTools.MeasurementObject(ConfigFile, OutputFile)

# MeasObj.updateValue('Alazar - Number of records', 3000)
# MeasObj.updateValue('Pump - Frequency', 291e6, itemType='START')
# MeasObj.updateValue('Pump - Frequency', 295e6, itemType='STOP')
# measure resonator
MeasObj.performMeasurement(return_data=False)