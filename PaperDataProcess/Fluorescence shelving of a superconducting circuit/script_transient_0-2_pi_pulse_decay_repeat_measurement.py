import MeasurementControlFunc as mcf

WaitHour = 0  # hr
Avg = 200e3
DutyCyclePoints = 3000e3

# Avg = 200
# DutyCyclePoints = 500e3
DataFolderName = '11112019_back to waveguide'

ItemDict = {
    'Alazar - Number of records': Avg,
    'Pulse Generator - Number of points': DutyCyclePoints,
}

mcf.Wait(int(WaitHour * 3600))
for i in range(10):
    MeasLabel = 'transient_P2_P1_interleaved'
    ConfigName = MeasLabel + '.hdf5'
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict, DataFolderName=DataFolderName)

    MeasLabel = 't1_P2_P1_interleaved'
    ConfigName = MeasLabel + '.hdf5'
    [OutPath, OutFile] = mcf.RunMeasurement(ConfigName, MeasLabel, ItemDict=ItemDict, DataFolderName=DataFolderName)