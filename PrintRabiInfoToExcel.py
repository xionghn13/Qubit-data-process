from openpyxl import Workbook
from QubitDataProcessPackages import *
from ReferencedTSweepPlot import plotReferencedTSweep
from QubitSpectrumFunc import printPercent


def checkBadFit(val, err):
    if err > abs(val) * 0.4:
        err = np.nan
        val = np.nan
    return [val, err]


def printRabiInfoFromFileList(FileList, OutputFolder, OutputFile, FixedFolder=None, FitDoubleExp=False,
                               LabberFolder='C:\\Users/admin\Labber\Data/'):
    wb = Workbook()
    ws = wb.active
    freqSL = []
    CurSL = []
    ASL = []
    AerrSL = []
    for ind, file in enumerate(FileList):
        file_folder = edf.getFolder(file)
        if not file.endswith('hdf5'):
            file += '.hdf5'
        FitDict = plotReferencedTSweep(file_folder, file, FitDoubleExponential=fit_double, ShowFig=False, SaveFig=False)
        if i == 0:
            freq = edf.readPumpFreqLabber(file_folder + file)
            cur = edf.readCurrentLabber(file_folder + file)
            freqSL += [freq]
            CurSL += [cur]
            A = FitDict['opt'][0]
            Aerr = np.sqrt(FitDict['cov'][0, 0])
            [A, Aerr] = checkBadFit(A, Aerr)
            ASL += [A]
            AerrSL += [Aerr]
        else:
            T2 = FitDict['opt'][1] / 1e3
            T2err = np.sqrt(FitDict['cov'][1, 1]) / 1e3
            if T2err > T2:
                T2 = np.nan
                T2err = np.nan
        T2SL += [T2]
        T2errSL += [T2err]
        printPercent(ind, FileList.__len__())

    table = [
        ['Freq'] + freqSL,
        ['Current'] + CurSL,
        ['A'] + ASL,
        ['AErr'] + AerrSL,
    ]
    for row in table:
        ws.append(row)
    wb.save(filename=OutputFolder + OutputFile)
    return


def getRabiNameList(NameFolder, NameFile):
    file = open(NameFolder + NameFile)
    names = file.readlines()
    FileList = []
    for n in names:
        if n.startswith('rabi'):
            if n.endswith('\n'):
                FileList += [n[:-1]]
            else:
                FileList += [n]
    return FileList

if __name__ == '__main__':
    NameFolder = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619/'
    NameFile = 'filename0720_2.txt'
    FileList = getRabiNameList(NameFolder, NameFile)
    # print(FileList)
    # FileList = [
    #     [
    #         't1_2019-06-25-12-34-18',
    #         't2_echo_2019-06-25-13-15-28',
    #     ],
    # ]
    FixedFolder = None
    LabberFolder = 'C:\\Users/admin\Labber\Data/'
    OutputFolder = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619/'
    OutputFileTag = 'wg5 in 8.5GHz cavity 0710 cd'

    for FitDoubleExp in [False, True]:
        OutputFile = OutputFileTag
        if FitDoubleExp:
            OutputFile += ' double exp'
        OutputFile += '.xlsx'
        printCoherenceFromFileList(FileList, OutputFolder, OutputFile, FitDoubleExp=FitDoubleExp)
