from openpyxl import Workbook
from QubitDataProcessPackages import *
from ReferencedTSweepPlot import plotReferencedTSweep
from QubitSpectrumFunc import printPercent


def checkBadFit(val, err):
    if err > abs(val) * 0.4:
        err = np.nan
        val = np.nan
    return [val, err]

def printCoherenceFromFileList(FileList, OutputFolder, OutputFile, FixedFolder=None, FitDoubleExp=False,
                               LabberFolder='C:\\Users/admin\Labber\Data/'):
    wb = Workbook()
    ws = wb.active
    freqSL = []
    CurSL = []
    T1SL = []
    T1errSL = []
    T2SL = []
    T2errSL = []
    TRSL = []
    TRerrSL = []
    TqpSL = []
    TqperrSL = []
    nqpSL = []
    nqperrSL = []
    for ind, t1t2 in enumerate(FileList):
        for i, file in enumerate(t1t2):
            if i == 0:
                fit_double = FitDoubleExp
            else:
                fit_double = False
            name_str_list = file.split('_')
            type_str = name_str_list[0]
            date_str_lsit = name_str_list[-1].split('-')
            year = date_str_lsit[0]
            month = date_str_lsit[1]
            day = date_str_lsit[2]
            file_folder = LabberFolder + year + '\\' + month + '\\' + 'Data_' + month + day + '\\'
            if not file.endswith('hdf5'):
                file += '.hdf5'
            FitDict = plotReferencedTSweep(file_folder, file, FitDoubleExponential=fit_double, ShowFig=False, SaveFig=False)
            if i == 0:
                freq = edf.readPumpFreqLabber(file_folder + file)
                cur = edf.readCurrentLabber(file_folder + file)
                freqSL += [freq]
                CurSL += [cur]
                if FitDoubleExp:
                    TR = FitDict['opt'][1] / 1e3
                    TRerr = np.sqrt(FitDict['cov'][1, 1]) / 1e3
                    [TR, TRerr] = checkBadFit(TR, TRerr)
                    Tqp = FitDict['opt'][3] / 1e3
                    Tqperr = np.sqrt(FitDict['cov'][3, 3]) / 1e3
                    [Tqp, Tqperr] = checkBadFit(Tqp, Tqperr)
                    nqp = FitDict['opt'][4]
                    nqperr = np.sqrt(FitDict['cov'][4, 4])
                    [nqp, nqperr] = checkBadFit(nqp, nqperr)
                    TRSL += [TR]
                    TRerrSL += [TRerr]
                    TqpSL += [Tqp]
                    TqperrSL += [Tqperr]
                    nqpSL += [nqp]
                    nqperrSL += [nqperr]
                else:
                    T1 = FitDict['opt'][1] / 1e3
                    T1err = np.sqrt(FitDict['cov'][1, 1]) / 1e3
                    [T1, T1err] = checkBadFit(T1, T1err)
                    T1SL += [T1]
                    T1errSL += [T1err]
                if t1t2.__len__() == 1:
                    T2 = np.nan
                    T2err = np.nan
            else:
                T2 = FitDict['opt'][1] / 1e3
                T2err = np.sqrt(FitDict['cov'][1, 1]) / 1e3
                if T2err > T2:
                    T2 = np.nan
                    T2err = np.nan
        T2SL += [T2]
        T2errSL += [T2err]
        printPercent(ind, FileList.__len__())
    if FitDoubleExp:
        table = [
            ['FreqSingle'] + freqSL,
            ['CurrentSingle'] + CurSL,
            ['TRSingle'] + TRSL,
            ['TRErrSingle'] + TRerrSL,
            ['TqpSingle'] + TqpSL,
            ['TqperrSingle'] + TqperrSL,
            ['nqpSingle'] + nqpSL,
            ['nqperrSingle'] + nqperrSL,
            ['T2Single'] + T2SL,
            ['T2ErrSingle'] + T2errSL,
            ['FreqRepeated'],
            ['CurrentRepeated'],
            ['TRRepeated'],
            ['TRErrRepeated'],
            ['TqpRepeated'],
            ['TqperrRepeated'],
            ['nqpRepeated'],
            ['nqperrRepeated'],
            ['T2Repeated'],
            ['T2ErrRepeated'],
        ]
    else:
        table = [
            ['FreqSingle'] + freqSL,
            ['CurrentSingle'] + CurSL,
            ['T1Single'] + T1SL,
            ['T1ErrSingle'] + T1errSL,
            ['T2Single'] + T2SL,
            ['T2ErrSingle'] + T2errSL,
            ['FreqRepeated'],
            ['CurrentRepeated'],
            ['T1Repeated'],
            ['T1ErrRepeated'],
            ['T2Repeated'],
            ['T2ErrRepeated'],
        ]
    for row in table:
        ws.append(row)

    wb.save(filename=OutputFolder + OutputFile)
    return


def getT1T2NameList(NameFolder, NameFile):
    file = open(NameFolder + NameFile)
    names = file.readlines()
    FileList = []
    for n in names:
        if n.startswith('t1'):
            FileList += [[n[:-1]]]
        elif n.startswith('t2'):
            FileList[-1] += [n[:-1]]
    return FileList

if __name__ == '__main__':
    NameFolder = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619/'
    NameFile = 'filename0705.txt'
    FileList = getT1T2NameList(NameFolder, NameFile)
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
    OutputFileTag = 'wg5 in 8.5GHz cavity 0628 cd_3'

    for FitDoubleExp in [False, True]:
        OutputFile = OutputFileTag
        if FitDoubleExp:
            OutputFile += ' double exp'
        OutputFile += '.xlsx'
        printCoherenceFromFileList(FileList, OutputFolder, OutputFile, FitDoubleExp=FitDoubleExp)
