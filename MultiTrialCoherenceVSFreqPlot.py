from QubitDataProcessPackages import *
import Single_small_junction as ssj
from xlrd import open_workbook
from CoherenceVSFreqPlot import getRowData, sortPlot


if __name__ == '__main__':
    I0 = -47.5e-3  # mA
    hI = 2.5975  # mA
    I_period = hI * 2  # mA

    Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
    FileList = {
        'wg5 in 8.5GHz cavity 0612 cd.xlsx',
        'wg5 in 8.5GHz cavity 0628 cd.xlsx',
        'wg5 in 8.5GHz cavity 0628 cd_2.xlsx',
        'wg5 in 8.5GHz cavity 0628 cd_3.xlsx',
    }
    ReadDataFromExcel = True

    FitDoubleExp = False

    for File in FileList:
        book = open_workbook(Folder + File, on_demand=True)
        sheet = book.sheet_by_index(0)
        FSrow = sheet.row(0)[1:]
        FreqSingle = getRowData(sheet, 0)
        Len = len(FreqSingle)

        CurrentSingle = getRowData(sheet, 1, Len)
        T1Single = getRowData(sheet, 2, Len)
        T1ErrSingle = getRowData(sheet, 3, Len)
        if FitDoubleExp:
            TqpSingle = getRowData(sheet, 4, Len)
            TqpErrSingle = getRowData(sheet, 5, Len)
            nqpSingle = getRowData(sheet, 6, Len)
            nqpErrSingle = getRowData(sheet, 7, Len)
            T2Single = getRowData(sheet, 8, Len)
            T2ErrSingle = getRowData(sheet, 9, Len)

            FreqRepeated = getRowData(sheet, 10, Len)
            Len = len(FreqRepeated)
            CurrentRepeated = getRowData(sheet, 11, Len)
            T1Repeated = getRowData(sheet, 12, Len)
            T1ErrRepeated = getRowData(sheet, 13, Len)
            TqpRepeated = getRowData(sheet, 14, Len)
            TqpErrRepeated = getRowData(sheet, 15, Len)
            nqpRepeated = getRowData(sheet, 16, Len)
            nqpErrRepeated = getRowData(sheet, 17, Len)
            T2Repeated = getRowData(sheet, 18, Len)
            T2ErrRepeated = getRowData(sheet, 19, Len)
        else:
            T2Single = getRowData(sheet, 4, Len)
            T2ErrSingle = getRowData(sheet, 5, Len)

            FreqRepeated = getRowData(sheet, 6)
            Len = len(FreqRepeated)
            CurrentRepeated = getRowData(sheet, 7, Len)
            T1Repeated = getRowData(sheet, 8, Len)
            T1ErrRepeated = getRowData(sheet, 9, Len)
            T2Repeated = getRowData(sheet, 10, Len)
            T2ErrRepeated = getRowData(sheet, 11, Len)

        FluxSingle = (CurrentSingle - I0) / I_period
        FluxRepeated = (CurrentRepeated - I0) / I_period

        FreqMerge = np.concatenate((FreqRepeated, FreqSingle))
        FluxMerge = np.concatenate((FluxRepeated, FluxSingle))
        T1Merge = np.concatenate((T1Repeated, T1Single))
        T2Merge = np.concatenate((T2Repeated, T2Single))
        if FitDoubleExp:
            TqpMerge = np.concatenate((TqpRepeated, TqpSingle))
            nqpMerge = np.concatenate((nqpRepeated, nqpSingle))

        fig, ax = plt.subplots()
        ax.errorbar(FreqSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
        ax.errorbar(FreqSingle, T2Single, yerr=T2ErrSingle, fmt='b^')
        if FitDoubleExp:
            ax.errorbar(FreqSingle, TqpSingle, yerr=TqpErrSingle, fmt='bs')

        ax.errorbar(FreqRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
        ax.errorbar(FreqRepeated, T2Repeated, yerr=T2ErrRepeated, fmt='r^')

        if FitDoubleExp:
            ax.errorbar(FreqRepeated, TqpRepeated, yerr=TqpErrRepeated, fmt='rs')
            if FreqRepeated.__len__() == 0 or np.all(FluxRepeated != FluxRepeated):
                legR = []
            else:
                legR = ['TR - Averge', 'T2echo - Averge', 'Tqp - Average']
            if T2Single.__len__() == 0 or np.all(T2Single != T2Single):
                plt.legend(['TR - Single', 'Tqp - Single'] + legR)
            else:
                plt.legend(['TR - Single', 'T2echo - Single', 'Tqp - Single'] + legR)
            ax.set_yscale('log')

        else:
            if FreqRepeated.__len__() == 0:
                legR = []
            else:
                legR = ['TR - Averge', 'T2echo - Averge']
            if T2Single.__len__() == 0 or np.all(T2Single != T2Single):
                plt.legend(['T1 - Single'] + legR)
            else:
                plt.legend(['T1 - Single', 'T2echo - Single'] + legR)

        sortPlot(FreqMerge, T1Merge, ':')
        sortPlot(FreqMerge, T2Merge, ':')
        if FitDoubleExp:
            sortPlot(FreqMerge, TqpMerge, ':')

        plt.xlabel('Freq(GHz)', fontsize='x-large')
        plt.ylabel('Decay time(us)', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

        fig, ax = plt.subplots()
        leg = []
        ax.errorbar(FluxSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
        if FitDoubleExp:
            ax.errorbar(FluxSingle, TqpSingle, yerr=TqpErrSingle, fmt='bs')

        ax.errorbar(FluxRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
        if FitDoubleExp:
            ax.errorbar(FluxRepeated, TqpRepeated, yerr=TqpErrRepeated, fmt='rs')
        sortPlot(FluxMerge, T1Merge, ':')
        sortPlot(FluxMerge, T2Merge, ':')
        if FitDoubleExp:
            sortPlot(FluxMerge, TqpMerge, ':')
        plt.xlabel('Flux/Phi_0', fontsize='x-large')
        plt.ylabel('Decay time(us)', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if FitDoubleExp:
            ax.set_yscale('log')

        if FitDoubleExp:
            fig, ax = plt.subplots()
            ax.errorbar(FreqSingle, nqpSingle, yerr=nqpErrSingle, fmt='bo')
            ax.errorbar(FreqRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='ro')


            sortPlot(FreqMerge, nqpMerge, ':')
            plt.xlabel('Freq(GHz)', fontsize='x-large')
            plt.ylabel('nqp', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

            fig, ax = plt.subplots()
            ax.errorbar(FluxSingle, nqpSingle, yerr=nqpErrSingle, fmt='bo')
            ax.errorbar(FluxRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='ro')

            if FluxRepeated.__len__() == 0 or np.all(FluxRepeated != FluxRepeated):  # check if no data or array of nan
                legR = []
            else:
                legR = ['nqp - Average']
            plt.legend(['nqp - Single'] + legR)

            sortPlot(FluxMerge, nqpMerge, ':')
            plt.xlabel('Flux/Phi_0', fontsize='x-large')
            plt.ylabel('nqp', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        fig, ax = plt.subplots()
        sortPlot(FluxMerge, FreqMerge, 'o:')
        plt.xlabel('Flux/Phi_0', fontsize='x-large')
        plt.ylabel('Freq(GHz)', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

    plt.show()
