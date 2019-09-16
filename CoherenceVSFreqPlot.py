from QubitDataProcessPackages import *
import Single_small_junction as ssj
from xlrd import open_workbook


def getRowData(sheet, rowInd, length=None):
    Row = sheet.row(rowInd)[1:]
    RowData = np.array([])
    if length == None:
        for ele in Row:
            if ele.ctype == 0:
                break
            RowData = np.insert(RowData, len(RowData), ele.value)
    else:
        for i in range(length):
            ele = Row[i]
            if ele.ctype == 0:
                data = np.nan
            else:
                data = ele.value
            RowData = np.insert(RowData, len(RowData), data)
    return RowData

def sortPlot(x, y, spec):
    ind = x.argsort()
    x_sort = x[ind]
    y_sort = y[ind]
    plt.plot(x_sort, y_sort, spec)


if __name__ == '__main__':
    I0 = 32.5e-3  # mA
    hI = 2.5895  # mA
    I_period = hI * 2  # mA
    SortWithFlux = False
    PlotT1diel = False

    # for dielectric loss estimation
    N = 50
    EL = 0.437
    EC = 2.265
    EJ = 6.487
    loss_tan = 3.33e-6
    Q_cap = 1 / loss_tan
    T = 0e-3
    f01 = 0.4864

    Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
    File = 'wg5 in 8.5GHz cavity 0830 cd 3 double exp.xlsx'
    ReadDataFromExcel = True

    # I0 = 2.199
    # I_period = 2.469 * 2
    # FreqSingle = np.array(
    # [6.6199, 4.5596, 4.2657, 3.9635, 3.6598, 3.3561, 1.8472, 1.5516, 1.2614, 0.9822, 0.7298, 0.54, 0.4864])
    # CurrentSingle = np.array([2.5, 3.2, 3.3, 3.4, 3.5, 3.6, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.674])
    # T1Single = np.array([8.81, 19.3, 11.9, 29.6, 30.2, 27.3, 25.2, 23, 28, 41.6, 42.3, 48.1, 11.5])
    # T1ErrSingle = np.array([0.17, 0.49, 0.57, 0.66, 0.37, 0.5, 1.1, 1.9, 3.9, 4.9, 1.2, 1.2, 2.4])
    # T2Single = np.array([3.18, 4.85, 3.42, 3.8, 3.79, 3.75, 3.27, 4.11, 3.34, 4.29, 5.66, 5.33, 7.73])
    # T2ErrSingle = np.array([0.11, 0.42, 0.25, 0.26, 0.26, 0.35, 0.26, 0.48, 0.3, 0.54, 0.64, 0.27, 0.95])
    # FluxSingle = (CurrentSingle - I0) / I_period
    #
    # FreqRepeated = np.array(
    # [7.18915, 7.1005, 6.8807, 6.338, 6.0534, 5.7595, 5.462, 5.169, 4.87, 3.0522, 2.7488, 2.4463, 2.1462])
    # CurrentRepeated = np.array([2.22, 2.3, 2.4, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.7, 3.8, 3.9, 4])
    # T1Repeated = np.array([2.87, 3.74, 1.61, 7.37, 17.3, 21.2, 21.4, 17.1, 11.5, 8.98, 34.9, 17.1, 45.4])
    # T1ErrRepeated = np.array([0.22, 0.4, 0.068, 1.1, 1.3, 4.6, 0.78, 0.82, 0.68, 1, 2.5, 3.1, 5.3])
    # T2Repeated = np.array([4.34, 3.62, 2.05, 3.01, 3.77, 3.56, 3.85, 3.17, 2.97, 3.24, 3.6, 2.72, 3.46])
    # T2ErrRepeated = np.array([0.32, 0.41, 0.075, 0.76, 0.072, 0.17, 0.088, 0.11, 0.075, 0.4, 0.07, 0.32, 0.077])
    # FluxRepeated = (CurrentRepeated - I0) / I_period

    FitDoubleExp = File.endswith('double exp.xlsx')
    SpuriousMode = [0.51373, 0.51491, 0.51508, 0.51596, 0.51883, 0.52205, 0.52262]
    PlotSpuriousMode = False

    if ReadDataFromExcel:
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

    # FluxSortInd = FluxMerge.argsort()
    # FreqSortInd = FreqMerge.argsort()
    # FluxSort_FreqMerge = FreqMerge[FluxSortInd]
    # FreqSort_FreqMerge = FreqMerge[FreqSortInd]
    # FluxSort_FluxMerge = FluxMerge[FluxSortInd]
    # FluxSort_T1Merge = T1Merge[FluxSortInd]
    # FluxSort_T2Merge = T2Merge[FluxSortInd]
    # T1Merge = T1Merge[SortInd]
    # T2Merge = T2Merge[SortInd]
    # if FitDoubleExp:
    #     TqpMerge = TqpMerge[SortInd]
    #     nqpMerge = nqpMerge[SortInd]

    if PlotT1diel:
        MinFlux = np.min(FluxMerge)
        MaxFlux = np.max(FluxMerge)
        FluxPlot = np.linspace(MinFlux, MaxFlux, 21)
        T1diel = FluxPlot * 0
        for i, flux in enumerate(FluxPlot):
            [pem, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
            print(pem)
            # print('flux=%.3G, pem=%.3G, freq=%.3G' % (flux, pem, freq))
            T1diel[i] = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap, freq, pem, T)
            qsf.printPercent(i, 21)

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    ax.errorbar(FreqSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
    if not (T2Single.__len__() == 0 or np.all(T2Single != T2Single)):
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
    if PlotSpuriousMode:
        for mode in SpuriousMode:
            plt.axvline(x=mode, color='g', linestyle='--')

    plt.xlabel('Freq(GHz)', fontsize='x-large')
    plt.ylabel('Decay time(us)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    if PlotT1diel:
        sortPlot(FluxPlot, T1diel)
        leg = ['Dielectric loss tangent = %.3G' % loss_tan]
    ax.errorbar(FluxSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
    if not (T2Single.__len__() == 0 or np.all(T2Single != T2Single)):
        ax.errorbar(FluxSingle, T2Single, yerr=T2ErrSingle, fmt='b^')
    if FitDoubleExp:
        ax.errorbar(FluxSingle, TqpSingle, yerr=TqpErrSingle, fmt='bs')

    ax.errorbar(FluxRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
    ax.errorbar(FluxRepeated, T2Repeated, yerr=T2ErrRepeated, fmt='r^')
    if FitDoubleExp:
        ax.errorbar(FluxRepeated, TqpRepeated, yerr=TqpErrRepeated, fmt='rs')
        if FluxRepeated.__len__() == 0 or np.all(FluxRepeated != FluxRepeated):  # check if no data or array of nan
            legR = []
        else:
            legR = ['TR - Averge', 'T2echo - Averge', 'Tqp - Average']
        if T2Single.__len__() == 0 or np.all(T2Single != T2Single):
            plt.legend(['TR - Single', 'Tqp - Single'] + legR)
        else:
            plt.legend(leg + ['TR - Single', 'T2echo - Single', 'Tqp - Single'] + legR)
    else:
        if FluxRepeated.__len__() == 0:
            legR = []
        else:
            legR = ['TR - Averge', 'T2echo - Averge']
        if T2Single.__len__() == 0 or np.all(T2Single != T2Single):
            plt.legend(['T1 - Single'] + legR)
        else:
            plt.legend(['T1 - Single', 'T2echo - Single'] + legR)
    sortPlot(FluxMerge, T1Merge, ':')
    sortPlot(FluxMerge, T2Merge, ':')
    if FitDoubleExp:
        sortPlot(FluxMerge, TqpMerge, ':')
    plt.xlabel('Flux/Phi_0', fontsize='x-large')
    plt.ylabel('Decay time(us)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    if PlotT1diel or FitDoubleExp:
        ax.set_yscale('log')

    if FitDoubleExp:
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        ax.errorbar(FreqSingle, nqpSingle, yerr=nqpErrSingle, fmt='bo')
        ax.errorbar(FreqRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='ro')

        if FluxRepeated.__len__() == 0 or np.all(FluxRepeated != FluxRepeated):  # check if no data or array of nan
            legR = []
        else:
            legR = ['nqp - Average']
        plt.legend(['nqp - Single'] + legR)

        sortPlot(FreqMerge, nqpMerge, ':')
        plt.xlabel('Freq(GHz)', fontsize='x-large')
        plt.ylabel('nqp', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        ax.errorbar(FluxSingle, nqpSingle, yerr=nqpErrSingle, fmt='bo')
        ax.errorbar(FluxRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='ro')

        if FluxRepeated.__len__() == 0 or np.all(FluxRepeated != FluxRepeated):  # check if no data or array of nan
            legR = []
        else:
            legR = ['nqp - Average']
        plt.legend(['nqp - Single'] + legR)

        sortPlot(FluxMerge, nqpMerge, ':')
        if PlotSpuriousMode:
            for mode in SpuriousMode:
                plt.axvline(x=mode, color='g', linestyle='--')
        plt.xlabel('Flux/Phi_0', fontsize='x-large')
        plt.ylabel('nqp', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    sortPlot(FluxMerge, FreqMerge, 'o:')
    plt.xlabel('Flux/Phi_0', fontsize='x-large')
    plt.ylabel('Freq(GHz)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

    plt.show()
