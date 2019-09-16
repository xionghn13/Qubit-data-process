from QubitDataProcessPackages import *
import Single_small_junction as ssj
from xlrd import open_workbook
from CoherenceVSFreqPlot import getRowData, sortPlot

if __name__ == '__main__':
    I0 = -47.5e-3  # mA
    hI = 2.5975  # mA
    I_period = hI * 2  # mA

    Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
    TrialList = [
        # ['wg5 in 8.5GHz cavity at corner.xlsx'],

        [
            'wg5 in 8.5GHz cavity 0612 cd.xlsx',
            'wg5 in 8.5GHz cavity 0628 cd.xlsx',
            # 'wg5 in 8.5GHz cavity 0628 cd_2.xlsx',
            'wg5 in 8.5GHz cavity 0628 cd_3.xlsx'
        ],

        ['wg5 in 8.5GHz cavity 0710 cd.xlsx'],

        ['wg5 in 8.5GHz cavity 0803 cd.xlsx',
         'wg5 in 8.5GHz cavity 0803 cd 2.xlsx'],

        ['wg5 in 8.5GHz cavity 0830 cd.xlsx']
    ]
    ExtraPoints = [
        # np.array([[]]),
        np.array([[7.84, 514.9e-3, 81.7, 14, 166, 37, 31.1, 23, 0.499, 0.21]]),
        np.array([[]]),
        np.array([[2.545, 524e-3, 52.3, 7.9, 121, 37, 43.2, 32, 0.906, 0.36]]),
        np.array([[2.63, 524.6e-3, 32.3, 4.2, 78, 47, 52.4, 23, 1.76, 0.85]]),
    ]
    # cur, freq, T1, TR, Tqp, nqp
    LabelList = [
        # 'CD 5/31',
        'CD 6/12',
        'CD 7/10',
        'CD 8/3',
        'CD 8/30'
    ]

    xlim = (0.514, 0.53)

    ReadDataFromExcel = True

    FitDoubleExp = True

    FreqSingleList = []
    FreqRepeatedList = []
    T1SingleList = []
    T1RepeatedList = []
    T1ErrSingleList = []
    T1ErrRepeatedList = []
    TqpSingleList = []
    TqpRepeatedList = []
    TqpErrSingleList = []
    TqpErrRepeatedList = []
    nqpSingleList = []
    nqpRepeatedList = []
    nqpErrSingleList = []
    nqpErrRepeatedList = []
    FluxSingleList = []
    FluxRepeatedList = []
    FreqMergeList = []
    FluxMergeList = []

    T1MergeList = []
    T2MergeList = []
    TqpMergeList = []
    nqpMergeList = []

    for i, FileList in enumerate(TrialList):
        FreqSingle = np.array([])
        CurrentSingle = np.array([])
        T1Single = np.array([])
        T1ErrSingle = np.array([])
        TqpSingle = np.array([])
        TqpErrSingle = np.array([])
        nqpSingle = np.array([])
        nqpErrSingle = np.array([])
        T2Single = np.array([])
        T2ErrSingle = np.array([])
        if len(ExtraPoints[i][0]) == 0:
            FreqRepeated = np.array([])
            CurrentRepeated = np.array([])
            T1Repeated = np.array([])
            T1ErrRepeated = np.array([])
            TqpRepeated = np.array([])
            TqpErrRepeated = np.array([])
            nqpRepeated = np.array([])
            nqpErrRepeated = np.array([])
            T2Repeated = np.array([])
            T2ErrRepeated = np.array([])
        else:
            pts = ExtraPoints[i]
            FreqRepeated = pts[:, 1]
            CurrentRepeated = pts[:, 0]
            if FitDoubleExp:
                T1Repeated = pts[:, 4]
                T1ErrRepeated = pts[:, 5]
            else:
                T1Repeated = pts[:, 2]
                T1ErrRepeated = pts[:, 3]
            TqpRepeated = pts[:, 6]
            TqpErrRepeated = pts[:, 7]
            nqpRepeated = pts[:, 8]
            nqpErrRepeated = pts[:, 9]
            T2Repeated = np.array([])
            T2ErrRepeated = np.array([])

        for File in FileList:
            if FitDoubleExp:
                File = File[:-5] + ' double exp.xlsx'
            book = open_workbook(Folder + File, on_demand=True)
            sheet = book.sheet_by_index(0)
            FSrow = sheet.row(0)[1:]
            ThisFreqSingle = getRowData(sheet, 0)
            Len = len(ThisFreqSingle)

            FreqSingle = np.concatenate((FreqSingle, ThisFreqSingle))
            CurrentSingle = np.concatenate((CurrentSingle, getRowData(sheet, 1, Len)))
            T1Single = np.concatenate((T1Single, getRowData(sheet, 2, Len)))
            T1ErrSingle = np.concatenate((T1ErrSingle, getRowData(sheet, 3, Len)))
            if FitDoubleExp:
                TqpSingle = np.concatenate((TqpSingle, getRowData(sheet, 4, Len)))
                TqpErrSingle = np.concatenate((TqpErrSingle, getRowData(sheet, 5, Len)))
                nqpSingle = np.concatenate((nqpSingle, getRowData(sheet, 6, Len)))
                nqpErrSingle = np.concatenate((nqpErrSingle, getRowData(sheet, 7, Len)))
                T2Single = np.concatenate((T2Single, getRowData(sheet, 8, Len)))
                T2ErrSingle = np.concatenate((T2ErrSingle, getRowData(sheet, 9, Len)))

                ThisFreqRepeated = getRowData(sheet, 10, Len)
                Len = len(ThisFreqRepeated)
                FreqRepeated = np.concatenate((FreqRepeated, ThisFreqRepeated))
                CurrentRepeated = np.concatenate((CurrentRepeated, getRowData(sheet, 11, Len)))
                T1Repeated = np.concatenate((T1Repeated, getRowData(sheet, 12, Len)))
                T1ErrRepeated = np.concatenate((T1ErrRepeated, getRowData(sheet, 13, Len)))
                TqpRepeated = np.concatenate((TqpRepeated, getRowData(sheet, 14, Len)))
                TqpErrRepeated = np.concatenate((TqpErrRepeated, getRowData(sheet, 15, Len)))
                nqpRepeated = np.concatenate((nqpRepeated, getRowData(sheet, 16, Len)))
                nqpErrRepeated = np.concatenate((nqpErrRepeated, getRowData(sheet, 17, Len)))
                T2Repeated = np.concatenate((T2Repeated, getRowData(sheet, 18, Len)))
                T2ErrRepeated = np.concatenate((T2ErrRepeated, getRowData(sheet, 19, Len)))
            else:
                T2Single = np.concatenate((T2Single, getRowData(sheet, 4, Len)))
                T2ErrSingle = np.concatenate((T2ErrSingle, getRowData(sheet, 5, Len)))

                ThisFreqRepeated = getRowData(sheet, 6)
                Len = len(ThisFreqRepeated)
                FreqRepeated = np.concatenate((FreqRepeated, ThisFreqRepeated))
                CurrentRepeated = np.concatenate((CurrentRepeated, getRowData(sheet, 7, Len)))
                T1Repeated = np.concatenate((T1Repeated, getRowData(sheet, 8, Len)))
                T1ErrRepeated = np.concatenate((T1ErrRepeated, getRowData(sheet, 9, Len)))
                T2Repeated = np.concatenate((T2Repeated, getRowData(sheet, 10, Len)))
                T2ErrRepeated = np.concatenate((T2ErrRepeated, getRowData(sheet, 11, Len)))

        FluxSingle = (CurrentSingle - I0) / I_period
        FluxRepeated = (CurrentRepeated - I0) / I_period
        FreqMerge = np.concatenate((FreqRepeated, FreqSingle))
        FluxMerge = np.concatenate((FluxRepeated, FluxSingle))
        T1Merge = np.concatenate((T1Repeated, T1Single))
        T2Merge = np.concatenate((T2Repeated, T2Single))
        if FitDoubleExp:
            TqpMerge = np.concatenate((TqpRepeated, TqpSingle))
            nqpMerge = np.concatenate((nqpRepeated, nqpSingle))

        FreqSingleList += [FreqSingle]
        T1SingleList += [T1Single]
        T1ErrSingleList += [T1ErrSingle]
        TqpSingleList += [TqpSingle]
        TqpErrSingleList += [TqpErrSingle]
        nqpSingleList += [nqpSingle]
        nqpErrSingleList += [nqpErrSingle]

        FreqRepeatedList += [FreqRepeated]
        T1RepeatedList += [T1Repeated]
        T1ErrRepeatedList += [T1ErrRepeated]
        TqpRepeatedList += [TqpRepeated]
        TqpErrRepeatedList += [TqpErrRepeated]
        nqpRepeatedList += [nqpRepeated]
        nqpErrRepeatedList += [nqpErrRepeated]
        FluxSingleList += [FluxSingle]
        FluxRepeatedList += [FluxRepeated]
        FreqMergeList += [FreqMerge]
        FluxMergeList += [FluxMerge]

        T1MergeList += [T1Merge]
        T2MergeList += [T2Merge]
        if FitDoubleExp:
            TqpMergeList += [TqpMerge]
            nqpMergeList += [nqpMerge]

    if FitDoubleExp:
        T1Name = 'TR'
    else:
        T1Name = 'T1'
    leg = []
    # print(len(TrialList))
    # print(len(FreqSingleList))
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    for i in range(len(TrialList)):
        # print(i)
        FreqSingle = FreqSingleList[i]
        T1Single = T1SingleList[i]
        T1ErrSingle = T1ErrSingleList[i]
        TqpSingle = TqpSingleList[i]
        TqpErrSingle = TqpErrSingleList[i]
        nqpSingle = nqpSingleList[i]
        nqpErrSingle = nqpErrSingleList[i]
        T2Single = np.array([])
        T2ErrSingle = np.array([])

        FreqRepeated = FreqRepeatedList[i]
        T1Repeated = T1RepeatedList[i]
        T1ErrRepeated = T1ErrRepeatedList[i]
        TqpRepeated = TqpRepeatedList[i]
        TqpErrRepeated = TqpErrRepeatedList[i]
        nqpRepeated = nqpRepeatedList[i]
        nqpErrRepeated = nqpErrRepeatedList[i]
        T2Repeated = np.array([])
        T2ErrRepeated = np.array([])
        ExistRepeatedMeas = len(FreqRepeated) > 0 and np.any(FreqRepeated == FreqRepeated)
        color = next(ax._get_lines.prop_cycler)['color']
        # ax.errorbar(FreqSingle, T1Single, yerr=T1ErrSingle, fmt='.')
        plt.plot(FreqSingle, T1Single, '.', color=color)
        leg += [LabelList[i] + ': ' + T1Name + ' - Single']
        if ExistRepeatedMeas:
            # print(FreqRepeated)
            # print(T1Repeated)
            # print(T1ErrRepeated)
            # print(i)
            leg += [LabelList[i] + ': ' + T1Name + ' - Average']
            # ax.errorbar(FreqRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='o', color=color)
            plt.plot(FreqRepeated, T1Repeated, 'o', color=color, markersize=10)
        # ax.errorbar(FreqRepeated, T2Repeated, yerr=T2ErrRepeated, fmt='r^')
        # ax.errorbar(FreqSingle, T2Single, yerr=T2ErrSingle, fmt='b^')
        if FitDoubleExp:
            color = next(ax._get_lines.prop_cycler)['color']
            # ax.errorbar(FreqSingle, TqpSingle, yerr=TqpErrSingle, fmt='s')
            plt.plot(FreqSingle, TqpSingle, 's', color=color, markersize=5)
            leg += [LabelList[i] + ': ' + 'Tqp - Single']
            if ExistRepeatedMeas:
                leg += [LabelList[i] + ': ' + 'Tqp - Average']
                # ax.errorbar(FreqRepeated, TqpRepeated, yerr=TqpErrRepeated, fmt='s')
                plt.plot(FreqRepeated, TqpRepeated, 's', color=color, markersize=10)
            ax.set_yscale('log')
    plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))

    # for i in range(len(TrialList)):
    #     FreqMerge = FreqMergeList[i]
    #     T1Merge = T1MergeList[i]
    #     sortPlot(FreqMerge, T1Merge, ':')
    #     # sortPlot(FreqMerge, T2Merge, ':')
    #     if FitDoubleExp:
    #         TqpMerge = TqpMergeList[i]
    #         sortPlot(FreqMerge, TqpMerge, ':')

    plt.xlabel('Freq(GHz)', fontsize='x-large')
    plt.ylabel('Decay time(us)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.xlim(xlim)
    plt.ylim((4, 600))
    plt.tight_layout()
    # fig, ax = plt.subplots()
    # leg = []
    # ax.errorbar(FluxSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
    # if FitDoubleExp:
    #     ax.errorbar(FluxSingle, TqpSingle, yerr=TqpErrSingle, fmt='bs')
    #
    # ax.errorbar(FluxRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
    # if FitDoubleExp:
    #     ax.errorbar(FluxRepeated, TqpRepeated, yerr=TqpErrRepeated, fmt='rs')
    # sortPlot(FluxMerge, T1Merge, ':')
    # sortPlot(FluxMerge, T2Merge, ':')
    # if FitDoubleExp:
    #     sortPlot(FluxMerge, TqpMerge, ':')
    # plt.xlabel('Flux/Phi_0', fontsize='x-large')
    # plt.ylabel('Decay time(us)', fontsize='x-large')
    # plt.tick_params(axis='both', which='major', labelsize='x-large')
    # plt.tight_layout()
    # if FitDoubleExp:
    #     ax.set_yscale('log')
    # fig, ax = plt.subplots()
    # sortPlot(FluxMerge, FreqMerge, 'o:')
    # plt.xlabel('Flux/Phi_0', fontsize='x-large')
    # plt.ylabel('Freq(GHz)', fontsize='x-large')
    # plt.tick_params(axis='both', which='major', labelsize='x-large')
    # plt.tight_layout()

    if FitDoubleExp:
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        leg = []
        for i in range(len(TrialList)):
            FreqSingle = FreqSingleList[i]
            nqpSingle = nqpSingleList[i]
            nqpErrSingle = nqpErrSingleList[i]

            FreqRepeated = FreqRepeatedList[i]
            nqpRepeated = nqpRepeatedList[i]
            nqpErrRepeated = nqpErrRepeatedList[i]

            # ax.errorbar(FreqSingle, nqpSingle, yerr=nqpErrSingle, fmt='o')
            color = next(ax._get_lines.prop_cycler)['color']
            plt.plot(FreqSingle, nqpSingle, '.', color=color)
            leg += [LabelList[i] + ': ' + 'nqp - Single']
            if ExistRepeatedMeas:
                leg += [LabelList[i] + ': ' + 'nqp - Average']
                # ax.errorbar(FreqRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='o')
                plt.plot(FreqRepeated, nqpRepeated, 'o', color=color, markersize=10)
        plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
        # sortPlot(FreqMerge, nqpMerge, ':')
        plt.xlim(xlim)
        plt.xlabel('Freq(GHz)', fontsize='x-large')
        plt.ylabel('nqp', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        # fig, ax = plt.subplots()
        # ax.errorbar(FluxSingle, nqpSingle, yerr=nqpErrSingle, fmt='bo')
        # ax.errorbar(FluxRepeated, nqpRepeated, yerr=nqpErrRepeated, fmt='ro')
        # sortPlot(FluxMerge, nqpMerge, ':')
        # plt.xlabel('Flux/Phi_0', fontsize='x-large')
        # plt.ylabel('nqp', fontsize='x-large')
        # plt.tick_params(axis='both', which='major', labelsize='x-large')
        # plt.tight_layout()

    plt.show()
