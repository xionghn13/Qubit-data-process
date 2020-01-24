from QubitDataProcessPackages import *
from QubitDecayFunc import T1_curve, rabi_curve, AutoRotate, DoubleExp_curve, TwoExp_curve
import os


def plotMultiPopulationTSweep(DataPath, RabiFile, BackgroundFolder='', BackgroundFile='', CircleCorrection=False,
                              CorrectionParam=[1, 0.027, 0.798, -0.0098], PopulationConversionConst=[1, 1.390],
                              IQModFreq=0.05, PhaseSlope=326.7041108065019, PhaseReferenceFreq=4.105, Calibration=False,
                              FitCorrectedR=False, LimitTimeRange=False, RotateComplex=True, FitDoubleExponential=False,
                              FitTwoExponential=True, StartTime=5000, EndTime=1e8, SaveFig=True, ShowFig=False,
                              LogScale=False):
    # read calibration file
    if Calibration:
        if BackgroundFile.startswith('transient'):
            [Time, ComplexLowerFreq] = edf.readRabiLabber(BackgroundFolder + BackgroundFile)
            BackPower = edf.readReadoutPowerLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            BackFreq = edf.readReadoutFreqLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackFreq = np.array([BackFreq])
            Plus50MHzBackComplex = np.array([ComplexLowerFreq.mean()])
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex
        elif BackgroundFile.startswith('calibration') or BackgroundFile.startswith('power'):
            [BackFreq, BackComplex] = edf.readFSweepLabber(BackgroundFolder + BackgroundFile)
            BackPower = edf.readReadoutPowerLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            Plus50MHzBackFreq = BackFreq
            Plus50MHzBackComplex = BackComplex
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex
        elif BackgroundFile.startswith('RefRabiCal'):
            # use rabi to calibrate rabi. Plus or minus are merged in one array.
            [BackFreq, BackPower, BackComplex] = edf.readRefRabiCalLabber(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            Plus50MHzBackFreq = BackFreq
            Plus50MHzBackComplex = BackComplex
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex

    # read data file
    ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
    ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
    if RabiFile.startswith('transient'):
        MeasurementType = 'transient no ref'
        [Time, Complex] = edf.readMultiRabiLabber(DataPath + RabiFile)
    elif RabiFile.startswith('t1'):
        MeasurementType = 't1 no ref'
        [Time, Complex] = edf.readMultiT1Labber(DataPath + RabiFile)

    if LimitTimeRange:
        TimeInd = (EndTime >= Time) == (Time >= StartTime)
        Time = Time[TimeInd]
        if len(Complex.shape) == 1:
            Complex = Complex[TimeInd]
        else:
            Complex = Complex[TimeInd, :]

    if MeasurementType in ('transient no ref', 't1 no ref'):
        ComplexNormalized = Complex * 10 ** (- ReadoutPower / 20)
        if Calibration:
            RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, Plus50MHzBackFreq,
                                                      Plus50MHzBackComplex, Plus50MHzBackPower)
        else:
            RComplex = ComplexNormalized

        if RotateComplex:
            RComplex = AutoRotate(RComplex)

        [amp_cor_re_fit, amp_cor_im_fit, P0_fit, P0_im_fit] = CorrectionParam
        if CircleCorrection:
            RComplex /= (amp_cor_re_fit + amp_cor_im_fit * 1j)
            RComplex = (RComplex - 1) / (P0_fit + P0_im_fit * 1j) * np.abs(
                P0_fit + P0_im_fit * 1j) + 1

    x_data = np.array(Time, dtype='float64')
    num_curve = RComplex.shape[1]
    for i in range(num_curve):
        y_data = np.array(RComplex[:, i].real, dtype='float64')
        # y_data = np.sum(np.array(RComplex[:, [0, 1, 3]].real, dtype='float64'), axis=1)
        if MeasurementType in ('transient no ref', 't1 no ref'):
            B_guess = y_data[-1]
            A_guess = y_data.max() - y_data.min()
            T1_guess = x_data[-1] / 5
            # bounds = (
            #     (-2, 1, -1),
            #     (2, 1e6, 1)
            # )
            bounds = (
                (- np.inf, 1, - np.inf),
                (np.inf, 1e6, np.inf)
            )

            if FitDoubleExponential:
                try:
                    opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                                         p0=[A_guess, T1_guess, B_guess, T1_guess * 0.1, 0.5],
                                         maxfev=300000)
                    # print('guess = %s' % str([A_guess, T1_guess, B_guess, T1_guess * 0.1, 1]))
                    # print('Double exp fit opt = %s' % str(opt))
                except RuntimeError:
                    print("Error - curve_fit failed")
                    opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
                    cov = np.zeros([len(opt), len(opt)])

                A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
                A_std, TR_std, B_std, Tqp_std, lamb_std = np.sqrt(cov.diagonal())

                TimeFit = np.linspace(Time.min(), Time.max(), 200)

                FitR = DoubleExp_curve(TimeFit, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                y_pred = DoubleExp_curve(Time, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
            elif FitTwoExponential:
                bounds = (
                    (- 1, 0.1, - 1, 0.01, -1),
                    (1, 1e6, 1, 1e6, 1)
                )
                try:
                    opt, cov = curve_fit(TwoExp_curve, x_data, y_data,
                                         p0=[A_guess, T1_guess, -A_guess, T1_guess * 0.1, B_guess], bounds=bounds,
                                         maxfev=300000)
                    print('guess = %s' % str([A_guess, T1_guess, -A_guess, T1_guess * 0.1, B_guess]))
                    # print('Double exp fit opt = %s' % str(opt))
                except (RuntimeError, TypeError, ValueError) as e:
                    print("Error - curve_fit failed")
                    opt = np.array([A_guess / 2, T1_guess, A_guess / 2, T1_guess * 0.1, B_guess])
                    cov = np.zeros([len(opt), len(opt)])

                A_fit, T1_fit, B_fit, T2_fit, C_fit = opt
                A_std, T1_std, B_std, T2_std, C_std = np.sqrt(cov.diagonal())

                TimeFit = np.linspace(Time.min(), Time.max(), 200)

                FitR = TwoExp_curve(TimeFit, A_fit, T1_fit, B_fit, T2_fit, C_fit)
                y_pred = TwoExp_curve(Time, A_fit, T1_fit, B_fit, T2_fit, C_fit)
                ParamList = ['A', 'T_Exp1(ns)', 'B', 'T_Exp2(ns)', 'C']
            else:
                try:
                    opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
                except (RuntimeError, TypeError, ValueError) as e:
                    print("Error - curve_fit failed")
                    opt = np.array([A_guess, T1_guess, B_guess])
                    cov = np.zeros([len(opt), len(opt)])
                A_fit, T1_fit, B_fit = opt
                A_std, T1_std, B_std = np.sqrt(cov.diagonal())
                TimeFit = np.linspace(Time.min(), Time.max(), 200)

                FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
                y_pred = T1_curve(Time, A_fit, T1_fit, B_fit)
                y_guess = T1_curve(Time, A_guess, T1_guess, B_guess)
                ParamList = ['A', 'Decay time/ns', 'B']
        # concatenate data
        if i == 0:
            TimeList = []
            TimeFitList = []
            RComplexList = []
            FitRList = []
            OptMatrix = np.reshape(opt, (len(opt), 1))
            ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
        else:
            OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
            ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)

        TimeList.append(Time)
        TimeFitList.append(TimeFit)
        # RComplexList.append(RComplex[:, i])
        RComplexList.append(y_data)
        FitRList.append(FitR)
        for i_p, par in enumerate(opt):
            if np.abs(ErrMatrix[i_p, -1]) > 0.4 * np.abs(par) or opt[1] > 5 * Time[-1]:
                # OptMatrix[i_p, -1] = np.nan
                OptMatrix[i_p, -1] = par
            else:
                OptMatrix[i_p, -1] = par

    # print(cov)
    limit = 1.7
    print(np.real(RComplex[-1, 3]))
    if ShowFig:
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        for i in range(num_curve):
            plt.plot(np.real(RComplex[:, i]), np.imag(RComplex[:, i]))
        if Calibration:
            plt.plot([-2, 2], [0, 0], '--')
            plt.plot([1], [0], 'ro')
        plt.xlabel('Re', fontsize='x-large')
        plt.ylabel('Im', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if Calibration:
            plt.xlim(-limit, limit)
            plt.ylim(-limit, limit)
        ax.set_aspect('equal')

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        for i in range(num_curve):
            plt.plot(Time, RComplex.real[:, i], 'o')
        plt.xlabel('Time(ns)', fontsize='x-large')
        plt.ylabel('Re', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        #
        # fig, ax = plt.subplots()
        # ax.grid(linestyle='--')
        # for i in range(num_curve):
        #     plt.plot(Time, (np.angle(ComplexNormalized[:, i]) - PhaseSlope * (ReadoutFreq - 4.105)) % (2 * np.pi), 'o')
        # plt.xlabel('Time/ns', fontsize='x-large')
        # plt.ylabel('Phase', fontsize='x-large')
        # plt.tick_params(axis='both', which='major', labelsize='x-large')
        # plt.tight_layout()

        Population = (PopulationConversionConst[0] - RComplex.real) * PopulationConversionConst[1]
        CorrectP2 = True
        PlotErrBar = False
        if CorrectP2:
            P2PiPulse = 108
            # P2RabiT1 = 604
            P2RabiT1 = 990
            k = np.exp(- P2PiPulse / P2RabiT1)
            P0 = np.mean(Population[:, [0, 2]], axis=1)
            P2 = Population[:, 3]
            print(P2[0])
            Population[:, 3] = 2 / (k + 1) * (P2 + (k - 1) / 2 * P0)
            print(P2[0])
            # print(Population[:, 3][0])

        if num_curve == 4:
            P0 = np.mean(Population[:, [0, 2]], axis=1)
            P1 = Population[:, 1]
            P2 = Population[:, 3]
            P_sum = P0 + P1 + P2
            P01 = P0 + P1
            print('std:', [P0.std(), P1.std(), P2.std(), P_sum.std(), P01.std()])
            print('avg:', [P0.mean(), P1.mean(), P2.mean(), P_sum.mean(), P01.mean()])
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            if PlotErrBar:
                stds = [0.003971263619235736, 0.007063671960122041, 0.004772838906966845, 0.009408918991044542,
                        0.007791820522548097]
                ax.errorbar(Time, P0, yerr=np.ones_like(P0) * stds[0], fmt='o', label='P0')
                ax.errorbar(Time, P1, yerr=np.ones_like(P0) * stds[1], fmt='o', label='P1')
                ax.errorbar(Time, P2, yerr=np.ones_like(P0) * stds[2], fmt='o', label='P2')
                ax.errorbar(Time, P_sum, yerr=np.ones_like(P0) * stds[3], fmt='o', label='P_sum')
                ax.errorbar(Time, P01, yerr=np.ones_like(P0) * stds[4], fmt='o', label='P0+P1')
            else:
                plt.plot(Time, P0, 'o', label='P0')
                plt.plot(Time, P1, 'o', label='P1')
                plt.plot(Time, P2, 'o', label='P2')
                plt.plot(Time, P_sum, 'o', label='P_sum')
                plt.plot(Time, P01, 'o', label='P0+P1')

            plt.legend()
            plt.xlabel('Time(ns)', fontsize='x-large')
            plt.ylabel('Population', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
        else:
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plt.plot(Time, Population, 'o')
            plt.xlabel('Time(ns)', fontsize='x-large')
            plt.ylabel('Population', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        plot_level = 'P0'

        if plot_level == 'P0':
            plot_ind = [0, 2]
        elif plot_level == 'P1':
            plot_ind = [1]
        elif plot_level == 'P2':
            plot_ind = [3]
        elif plot_level == 'P_sum':
            plot_ind = [0, 1, 3]
        P_plot = (PopulationConversionConst[0] - np.mean(np.array(RComplexList).transpose()[:, plot_ind].real,
                                                         axis=1)) * PopulationConversionConst[1]
        fit_plot = (PopulationConversionConst[0] - np.mean(np.array(FitRList).transpose()[:, plot_ind], axis=1)) * \
                   PopulationConversionConst[1]

        def r_to_P(x):
            return (PopulationConversionConst[0] - x) * PopulationConversionConst[1]

        opt = np.mean(OptMatrix[:, plot_ind], axis=1)
        err = ErrMatrix[:, plot_ind[0]]
        if FitDoubleExponential:
            A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
            A_std, TR_std, B_std, Tqp_std, lamb_std = err
        elif FitTwoExponential:
            A_fit, T1_fit, B_fit, T2_fit, C_fit = opt
            A_std, T1_std, B_std, T2_std, C_std = err
        else:
            A_fit, T1_fit, B_fit = opt
            A_std, T1_std, B_std = err

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        plt.plot(TimeList[plot_ind[0]], P_plot, 'o')
        plt.plot(TimeFitList[plot_ind[0]], fit_plot)
        if MeasurementType in ('t1', 't1 no ref'):
            if FitDoubleExponential:
                plt.title('TR=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                          'Tqp=%.3G$\pm$%.2Gus, lambda=%.3G$\pm$%.2G' % (
                              TR_fit / 1000, TR_std / 1000, A_fit, B_fit, Tqp_fit / 1000, Tqp_std / 1000, lamb_fit,
                              lamb_std))
            elif FitTwoExponential:
                plt.title('T_Exp1=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                          'T_Exp2=%.3G$\pm$%.2Gus, C=%.3G$\pm$%.2G' % (
                              T1_fit / 1000, T1_std / 1000, -A_fit * PopulationConversionConst[1],
                              -B_fit * PopulationConversionConst[1], T2_fit / 1000, T2_std / 1000, r_to_P(C_fit),
                              C_std))
            else:
                plt.title('T1=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                    T1_fit / 1000, T1_std / 1000, -A_fit * PopulationConversionConst[1], r_to_P(B_fit)))
        elif MeasurementType in ('transient', 'transient no ref'):
            if FitDoubleExponential:
                plt.title('TR=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                          'Tqp=%.3G$\pm$%.2Gus, lambda=%.3G$\pm$%.2G' % (
                              TR_fit / 1000, TR_std / 1000, A_fit, B_fit, Tqp_fit / 1000, Tqp_std / 1000, lamb_fit,
                              lamb_std))
            elif FitTwoExponential:
                plt.title('T_Exp1=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                          'T_Exp2=%.3G$\pm$%.2Gus, C=%.3G$\pm$%.2G' % (
                              T1_fit / 1000, T1_std / 1000, A_fit, B_fit, T2_fit / 1000, T2_std / 1000, C_fit,
                              C_std))
            else:
                plt.title('T_transient=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                    T1_fit / 1000, T1_std / 1000, -A_fit * PopulationConversionConst[1], r_to_P(B_fit)))
        plt.xlabel('Time(ns)', fontsize='x-large')
        plt.ylabel(plot_level, fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if LogScale:
            ax.set_yscale('log')

    if ShowFig:
        plt.show()
    else:
        plt.close('all')

    return {'opt': opt, 'cov': cov, 'ParamList': ParamList}


if __name__ == '__main__':
    DataFolderName = '11112019_back to waveguide'
    DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2020/01\Data_0122\\'
    BackgroundFolder = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
    BackgroundFile = []
    Plus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    Minus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    BackgroundFile = 'power spectroscopy_105.hdf5'
    RabiFile = 't1_P2_P1_30.hdf5'
    IQModFreq = 0.05
    CircleCorrection = True
    CorrectionParam = [1, 0.044, 0.737, 0.037]
    PhaseSlope = 326.7041108065019
    PhaseReferenceFreq = 4.105
    Calibration = True
    FitCorrectedR = False
    LimitTimeRange = False
    RotateComplex = False
    FitDoubleExponential = False
    FitTwoExponential = False
    LogScale = False
    SaveFig = False
    ShowFig = True
    StartTime = 0.5e3
    EndTime = 40e3
    PopulationConversionConst = [1., 1. / 0.8694655525612803]
    FitDict = plotMultiPopulationTSweep(DataPath, RabiFile, BackgroundFolder=BackgroundFolder,
                                        BackgroundFile=BackgroundFile,
                                        IQModFreq=IQModFreq, PopulationConversionConst=PopulationConversionConst,
                                        FitTwoExponential=FitTwoExponential,
                                        CircleCorrection=CircleCorrection,
                                        CorrectionParam=CorrectionParam,
                                        PhaseSlope=PhaseSlope, PhaseReferenceFreq=PhaseReferenceFreq,
                                        Calibration=Calibration,
                                        FitCorrectedR=FitCorrectedR, LimitTimeRange=LimitTimeRange,
                                        RotateComplex=RotateComplex,
                                        StartTime=StartTime, EndTime=EndTime, FitDoubleExponential=FitDoubleExponential,
                                        SaveFig=SaveFig, ShowFig=ShowFig, LogScale=LogScale)
    print(FitDict)
