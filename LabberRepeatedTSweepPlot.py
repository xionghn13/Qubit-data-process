import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, DoubleExp_curve, rabi_curve, FitTransientTime, AutoRotate
import ExtractDataFunc as edf
import h5py
import sklearn.metrics

def plotLabberRepeatedTSweepPlot(DataPath, RabiFileList, BackgroundFile='calibration_5.hdf5',
                                 IQModFreq=0.05, PhaseSlope=326.7, PhaseReferenceFreq=4.105,
                                 Calibration=False, FitCorrectedR=True, RotateComplex=True,
                                 LogScale=False, FitDoubleExponential=False, PlotNumber=11, MinPlotInd=0, MaxPlotInd=11,
                                 PlotInd=[0, 1, 2, 3], T2MaxTime=2e4):
    if not isinstance(RabiFileList, list):
        RabiFileList = [RabiFileList]
    NumFile = len(RabiFileList)
    CounterArray = np.zeros([NumFile, ])
    # analyze background file

    if Calibration:
        [BackFreq, BackComplex] = edf.readFSweepLabber(DataPath + BackgroundFile)
        BackPower = edf.readQubitPowerLabber(DataPath + BackgroundFile)
        BackPowerStr = str(BackPower)

    LastCount = 0
    for i, RabiFile in enumerate(RabiFileList):
        RabiFileStrList = RabiFile[:-5].split('_')
        MeasurementType = RabiFileStrList[0]
        if MeasurementType == 't1' and RabiFileStrList[1] == 't2':
            MeasurementType = 't1t2interleaved'
        ReadoutFreq = edf.readQubitFreqLabber(DataPath + RabiFile)
        ReadoutPower = edf.readQubitPowerLabber(DataPath + RabiFile)
        QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)
        # Counter = edf.readPumpPowerLabber(DataPath + RabiFile)

        if MeasurementType in ('rabi', 'transient'):
            [Time, Counter, ComplexVoltage] = edf.readRepeatedRabiSweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't1':
            [Time, Counter, ComplexVoltage] = edf.readRepeatedT1SweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't2':
            [Time, Counter, ComplexVoltage] = edf.readRepeatedT2SweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't1t2interleaved':
            [Time, Counter, ComplexVoltageT1, ComplexVoltageT2] = edf.readRepeatedT1T2InterleavedSweepLabber(
                DataPath + RabiFile)
        Counter += LastCount + 1
        LastCount = Counter[-1]
        if MeasurementType == 't1t2interleaved':
            ComplexVoltageT1Normalized = ComplexVoltageT1 * 10 ** (- ReadoutPower / 20)
            ComplexVoltageT2Normalized = ComplexVoltageT2 * 10 ** (- ReadoutPower / 20)
            if Calibration:
                RComplexT1 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT1, BackFreq,
                                                            BackComplex,
                                                            BackPower)
                RComplexT2 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT2, BackFreq,
                                                            BackComplex,
                                                            BackPower)
            else:
                RComplexT1 = ComplexVoltageT1Normalized
                RComplexT2 = ComplexVoltageT2Normalized
            if T2MaxTime > 0:
                TimeInd = Time < T2MaxTime
                TimeT2 = Time[TimeInd]
                RComplexT2 = RComplexT2[TimeInd]
            else:
                TimeT2 = Time
        else:
            ComplexVoltageNormalized = ComplexVoltage * 10 ** (- ReadoutPower / 20)
            if Calibration:
                RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltage, BackFreq,
                                                          BackComplex,
                                                          BackPower)
            else:
                RComplex = ComplexVoltageNormalized

        for j, trial in enumerate(Counter):
            if RotateComplex:
                if MeasurementType == 't1t2interleaved':
                    RComplexT1[:, j] = AutoRotate(RComplexT1[:, j])
                    RComplexT2[:, j] = AutoRotate(RComplexT2[:, j])
                else:
                    RComplex[:, j] = AutoRotate(RComplex[:, j])
            if MeasurementType == 't1t2interleaved':
                y_dataList = [np.array(RComplexT1[:, j].real, dtype='float64'),
                              np.array(RComplexT2[:, j].real, dtype='float64')]
            else:
                y_data = np.array(RComplex[:, j].real, dtype='float64')
            x_data = np.array(Time, dtype='float64')

            if MeasurementType in ('t1', 'transient', 't2'):
                B_guess = y_data[-1]
                A_guess = y_data[0] - B_guess
                T1_guess = x_data[-1] / 2
                bounds = (
                    (-2, 1, -1),
                    (2, 1e6, 1)
                )
                if FitDoubleExponential:
                    try:
                        opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                                             p0=[A_guess, T1_guess, B_guess, T1_guess * 0.1, 0.5],
                                             maxfev=300000)
                        print('guess = %s' % str([A_guess, T1_guess, B_guess, T1_guess * 0.1, 1]))
                        print('Double exp fit opt = %s' % str(opt))
                    except RuntimeError:
                        print("Error - curve_fit failed")
                        opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
                        cov = np.zeros([len(opt), len(opt)])

                    A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
                    FitR = DoubleExp_curve(Time, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                    ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
                else:
                    opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
                    A_fit, T1_fit, B_fit = opt
                    FitR = T1_curve(Time, A_fit, T1_fit, B_fit)
                    ParamList = ['A', 'Decay time/ns', 'B']
                R2 = sklearn.metrics.r2_score(y_data, FitR)
            elif MeasurementType == 't1t2interleaved':
                FitRt1t2List = []
                optList = []
                covList = []
                for ind in range(2):
                    if ind == 1:
                        x_data = TimeT2
                    y_data = y_dataList[ind]
                    B_guess = y_data[-1]
                    A_guess = y_data[0].real - B_guess
                    T1_guess = x_data[-1] / 10
                    bounds = (
                        (-2, 1, -1),
                        (2, 1e6, 1)
                    )
                    opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
                    A_fit, T1_fit, B_fit = opt
                    FitR = T1_curve(Time, A_fit, T1_fit, B_fit)
                    FitRt1t2List += [FitR]
                    optList += [opt]
                    covList += [cov]
                ParamList = ['A', 'Decay time/ns', 'B']
            elif MeasurementType in ('rabi'):
                B_guess = y_data.mean()
                A_guess = y_data[0] - B_guess
                T1_guess = x_data[-1]
                MaxInd = y_data.argmax()
                MinInd = y_data.argmin()
                Tpi_guess = np.abs(x_data[MaxInd] - x_data[MinInd])
                phi0_guess = 0
                guess = ([A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess])
                bounds = (
                    (-2, 1, -1, 1, - np.pi / 2),
                    (2, np.inf, 1, np.inf, np.pi / 2)
                )
                opt, cov = curve_fit(rabi_curve, x_data, y_data, p0=guess, bounds=bounds)
                A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit = opt
                FitR = rabi_curve(Time, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
                ParamList = ['A', 'Decay time/ns', 'B', 'Tpi', 'phi0']

            if i == 0 and j == 0:
                TimeList = []
                RComplexList = []
                FitRList = []
                CounterArray = np.array([trial])
                if MeasurementType == 't1t2interleaved':
                    OptMatrixList = [np.reshape(optList[0], (len(optList[0]), 1)),
                                     np.reshape(optList[1], (len(optList[1]), 1))]
                    ErrMatrixList = [np.reshape(np.sqrt(covList[0].diagonal()), (len(optList[0]), 1)),
                                     np.reshape(np.sqrt(covList[1].diagonal()), (len(optList[0]), 1))]
                else:
                    R2Array = np.array([R2])
                    OptMatrix = np.reshape(opt, (len(opt), 1))
                    ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
            else:
                if MeasurementType == 't1t2interleaved':
                    for ind in range(2):
                        OptMatrixList[ind] = np.concatenate(
                            (OptMatrixList[ind], np.reshape(optList[ind], (len(opt), 1))),
                            axis=1)
                        ErrMatrixList[ind] = np.concatenate(
                            (ErrMatrixList[ind], np.reshape(np.sqrt(covList[ind].diagonal()), (len(opt), 1))), axis=1)
                else:
                    OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
                    ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)
                    R2Array = np.concatenate((R2Array, np.array([R2])))
                CounterArray = np.concatenate((CounterArray, np.array([trial])))


            TimeList.append(Time)
            if MeasurementType == 't1t2interleaved':
                RComplexList.append([RComplexT1[:, j], RComplexT2[:, j]])
                FitRList.append(FitRt1t2List)
                for ind in range(2):
                    for i_p, par in enumerate(optList[ind]):
                        if np.abs(ErrMatrixList[ind][i_p, -1]) > np.abs(par) or optList[ind][1] > 5 * Time[-1]:
                            OptMatrixList[ind][i_p, -1] = np.nan
                        else:
                            OptMatrixList[ind][i_p, -1] = par
            else:
                RComplexList.append(RComplex[:, j])
                FitRList.append(FitR)
                for i_p, par in enumerate(opt):
                    if np.abs(ErrMatrix[i_p, -1]) > np.abs(par) or opt[1] > 5 * Time[-1]:
                        OptMatrix[i_p, -1] = np.nan
                    else:
                        OptMatrix[i_p, -1] = par

    limit = 1.7

    NumPoints = len(TimeList)
    fig, ax = plt.subplots()
    for i in range(NumPoints):
        if MeasurementType == 't1t2interleaved':
            sym = ['-', '--']
            for ind in range(2):
                plt.plot(np.real(RComplexList[i][ind]), np.imag(RComplexList[i][ind]), sym[ind])
        else:
            plt.plot(np.real(RComplexList[i]), np.imag(RComplexList[i]))
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
    for i in np.linspace(0, NumPoints - 1, min(PlotNumber, NumPoints)):
        i = int(i)
        if MeasurementType == 't1t2interleaved':
            for ind in range(2):
                if ind == 1 and T2MaxTime > 0:
                    plt.plot(TimeT2, np.real(RComplexList[i][ind]))
                else:
                    plt.plot(TimeList[i], np.real(RComplexList[i][ind]))
                plt.plot(TimeList[i], FitRList[i][ind], '--')
        else:
            if LogScale:
                plt.plot(TimeList[i], FitRList[i] - opt[2])
                plt.plot(TimeList[i], np.real(RComplexList[i]) - opt[2], 'o')
            else:
                plt.plot(TimeList[i], FitRList[i] - opt[2])
                plt.plot(TimeList[i], np.real(RComplexList[i]) - opt[2], 'o')

    plt.xlabel('Time/ns', fontsize='x-large')
    plt.ylabel('Re', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    if LogScale:
        ax.set_yscale('log')

    if max(PlotIndex) < NumPoints:
        for i in PlotIndex:
            fig, ax = plt.subplots()
            if MeasurementType == 't1t2interleaved':
                for ind in range(2):
                    if ind == 1 and T2MaxTime > 0:
                        plt.plot(TimeT2, np.real(RComplexList[i][ind]))
                    else:
                        plt.plot(TimeList[i], np.real(RComplexList[i][ind]))
                    plt.plot(TimeList[i], FitRList[i][ind], '--')
            else:
                plt.plot(TimeList[i], np.real(RComplexList[i]), 'o')
                plt.plot(TimeList[i], FitRList[i])
            plt.xlabel('Time/ns', fontsize='x-large')
            plt.ylabel('Re', fontsize='x-large')
            plt.title('Plot index=%d' % i)
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

    fig, ax = plt.subplots()
    plotInd = 1
    # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
    if MeasurementType == 't1t2interleaved':
        avgList = []
        stdList = []
        for ind in range(2):
            avgList += [np.mean(OptMatrixList[ind][plotInd, :] / 1000)]
            if NumPoints > 1:
                stdList += [np.std(OptMatrixList[ind][plotInd, :] / 1000)]
            else:
                stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
            ax.errorbar(CounterArray, OptMatrixList[ind][plotInd, :] / 1000, yerr=ErrMatrixList[ind][plotInd, :] / 1000,
                        fmt='o')
        plt.legend(['T1', 'T2echo'])
        plt.title('T1=%.3G$\pm$%.2Gus, T2=%.3G$\pm$%.2Gus' % (avgList[0], stdList[0], avgList[1], stdList[1]))
    else:
        ax.errorbar(CounterArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
    plt.xlabel('Trial #', fontsize='x-large')
    plt.ylabel('Decay time/us', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    if LogScale:
        ax.set_yscale('log')

    if MaxPlotInd != 0 or MinPlotInd != 0:
        fig, ax = plt.subplots()
        plotInd = 1
        # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
        if MeasurementType == 't1t2interleaved':
            if FitDoubleExponential:
                print('Error: t1t2interleaved with double exponential fit has not been completed yet.')
            avgList = []
            stdList = []
            for ind in range(2):
                avgList += [np.mean(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                if NumPoints > 1:
                    stdList += [np.std(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                else:
                    stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
                ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd],
                            OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000,
                            yerr=ErrMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, fmt='o')
            plt.legend(['T1', 'T2echo'])
            for ind in range(2):
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
            plt.title('T1=%.3G$\pm$%.2Gus, T2=%.3G$\pm$%.2Gus' % (avgList[0], stdList[0], avgList[1], stdList[1]))
        else:
            avg = np.mean(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
            std = np.std(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
            ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                        yerr=ErrMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                        fmt='o')
            plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
            if FitDoubleExponential:
                plt.title('TR=%.3G$\pm$%.2Gus' % (avg, std))
            else:
                plt.title('T1=%.3G$\pm$%.2Gus' % (avg, std))

        plt.xlabel('Trial #', fontsize='x-large')
        plt.ylabel('Decay time/us', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if LogScale:
            ax.set_yscale('log')

        if FitDoubleExponential:
            fig, ax = plt.subplots()
            plotInd = 3
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            if MeasurementType == 't1t2interleaved':
                avgList = []
                stdList = []
                for ind in range(1):
                    avgList += [np.mean(OptMatrixList[ind][plotInd, :MaxPlotInd] / 1000)]
                    if NumPoints > 1:
                        stdList += [np.std(OptMatrixList[ind][plotInd, :MaxPlotInd] / 1000)]
                    else:
                        stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
                    ax.errorbar(CounterArray[:MaxPlotInd], OptMatrixList[ind][plotInd, :MaxPlotInd] / 1000,
                                yerr=ErrMatrixList[ind][plotInd, :MaxPlotInd] / 1000, fmt='o')
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
                # plt.legend(['T1', 'T2echo'])
                plt.title('Tqp=%.3G$\pm$%.2Gus' % (avgList[0], stdList[0]))
            else:
                avg = np.mean(OptMatrix[plotInd, :MaxPlotInd]) / 1000
                std = np.std(OptMatrix[plotInd, :MaxPlotInd]) / 1000
                ax.errorbar(CounterArray[:MaxPlotInd], OptMatrix[plotInd, :MaxPlotInd] / 1000,
                            yerr=ErrMatrix[plotInd, :MaxPlotInd] / 1000,
                            fmt='o')
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
                plt.title('Tqp=%.3G$\pm$%.2Gus' % (avg, std))

            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('Decay time/us', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            if LogScale:
                ax.set_yscale('log')

            fig, ax = plt.subplots()
            plotInd = 4
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            if MeasurementType == 't1t2interleaved':
                avgList = []
                stdList = []
                for ind in range(1):
                    avgList += [np.mean(OptMatrixList[ind][plotInd, :MaxPlotInd])]
                    if NumPoints > 1:
                        stdList += [np.std(OptMatrixList[ind][plotInd, :MaxPlotInd])]
                    else:
                        stdList += [ErrMatrixList[ind][plotInd, 0]]
                    ax.errorbar(CounterArray[:MaxPlotInd], OptMatrixList[ind][plotInd, :MaxPlotInd],
                                yerr=ErrMatrixList[ind][plotInd, :MaxPlotInd], fmt='o')
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd], '--')
                # plt.legend(['T1', 'T2echo'])
                plt.title('nqp=%.3G$\pm$%.2G' % (avgList[0], stdList[0]))
            else:
                avg = np.mean(OptMatrix[plotInd, :MaxPlotInd])
                std = np.std(OptMatrix[plotInd, :MaxPlotInd])
                ax.errorbar(CounterArray[:MaxPlotInd], OptMatrix[plotInd, :MaxPlotInd],
                            yerr=ErrMatrix[plotInd, :MaxPlotInd],
                            fmt='o')
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd], '--')
                plt.title('nqp=%.3G$\pm$%.2G' % (avg, std))

            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('lambda', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            if LogScale:
                ax.set_yscale('log')

        if not MeasurementType == 't1t2interleaved':
            fig, ax = plt.subplots()
            if MeasurementType == 't1t2interleaved':
                print('Not done yet for t1 t2 R^2')
            else:
                plt.plot(CounterArray, R2Array, 'o')
            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('R^2', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()


    if MeasurementType == 'rabi':
        fig, ax = plt.subplots()
        plotInd = 2
        # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
        ax.errorbar(CounterArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
        plt.xlabel('Trial #', fontsize='x-large')
        plt.ylabel('T_pi (ns)', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        # plt.ylim(0, 20000)
    plt.show()


if __name__ == '__main__':
    # DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium042619/'
    DataPath = 'C:/Users/admin\Labber\Data/2019/08\Data_0812/'
    BackgroundFile = 'calibration_5.hdf5'

    # RabiFileList = [
    #     'transient_9.hdf5',
    #     'transient_8.hdf5',
    #     'transient_7.hdf5',
    #     'transient_6.hdf5',
    #
    # ]
    RabiFileList = [
        't1_t2_interleaved_2019-08-12-11-55-08.hdf5',
        # 't1_2019-06-17-20-53-22.hdf5',
    ]

    IQModFreq = 0.05
    # FitForGamma = True
    Gamma_r = 2.5 * np.pi * 2
    FitCorrectedR = True
    LogScale = False
    Calibration = False
    RotateComplex = True
    FitDoubleExponential = False
    PlotNumber = 50  # fit plot
    MinPlotInd = 0
    MaxPlotInd = 50
    PlotIndex = [0, 9]
    T2MaxTime = 300e3  # ns

    PhaseSlope = 326.7041108065019
    PhaseReferenceFreq = 4.105
    plotLabberRepeatedTSweepPlot(DataPath, RabiFileList, BackgroundFile=BackgroundFile,
                                 IQModFreq=IQModFreq, PhaseSlope=PhaseSlope, PhaseReferenceFreq=PhaseReferenceFreq,
                                 Calibration=Calibration, FitCorrectedR=FitCorrectedR, RotateComplex=RotateComplex,
                                 LogScale=LogScale, FitDoubleExponential=FitDoubleExponential, PlotNumber=PlotNumber,
                                 MinPlotInd=MinPlotInd, MaxPlotInd=MaxPlotInd, PlotInd=PlotIndex,
                                 T2MaxTime=T2MaxTime)
