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
import ThermalPhotonDephasing
import pickle


def plotLabberRepeatedTSweepPlot(DataPath, RabiFileList, BackgroundFolder='', BackgroundFile='calibration_5.hdf5',
                                 IQModFreq=0.05, PhaseSlope=326.7, PhaseReferenceFreq=4.105,
                                 Calibration=False, FitCorrectedR=True, RotateComplex=True,
                                 LogScale=False, FitDoubleExponential=False, PlotNumber=11, MinPlotInd=0, MaxPlotInd=11,
                                 PlotIndex=[0, 1, 2, 3], T2MaxTime=2e4, ShowFig=True,
                                 chi_kappa_f=[5.207 * 2 * np.pi, 3.86 * 2 * np.pi, 7.96e9]):
    if not isinstance(RabiFileList, list):
        RabiFileList = [RabiFileList]
    NumFile = len(RabiFileList)
    CounterArray = np.zeros([NumFile, ])
    # analyze background file

    if Calibration:
        if BackgroundFile.startswith('RefRabiCal'):
            # use rabi to calibrate rabi. Plus or minus are merged in one array.
            [BackFreq, BackPower, BackComplex] = edf.readRefRabiCalLabber(BackgroundFolder + BackgroundFile)
        else:
            [BackFreq, BackComplex] = edf.readFSweepLabber(BackgroundFolder + BackgroundFile)
            BackPower = edf.readReadoutPowerLabber(BackgroundFolder + BackgroundFile)
            BackPowerStr = str(BackPower)

    LastCount = 0
    for i, RabiFile in enumerate(RabiFileList):
        RabiFileStrList = RabiFile[:-5].split('_')
        MeasurementType = RabiFileStrList[0]
        if MeasurementType == 't1' and RabiFileStrList[1] == 't2':
            MeasurementType = 't1t2interleaved'
        elif MeasurementType == 't1' and RabiFileStrList[1] == 'ramsey' and RabiFileStrList[2] == 'echo':
            MeasurementType = 't1_ramsey_echo_interleaved'
        elif MeasurementType == 'int':
            MeasurementType = 't1'
        elif RabiFileStrList[0:4] == ['ref', 't1', 't2', 'interleaved']:
            MeasurementType = 'ref_t1_t2_interleaved'
        elif RabiFileStrList[0:4] == ['ref', 't1', 'ramsey', 'echo']:
            MeasurementType = 'ref_t1_ramsey_echo_interleaved'
        if RabiFileStrList[0] == 'T1':
            ReadoutPower = edf.readSetup2ReadoutPowerLabber(DataPath + RabiFile)
        else:
            ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
            ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
            QubitFreq = edf.readPumpFreqLabber(DataPath + RabiFile)

        # Counter = edf.readPumpPowerLabber(DataPath + RabiFile)
        if MeasurementType in ('rabi', 'transient'):
            [Time, Counter, ComplexVoltage] = edf.readRepeatedRabiSweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't1':
            if RabiFileStrList[0] == 'int':
                [Time, Counter, ComplexVoltage] = edf.readRepeatedIntegratedT1SweepLabber(DataPath + RabiFile)
            else:
                [Time, Counter, ComplexVoltage] = edf.readRepeatedT1SweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't2':
            [Time, Counter, ComplexVoltage] = edf.readRepeatedT2SweepLabber(DataPath + RabiFile)
        elif MeasurementType == 't1t2interleaved':
            [Time, Counter, ComplexVoltageT1, ComplexVoltageT2] = edf.readRepeatedT1T2InterleavedSweepLabber(
                DataPath + RabiFile)
        elif MeasurementType == 't1_ramsey_echo_interleaved':
            [Time, Counter, ComplexVoltageT1, ComplexVoltageRamsey,
             ComplexVoltageEcho] = edf.readRepeatedT1RamseyEchoInterleavedSweepLabber(DataPath + RabiFile)
        elif MeasurementType == 'ref_t1_t2_interleaved':
            [Time, Counter, ReadoutFreq, ComplexVoltageT1,
             ComplexVoltageT2] = edf.readReferencedRepeatedT1T2InterleavedSweepLabber(
                DataPath + RabiFile)
        elif MeasurementType == 'ref_t1_ramsey_echo_interleaved':
            [Time, Counter, ReadoutFreq, ComplexVoltageT1,
             ComplexVoltageRamsey, ComplexVoltageEcho] = edf.readReferencedRepeatedT1RamseyEchoInterleavedSweepLabber(
                DataPath + RabiFile)
        elif RabiFileStrList[0] == 'T1' and RabiFileStrList[1] == 'T2E':
            MeasurementType = 't1t2interleaved'
            [Time, Counter, ComplexVoltageT1, ComplexVoltageT2] = edf.readSetup2RepeatedT1T2InterleavedSweepLabber(
                DataPath + RabiFile)
        elif MeasurementType == 'T1':
            MeasurementType = 't1'
            [Time, Counter, ComplexVoltage] = edf.readSetup2RepeatedT1SweepLabber(DataPath + RabiFile)
        Counter += LastCount + 1
        LastCount = Counter[-1]
        # print(ComplexVoltageT1[:, 0, 0])
        if MeasurementType in ('t1t2interleaved', 'ref_t1_t2_interleaved'):
            ComplexVoltageT1Normalized = ComplexVoltageT1 * 10 ** (- ReadoutPower / 20)
            ComplexVoltageT2Normalized = ComplexVoltageT2 * 10 ** (- ReadoutPower / 20)
            if Calibration:
                RComplexT1 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT1, BackFreq,
                                                            BackComplex,
                                                            BackPower)
                RComplexT2 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT2, BackFreq,
                                                            BackComplex,
                                                            BackPower)
                # print('T2')
                # print(RComplexT2[1, 0, 0])
            else:
                RComplexT1 = ComplexVoltageT1Normalized
                RComplexT2 = ComplexVoltageT2Normalized
            if T2MaxTime > 0:
                TimeInd = Time < T2MaxTime
                TimeT2 = Time[TimeInd]
                if MeasurementType.startswith('ref'):
                    RComplexT2 = RComplexT2[TimeInd, :, :]
                else:
                    RComplexT2 = RComplexT2[TimeInd]
            else:
                TimeT2 = Time
        elif MeasurementType in ('t1_ramsey_echo_interleaved', 'ref_t1_ramsey_echo_interleaved'):
            ComplexVoltageT1Normalized = ComplexVoltageT1 * 10 ** (- ReadoutPower / 20)
            ComplexVoltageRamseyNormalized = ComplexVoltageRamsey * 10 ** (- ReadoutPower / 20)
            ComplexVoltageEchoNormalized = ComplexVoltageEcho * 10 ** (- ReadoutPower / 20)
            if Calibration:
                RComplexT1 = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageT1, BackFreq,
                                                            BackComplex, BackPower)
                RComplexRamsey = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageRamsey,
                                                                BackFreq, BackComplex, BackPower)
                RComplexEcho = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, ComplexVoltageEcho,
                                                              BackFreq, BackComplex, BackPower)
                # print(RComplexT1[:, 0, 0])
            else:
                RComplexT1 = ComplexVoltageT1Normalized
                RComplexRamsey = ComplexVoltageRamseyNormalized
                RComplexEcho = ComplexVoltageEchoNormalized
            if T2MaxTime > 0:
                TimeInd = Time < T2MaxTime
                TimeT2 = Time[TimeInd]
                if MeasurementType.startswith('ref'):
                    RComplexRamsey = RComplexRamsey[TimeInd, :, :]
                    RComplexEcho = RComplexEcho[TimeInd, :, :]
                else:
                    RComplexRamsey = RComplexRamsey[TimeInd]
                    RComplexEcho = RComplexEcho[TimeInd]
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

        # fit decay
        if Calibration:
            RotateComplex = False
        TimeFit = np.linspace(Time.min(), Time.max(), 200)
        for j, trial in enumerate(Counter):
            if RotateComplex:
                if MeasurementType == 't1t2interleaved':
                    RComplexT1[:, j] = AutoRotate(RComplexT1[:, j])
                    RComplexT2[:, j] = AutoRotate(RComplexT2[:, j])
                elif MeasurementType == 't1_ramsey_echo_interleaved':
                    RComplexT1[:, j] = AutoRotate(RComplexT1[:, j])
                    RComplexRamsey[:, j] = AutoRotate(RComplexRamsey[:, j])
                    RComplexEcho[:, j] = AutoRotate(RComplexEcho[:, j])
                else:
                    RComplex[:, j] = AutoRotate(RComplex[:, j])
            if MeasurementType == 't1t2interleaved':
                y_dataList = [np.array(RComplexT1[:, j].real, dtype='float64'),
                              np.array(RComplexT2[:, j].real, dtype='float64')]
            elif MeasurementType == 't1_ramsey_echo_interleaved':
                y_dataList = [np.array(RComplexT1[:, j].real, dtype='float64'),
                              np.array(RComplexRamsey[:, j].real, dtype='float64'),
                              np.array(RComplexEcho[:, j].real, dtype='float64')]
            elif MeasurementType == 'ref_t1_t2_interleaved':
                y_dataList = [np.array(RComplexT1[:, j, 0].real, dtype='float64'),
                              np.array(RComplexT2[:, j, 0].real, dtype='float64')]
            elif MeasurementType == 'ref_t1_ramsey_echo_interleaved':
                y_dataList = [np.array(RComplexT1[:, j, 0].real, dtype='float64'),
                              np.array(RComplexRamsey[:, j, 0].real, dtype='float64'),
                              np.array(RComplexEcho[:, j, 0].real, dtype='float64')]
            else:
                y_data = np.array(RComplex[:, j].real, dtype='float64')
            x_data = np.array(Time, dtype='float64')

            if MeasurementType in ('t1', 'transient', 't2'):
                B_guess = y_data[-1]
                A_guess = y_data[0] - B_guess
                T1_guess = x_data[-1] / 10
                bounds = (
                    (-2, 1, -1),
                    (2, 1e6, 1)
                )
                if FitDoubleExponential:
                    try:
                        opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                                             p0=[A_guess, T1_guess, B_guess, T1_guess * 0.1, 0.5],
                                             maxfev=300000)
                    except RuntimeError:
                        print("Error - curve_fit failed")
                        opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
                        cov = np.zeros([len(opt), len(opt)])

                    A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
                    FitR = DoubleExp_curve(TimeFit, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                    y_pred = DoubleExp_curve(Time, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                    ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
                else:
                    opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
                    A_fit, T1_fit, B_fit = opt
                    FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
                    y_pred = T1_curve(Time, A_fit, T1_fit, B_fit)
                    ParamList = ['A', 'Decay time/ns', 'B']
                R2 = sklearn.metrics.r2_score(y_data, y_pred)
            elif MeasurementType in ('t1t2interleaved', 'ref_t1_t2_interleaved'):
                FitRt1t2List = []
                optList = []
                covList = []
                for ind in range(2):
                    if ind == 1:
                        x_data = TimeT2
                    y_data = y_dataList[ind]
                    B_guess = y_data[-1]
                    A_guess = y_data[0] - B_guess
                    T1_guess = x_data[-1] / 10
                    Tqp_guess = T1_guess / 10
                    bounds = (
                        (-10 * abs(A_guess), 0, -10 * abs(B_guess)),
                        (10 * abs(A_guess), x_data[-1] * 10, 10 * abs(B_guess))
                    )
                    # T1_guess = 180e3
                    # Tqp_guess = 30e3
                    if ind == 0 and FitDoubleExponential:
                        try:
                            opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                                                 p0=[A_guess, T1_guess, B_guess, Tqp_guess, 0.5],
                                                 maxfev=300000)
                        except RuntimeError:
                            print("Error - curve_fit failed")
                            opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
                            cov = np.zeros([len(opt), len(opt)])

                        A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
                        FitR = DoubleExp_curve(TimeFit, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                        FitRt1t2List += [FitR]
                        optList += [opt]
                        covList += [cov]
                        ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
                    if ind == 1 or not FitDoubleExponential:
                        print(y_data[0])
                        opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], bounds=bounds,
                                             maxfev=30000)
                        A_fit, T1_fit, B_fit = opt
                        # print(T1_fit)
                        FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
                        FitRt1t2List += [FitR]
                        optList += [opt]
                        print(opt)
                        covList += [cov]
                        ParamList = ['A', 'Decay time/ns', 'B']
            elif MeasurementType in ('t1_ramsey_echo_interleaved', 'ref_t1_ramsey_echo_interleaved'):
                FitRt1t2List = []
                optList = []
                covList = []
                for ind in range(3):
                    if ind > 0:
                        x_data = TimeT2
                    y_data = y_dataList[ind]
                    B_guess = y_data[-1]
                    A_guess = y_data[0] - B_guess
                    T1_guess = x_data[-1] / 10
                    Tqp_guess = T1_guess / 10
                    MaxInd = y_data.argmax()
                    MinInd = y_data.argmin()
                    Tpi_guess = np.abs(x_data[MaxInd] - x_data[MinInd])
                    phi0_guess = 0
                    # T1_guess = 180e3
                    # Tqp_guess = 30e3

                    if ind == 0 and FitDoubleExponential:
                        try:
                            opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                                                 p0=[A_guess, T1_guess, B_guess, Tqp_guess, 0.5],
                                                 maxfev=300000)
                        except RuntimeError:
                            print("Error - curve_fit failed")
                            opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
                            cov = np.zeros([len(opt), len(opt)])

                        A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit = opt
                        FitR = DoubleExp_curve(TimeFit, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
                        FitRt1t2List += [FitR]
                        optList += [opt]
                        covList += [cov]
                        ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
                    elif ind == 1:
                        bounds = (
                            (A_guess - 0.5 * np.abs(A_guess), 0, - np.inf, 0, - np.inf),
                            (A_guess + 0.5 * np.abs(A_guess), T1_guess * 1e6, np.inf, T1_guess * 1e6, np.inf)
                        )
                        opt, cov = curve_fit(rabi_curve, x_data, y_data,
                                             p0=[A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess],
                                             bounds=bounds, maxfev=300000)
                        A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit = opt
                        FitR = rabi_curve(TimeFit, A_fit, T1_fit, B_fit, Tpi_fit, phi0_guess)
                        FitRt1t2List += [FitR]
                        optList += [opt]
                        covList += [cov]
                        ParamList = ['A', 'Decay time/ns', 'B', 'Tpi(ns)']
                    elif ind == 2 or not FitDoubleExponential:
                        bounds = (
                            (- np.inf, 1, - np.inf),
                            (np.inf, 1e6, np.inf)
                        )
                        opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], bounds=bounds,
                                             maxfev=30000)
                        A_fit, T1_fit, B_fit = opt
                        FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
                        FitRt1t2List += [FitR]
                        optList += [opt]
                        covList += [cov]
                        ParamList = ['A', 'Decay time(ns)', 'B']
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
                FitR = rabi_curve(TimeFit, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
                ParamList = ['A', 'Decay time(ns)', 'B', 'Tpi', 'phi0']

            # concatenate data
            if i == 0 and j == 0:
                TimeList = []
                TimeFitList = []
                RComplexList = []
                RefRComplexList = []
                FitRList = []
                CounterArray = np.array([trial])
                if MeasurementType in ('t1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                                       'ref_t1_ramsey_echo_interleaved'):
                    OptMatrixList = [np.reshape(opt_arr, (len(opt_arr), 1)) for opt_arr in optList]
                    ErrMatrixList = [np.reshape(np.sqrt(cov_mat.diagonal()), (len(cov_mat.diagonal()), 1)) for cov_mat
                                     in covList]
                else:
                    R2Array = np.array([R2])
                    OptMatrix = np.reshape(opt, (len(opt), 1))
                    ErrMatrix = np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))
            else:
                if MeasurementType in ('t1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                                       'ref_t1_ramsey_echo_interleaved'):
                    num_meas = len(OptMatrixList)
                    for ind in range(num_meas):
                        # print(OptMatrixList)
                        # print(optList[ind])
                        OptMatrixList[ind] = np.concatenate(
                            (OptMatrixList[ind], np.reshape(optList[ind], (len(optList[ind]), 1))), axis=1)
                        ErrMatrixList[ind] = np.concatenate(
                            (ErrMatrixList[ind], np.reshape(np.sqrt(covList[ind].diagonal()), (len(optList[ind]), 1))),
                            axis=1)
                else:
                    OptMatrix = np.concatenate((OptMatrix, np.reshape(opt, (len(opt), 1))), axis=1)
                    ErrMatrix = np.concatenate((ErrMatrix, np.reshape(np.sqrt(cov.diagonal()), (len(opt), 1))), axis=1)
                    R2Array = np.concatenate((R2Array, np.array([R2])))
                CounterArray = np.concatenate((CounterArray, np.array([trial])))

            TimeList.append(Time)
            TimeFitList.append(TimeFit)
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                if MeasurementType == 't1t2interleaved':
                    RComplexList.append([RComplexT1[:, j], RComplexT2[:, j]])
                elif MeasurementType == 'ref_t1_t2_interleaved':
                    RComplexList.append([RComplexT1[:, j, 0], RComplexT2[:, j, 0]])
                    RefRComplexList.append([RComplexT1[:, j, 1], RComplexT2[:, j, 1]])
                elif MeasurementType == 'ref_t1_ramsey_echo_interleaved':
                    RComplexList.append([RComplexT1[:, j, 0], RComplexRamsey[:, j, 0], RComplexEcho[:, j, 0]])
                    RefRComplexList.append([RComplexT1[:, j, 1], RComplexRamsey[:, j, 1], RComplexEcho[:, j, 1]])
                else:
                    RComplexList.append([RComplexT1[:, j], RComplexRamsey[:, j], RComplexEcho[:, j]])
                FitRList.append(FitRt1t2List)
                num_meas = len(OptMatrixList)
                for ind in range(num_meas):
                    for i_p, par in enumerate(optList[ind]):
                        if np.abs(ErrMatrixList[ind][i_p, -1]) > 0.4 * np.abs(par) or optList[ind][1] > 5 * Time[-1]:
                            OptMatrixList[ind][i_p, -1] = np.nan
                            # OptMatrixList[ind][i_p, -1] = par
                        else:
                            OptMatrixList[ind][i_p, -1] = par
                if MeasurementType == 't1_ramsey_echo_interleaved' and optList[1][1] > optList[2][1]:
                    OptMatrixList[1][1, -1] = np.nan

            else:
                RComplexList.append(RComplex[:, j])
                FitRList.append(FitR)
                for i_p, par in enumerate(opt):
                    if np.abs(ErrMatrix[i_p, -1]) > 0.4 * np.abs(par) or opt[1] > 5 * Time[-1]:
                        OptMatrix[i_p, -1] = np.nan
                        # OptMatrix[i_p, -1] = par
                    else:
                        OptMatrix[i_p, -1] = par

    # plot
    limit = 1.7

    NumPoints = len(TimeList)
    if ShowFig:
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        for i in range(NumPoints):
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                sym = ['-', '--', ':']
                num_meas = len(OptMatrixList)
                for ind in range(num_meas):
                    plt.plot(np.real(RComplexList[i][ind]), np.imag(RComplexList[i][ind]), sym[ind])
                    if MeasurementType.startswith('ref'):
                        plt.plot(np.real(RefRComplexList[i][ind]), np.imag(RefRComplexList[i][ind]), sym[ind])
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
        ax.grid(linestyle='--')
        for i in np.linspace(0, NumPoints - 1, min(PlotNumber, NumPoints)):
            i = int(i)
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                num_meas = len(OptMatrixList)
                for ind in range(num_meas):
                    if ind == 1 and T2MaxTime > 0:
                        plt.plot(TimeT2, np.real(RComplexList[i][ind]))
                        if MeasurementType.startswith('ref'):
                            plt.plot(TimeT2, np.real(RefRComplexList[i][ind]))
                    else:
                        plt.plot(TimeList[i], np.real(RComplexList[i][ind]))
                        if MeasurementType.startswith('ref'):
                            plt.plot(TimeList[i], np.real(RefRComplexList[i][ind]))
                    plt.plot(TimeFitList[i], FitRList[i][ind], '--')
            else:
                if LogScale:
                    plt.plot(TimeFitList[i], FitRList[i] - opt[2])
                    plt.plot(TimeList[i], np.real(RComplexList[i]) - opt[2], 'o')
                else:
                    plt.plot(TimeFitList[i], FitRList[i] - opt[2])
                    plt.plot(TimeList[i], np.real(RComplexList[i]) - opt[2], 'o')

        plt.xlabel('Time(ns)', fontsize='x-large')
        plt.ylabel('Re', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if LogScale:
            ax.set_yscale('log')

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        amp_correction = 10 ** (- ReadoutPower / 20) * 1e3
        if Calibration:
            amp_correction = 1
        for i in np.linspace(MinPlotInd, min(MaxPlotInd, NumPoints - 1), min(PlotNumber, NumPoints)):
            i = int(i)
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                ind = 0
                plt.plot(TimeList[i] / 1000, np.real(RComplexList[i][ind]) / amp_correction, 'o')
                if MeasurementType.startswith('ref'):
                    plt.plot(TimeList[i] / 1000, np.real(RefRComplexList[i][ind]) / amp_correction,
                             'o')
        # plt.legend(['base temperature', '50mK'])
        for i in np.linspace(MinPlotInd, min(MaxPlotInd, NumPoints - 1), min(PlotNumber, NumPoints)):
            i = int(i)
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                ind = 0
                plt.plot(TimeFitList[i] / 1000, FitRList[i][ind] / amp_correction, '--')
        plt.xlabel('Time(us)', fontsize='x-large')
        plt.ylabel('Re', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if LogScale:
            ax.set_yscale('log')

        if max(PlotIndex) < NumPoints:
            for i in PlotIndex:
                fig, ax = plt.subplots()
                ax.grid(linestyle='--')
                if MeasurementType in ('t1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                                       'ref_t1_ramsey_echo_interleaved'):
                    num_meas = len(OptMatrixList)
                    for ind in range(num_meas):
                        if ind == 1 and T2MaxTime > 0:
                            plt.plot(TimeT2, np.real(RComplexList[i][ind]))
                            if MeasurementType.startswith('ref'):
                                plt.plot(TimeT2, np.real(RefRComplexList[i][ind]))
                        else:
                            plt.plot(TimeList[i], np.real(RComplexList[i][ind]))
                            if MeasurementType.startswith('ref'):
                                plt.plot(TimeList[i], np.real(RefRComplexList[i][ind]))
                        plt.plot(TimeFitList[i], FitRList[i][ind], '--')
                else:
                    plt.plot(TimeList[i], np.real(RComplexList[i]), 'o')
                    plt.plot(TimeFitList[i], FitRList[i])
                plt.xlabel('Time(ns)', fontsize='x-large')
                plt.ylabel('Re', fontsize='x-large')
                plt.title('Plot index=%d' % i)
                plt.tick_params(axis='both', which='major', labelsize='x-large')
                plt.tight_layout()

        if MeasurementType in ('t1t2interleaved', 'ref_t1_t2_interleaved'):
            t2leg = ['T2echo']
            t2title = 'T2=%.3G$\pm$%.2Gus'
        else:
            t2leg = ['T2Ramsey', 'T2echo']
            t2title = 'T2Ramsey=%.3G$\pm$%.2Gus, T2echo=%.3G$\pm$%.2Gus'

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        plotInd = 1
        if MeasurementType in (
                't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                'ref_t1_ramsey_echo_interleaved'):
            avgList = []
            stdList = []
            avg_std_display = []
            num_meas = len(OptMatrixList)
            for ind in range(num_meas):
                avgList += [np.nanmean(OptMatrixList[ind][plotInd, :] / 1000)]
                if NumPoints > 1:
                    stdList += [np.nanstd(OptMatrixList[ind][plotInd, :] / 1000)]
                else:
                    stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
                ax.errorbar(CounterArray, OptMatrixList[ind][plotInd, :] / 1000,
                            yerr=ErrMatrixList[ind][plotInd, :] / 1000,
                            fmt='o')
                avg_std_display += [avgList[ind], stdList[ind]]
            plt.legend(['T1'] + t2leg)
            tit = 'T1=%.3G$\pm$%.2Gus, ' + t2title
            plt.title(tit % tuple(avg_std_display))
        else:
            ax.errorbar(CounterArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
        plt.xlabel('Trial #', fontsize='x-large')
        plt.ylabel('Decay time(us)', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if LogScale:
            ax.set_yscale('log')

        # plot t1 t2
        if MaxPlotInd != 0 or MinPlotInd != 0:
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plotInd = 1
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                avgList = []
                stdList = []
                avg_std_display = []
                num_meas = len(OptMatrixList)
                for ind in range(num_meas):
                    avgList += [np.nanmean(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                    if NumPoints > 1:
                        stdList += [np.nanstd(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                    else:
                        stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
                    ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd],
                                OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000,
                                yerr=ErrMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, fmt='o')
                    avg_std_display += [avgList[ind], stdList[ind]]
                if FitDoubleExponential:
                    T1Name = 'TR'
                else:
                    T1Name = 'T1'
                plt.legend([T1Name] + t2leg)
                tit = T1Name + '=%.3G$\pm$%.2Gus, ' + t2title
                plt.title(tit % tuple(avg_std_display))
                num_meas = len(OptMatrixList)
                for ind in range(num_meas):
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd],
                             OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
            else:
                avg = np.nanmean(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
                std = np.nanstd(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
                ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                            yerr=ErrMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                            fmt='o')
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
                if FitDoubleExponential:
                    plt.title('TR=%.3G$\pm$%.2Gus' % (avg, std))
                else:
                    plt.title('T1=%.3G$\pm$%.2Gus' % (avg, std))

            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('Decay time(us)', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            if LogScale:
                ax.set_yscale('log')

        # plot Ramsey frequency
        if MeasurementType in ('t1_ramsey_echo_interleaved', 'ref_t1_ramsey_echo_interleaved'):
            if MaxPlotInd != 0 or MinPlotInd != 0:
                fig, ax = plt.subplots()
                ax.grid(linestyle='--')
                plotInd = 3
                ind = 1
                Tpi_us = OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000
                Tpi_us_err = ErrMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000
                delta_f_kHz = 1000 / 2 / Tpi_us
                delta_f_kHz_err = Tpi_us_err / 2 / Tpi_us ** 2 * 1000
                avgList = []
                stdList = []
                avg_std_display = []
                avgList += [np.nanmean(delta_f_kHz)]
                if NumPoints > 1:
                    stdList += [np.nanstd(delta_f_kHz)]
                else:
                    stdList += [delta_f_kHz_err[0]]
                ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd], delta_f_kHz, yerr=delta_f_kHz_err, fmt='o')
                avg_std_display += [avgList[0], stdList[0]]
                tit = '$\delta f$=%.3G$\pm$%.2GkHz'
                plt.title(tit % tuple(avg_std_display))
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd],
                         delta_f_kHz, '--')
                plt.xlabel('Trial #', fontsize='x-large')
                plt.ylabel('$\delta f$(kHz)', fontsize='x-large')
                plt.tick_params(axis='both', which='major', labelsize='x-large')
                plt.tight_layout()
                if LogScale:
                    ax.set_yscale('log')

            # thermal photon
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved',
                    'ref_t1_t2_interleaved', 'ref_t1_ramsey_echo_interleaved') and not FitDoubleExponential:
                if MaxPlotInd != 0 or MinPlotInd != 0:
                    fig, ax = plt.subplots()
                    ax.grid(linestyle='--')
                    plotInd = 1
                    T1Array = OptMatrixList[0][plotInd, MinPlotInd:MaxPlotInd] / 1000
                    T2Array = OptMatrixList[-1][plotInd, MinPlotInd:MaxPlotInd] / 1000
                    TemperatureArray = np.array([])
                    [chi, kappa, f] = chi_kappa_f
                    for T1Ind in range(len(T1Array)):
                        T1 = T1Array[T1Ind]
                        T2 = T2Array[T1Ind]
                        Temp = 1e3 * ThermalPhotonDephasing.cavityThermalPhotonTemperature(chi, T1, T2, kappa, f,
                                                                                           PlotFigures=False)
                        TemperatureArray = np.append(TemperatureArray, Temp)
                    avg = np.nanmean(TemperatureArray)
                    std = np.nanstd(TemperatureArray)
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd], TemperatureArray, 'o')
                    plt.title('Photon temperature=%.3G$\pm$%.2GmK' % (avg, std))
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd], TemperatureArray, '--')
                    plt.xlabel('Trial #', fontsize='x-large')
                    plt.ylabel('Temperature(mK)', fontsize='x-large')
                    plt.tick_params(axis='both', which='major', labelsize='x-large')
                    plt.tight_layout()
                    if LogScale:
                        ax.set_yscale('log')

        if FitDoubleExponential:
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plotInd = 3
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                avgList = []
                stdList = []
                for ind in range(1):
                    avgList += [np.nanmean(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                    if NumPoints > 1:
                        stdList += [np.nanstd(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000)]
                    else:
                        stdList += [ErrMatrixList[ind][plotInd, 0] / 1000]
                    ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd],
                                OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000,
                                yerr=ErrMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, fmt='o')
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd],
                             OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd] / 1000, '--')
                # plt.legend(['T1', 'T2echo'])
                plt.title('Tqp=%.3G$\pm$%.2Gus' % (avgList[0], stdList[0]))
            else:
                avg = np.nanmean(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
                std = np.nanstd(OptMatrix[plotInd, MinPlotInd:MaxPlotInd]) / 1000
                ax.errorbar(CounterArray[:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                            yerr=ErrMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                            fmt='o')
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd] / 1000,
                         '--')
                plt.title('Tqp=%.3G$\pm$%.2Gus' % (avg, std))

            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('Decay time(us)', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            if LogScale:
                ax.set_yscale('log')

            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plotInd = 4
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            if MeasurementType in (
                    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                    'ref_t1_ramsey_echo_interleaved'):
                avgList = []
                stdList = []
                for ind in range(1):
                    avgList += [np.nanmean(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd])]
                    if NumPoints > 1:
                        stdList += [np.nanstd(OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd])]
                    else:
                        stdList += [ErrMatrixList[ind][plotInd, 0]]
                    ax.errorbar(CounterArray[MinPlotInd:MaxPlotInd], OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd],
                                yerr=ErrMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd], fmt='o')
                    plt.plot(CounterArray[MinPlotInd:MaxPlotInd],
                             OptMatrixList[ind][plotInd, MinPlotInd:MaxPlotInd], '--')
                # plt.legend(['T1', 'T2echo'])
                plt.title('nqp=%.3G$\pm$%.2G' % (avgList[0], stdList[0]))
            else:
                avg = np.nanmean(OptMatrix[plotInd, MinPlotInd:MaxPlotInd])
                std = np.nanstd(OptMatrix[plotInd, MinPlotInd:MaxPlotInd])
                ax.errorbar(CounterArray[:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd],
                            yerr=ErrMatrix[plotInd, MinPlotInd:MaxPlotInd],
                            fmt='o')
                plt.plot(CounterArray[MinPlotInd:MaxPlotInd], OptMatrix[plotInd, MinPlotInd:MaxPlotInd], '--')
                plt.title('nqp=%.3G$\pm$%.2G' % (avg, std))

            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('lambda', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            if LogScale:
                ax.set_yscale('log')

        if not MeasurementType in (
                't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved',
                'ref_t1_ramsey_echo_interleaved'):
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plt.plot(CounterArray, R2Array, 'o')
            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('R^2', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        if MeasurementType == 'rabi':
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plotInd = 2
            # plt.plot(CounterArray, OptMatrix[plotInd, :]/1000, 'o')
            ax.errorbar(CounterArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
            plt.xlabel('Trial #', fontsize='x-large')
            plt.ylabel('T_pi (ns)', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()
            # plt.ylim(0, 20000)

        plt.show()
    if MeasurementType in (
    't1t2interleaved', 't1_ramsey_echo_interleaved', 'ref_t1_t2_interleaved', 'ref_t1_ramsey_echo_interleaved'):
        return [CounterArray, OptMatrixList, ErrMatrixList]
    else:
        return [CounterArray, OptMatrix, ErrMatrix]


if __name__ == '__main__':
    # DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium042619/'
    DataFolderName = '11112019_back to waveguide'
    # DataFolderName = 'Data'
    DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2019/12\Data_1204\\'
    # DataPath = 'Z:\Projects\Transmon_Palmer\\2019\\10\Data_1017\\'
    BackgroundFolder = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
    BackgroundFile = 'power spectroscopy_83.hdf5'
    RabiFileList = [
        'transient_24.hdf5',
    ]

    IQModFreq = 0.05
    # FitForGamma = True
    Gamma_r = 2.5 * np.pi * 2
    FitCorrectedR = False
    LogScale = False
    Calibration = True
    RotateComplex = False
    FitDoubleExponential = False
    PlotNumber = 2  # fit plot
    MinPlotInd = 0
    MaxPlotInd = 8000
    PlotIndex = [0, 1]
    T2MaxTime = 1100e3  # ns

    # thermal photon
    chi = 6.281 * 2 * np.pi
    kappa = 4.08 * 2 * np.pi
    f = 7.97e9

    PhaseSlope = 326.7041108065019
    PhaseReferenceFreq = 4.105
    OutputList = plotLabberRepeatedTSweepPlot(DataPath, RabiFileList,
                                              BackgroundFolder=BackgroundFolder,
                                              BackgroundFile=BackgroundFile,
                                              IQModFreq=IQModFreq,
                                              PhaseSlope=PhaseSlope,
                                              PhaseReferenceFreq=PhaseReferenceFreq,
                                              Calibration=Calibration,
                                              FitCorrectedR=FitCorrectedR,
                                              RotateComplex=RotateComplex,
                                              LogScale=LogScale,
                                              FitDoubleExponential=FitDoubleExponential,
                                              PlotNumber=PlotNumber,
                                              MinPlotInd=MinPlotInd,
                                              MaxPlotInd=MaxPlotInd,
                                              PlotIndex=PlotIndex,
                                              T2MaxTime=T2MaxTime,
                                              chi_kappa_f=[chi, kappa, f])

    with open('FittingOutput', 'wb') as fp:
        pickle.dump(OutputList, fp)
