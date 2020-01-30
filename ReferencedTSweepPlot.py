import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, rabi_curve, AutoRotate, DoubleExp_curve, rabi_two_exp_curve
import ExtractDataFunc as edf
import os


def plotReferencedTSweep(DataPath, RabiFile, BackgroundFolder='', BackgroundFile='', Plus50MHzBackgroundFile='',
                         Minus50MHzBackgroundFile='', CircleCorrection=False,
                         CorrectionParam=[1, 0.027, 0.798, -0.0098], PopulationConversionConst=[1, 1.390],
                         IQModFreq=0.05, PhaseSlope=326.7041108065019, PhaseReferenceFreq=4.105, Calibration=False,
                         FitCorrectedR=False, LimitTimeRange=False, RotateComplex=True, FitDoubleExponential=False,
                         StartTime=5000, EndTime=1e8, SaveFig=True, ShowFig=False, LogScale=False):
    # read calibration file
    if Calibration:
        if BackgroundFile == []:
            [Plus50MHzBackFreq, Plus50MHzBackComplex] = edf.readFSweepDat(DataPath + Plus50MHzBackgroundFile)
            Plus50MHzBackPowerStr = Plus50MHzBackgroundFile.split('_')[5][:-3]
            Plus50MHzBackPower = float(Plus50MHzBackPowerStr)

            [Minus50MHzBackFreq, Minus50MHzBackComplex] = edf.readFSweepDat(DataPath + Minus50MHzBackgroundFile)
            Minus50MHzBackPowerStr = Minus50MHzBackgroundFile.split('_')[5][:-3]
            Minus50MHzBackPower = float(Minus50MHzBackPowerStr)
        elif BackgroundFile.startswith('one_tone'):
            [Plus50MHzBackFreq, Plus50MHzBackComplex] = edf.readFSweepDat(BackgroundFolder + BackgroundFile)
            Plus50MHzBackPowerStr = BackgroundFile.split('_')[5][:-3]
            Plus50MHzBackPower = float(Plus50MHzBackPowerStr)

            [Minus50MHzBackFreq, Minus50MHzBackComplex] = [Plus50MHzBackFreq, Plus50MHzBackComplex]
            Minus50MHzBackPowerStr = Plus50MHzBackPowerStr
            Minus50MHzBackPower = Plus50MHzBackPower
        elif BackgroundFile.endswith('dat'):
            BackgroundFileStrList = BackgroundFile.split('_')
            MeasurementType = BackgroundFileStrList[1]
            if MeasurementType == 'rabi' and BackgroundFileStrList[3] != 'no pump':
                MeasurementType = 'pump rabi'
            BackFreqDict = {'t1': 3, 'rabi': 5, 'pump rabi': 6}
            BackFreqInd = BackFreqDict[MeasurementType]
            BackPowerDict = {'t1': 5, 'rabi': 7, 'pump rabi': 8}
            BackPowerInd = BackPowerDict[MeasurementType]
            BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
            BackLowerFreq = round(float(BackgroundFileStrList[BackFreqInd][:-3]) - IQModFreq, 4)
            BackHigherFreq = round(BackLowerFreq + IQModFreq * 2, 15)
            BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
            Plus50MHzBackPower = BackPower
            Minus50MHzBackPower = BackPower
            if MeasurementType in ('rabi', 'transient'):
                [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(BackgroundFolder + BackgroundFile)
            elif MeasurementType == 't1':
                [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(BackgroundFolder + BackgroundFile)
            Plus50MHzBackFreq = np.array([BackLowerFreq, BackHigherFreq])
            Plus50MHzBackComplex = np.array([ComplexLowerFreq.mean(), ComplexHigherFreq.mean()])
            Minus50MHzBackFreq = Plus50MHzBackFreq
            Minus50MHzBackComplex = Plus50MHzBackComplex
        elif BackgroundFile.endswith('hdf5'):
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
    if RabiFile.endswith('dat'):
        RabiFileStrList = RabiFile.split('_')
        MeasurementType = RabiFileStrList[1]
        if MeasurementType == 'rabi' and RabiFileStrList[3] != 'no pump':
            MeasurementType = 'pump rabi'
        if MeasurementType == 'Ch1':
            MeasurementType = 'Ch1 rabi'
        if MeasurementType == 'Ch2 rabi 2Vpp':
            MeasurementType = 'Ch1 pump rabi'

        ReadoutFreqDict = {'t1': 3, 'rabi': 5, 'Ch1 rabi': 6, 'pump rabi': 6, 'Ch1 pump rabi': 4, 't2': 3, 't2echo': 3}
        ReadoutFreqInd = ReadoutFreqDict[MeasurementType]
        ReadoutPowerDict = {'t1': 5, 'rabi': 7, 'Ch1 rabi': 8, 'pump rabi': 8, 'Ch1 pump rabi': 6, 't2': 5, 't2echo': 5}
        ReadoutPowerInd = ReadoutPowerDict[MeasurementType]
        ReadoutFreq = round(float(RabiFileStrList[ReadoutFreqInd][:-3]), 4)
        ReadoutLowerFreq = round(float(RabiFileStrList[ReadoutFreqInd][:-3]) - IQModFreq, 4)
        ReadoutHigherFreq = round(ReadoutLowerFreq + IQModFreq * 2, 15)
        ReadoutPower = float(RabiFileStrList[ReadoutPowerInd][:-3])
        QubitFreqDict = {'rabi': 8, 'Ch1 rabi': 9, 't1': 6, 'Ch1 pump rabi': 7, 't2': 6, 't2echo': 7}
        QubitFreqInd = QubitFreqDict[MeasurementType]
        QubitFreq = float(RabiFileStrList[QubitFreqInd][5:-3])
        if QubitFreq in (ReadoutLowerFreq, ReadoutHigherFreq) and MeasurementType == 'rabi':
            MeasurementType = 'transient'

        if MeasurementType in ('rabi', 'Ch1 rabi', 'transient', 'Ch1 pump rabi'):
            [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(DataPath + RabiFile)
            if LimitTimeRange:
                TimeInd = (EndTime >= Time) == (Time >= StartTime)
                Time = Time[TimeInd]
                ComplexLowerFreq = ComplexLowerFreq[TimeInd]
                ComplexHigherFreq = ComplexHigherFreq[TimeInd]
        elif MeasurementType == 't1':
            [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(DataPath + RabiFile)
            if LimitTimeRange:
                TimeInd = (EndTime >= Time) == (Time >= StartTime)
                Time = Time[TimeInd]
                ComplexLowerFreq = ComplexLowerFreq[TimeInd]
                ComplexHigherFreq = ComplexHigherFreq[TimeInd]
        elif MeasurementType in ('t2', 't2echo'):
            [Time, Complex] = edf.readT2H5(DataPath + RabiFile)
            if LimitTimeRange:
                TimeInd = (EndTime >= Time) == (Time >= StartTime)
                Time = Time[TimeInd]
                Complex = Complex[TimeInd]
    elif RabiFile.endswith('hdf5'):
        ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
        ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
        if RabiFile.startswith('transient'):
            MeasurementType = 'transient no ref'
            [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
        elif RabiFile.startswith('rabi'):
            MeasurementType = 'rabi no ref'
            if RabiFile.startswith('rabi CH1 drive'):
                [Time, Complex] = edf.readRabiCH1DriveLabber(DataPath + RabiFile)
                # [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
            elif RabiFile.startswith('rabi CH1 pumped'):
                [Time, Complex] = edf.readRabiCH1PumpedLabber(DataPath + RabiFile)
            else:
                [Time, Complex] = edf.readRabiLabber(DataPath + RabiFile)
                # print(Complex)
        elif RabiFile.startswith('int_t1'):
            MeasurementType = 't1 no ref'
            [Time, Complex] = edf.readIntegratedT1Labber(DataPath + RabiFile)
        elif RabiFile.startswith('t1'):
            MeasurementType = 't1 no ref'
            [Time, Complex] = edf.readT1Labber(DataPath + RabiFile)
        elif RabiFile.startswith('t2_ramsey'):
            MeasurementType = 't2'
            [Time, Complex] = edf.readT2Labber(DataPath + RabiFile)
        elif RabiFile.startswith('t2_echo'):
            MeasurementType = 't2echo'
            [Time, Complex] = edf.readT2Labber(DataPath + RabiFile)
        elif RabiFile.startswith('RefRabi'):
            MeasurementType = 'rabi'
            [Time, ReadoutLowerFreq, ReadoutHigherFreq, ComplexLowerFreq, ComplexHigherFreq] = edf.readRefRabiLabber(
                DataPath + RabiFile)

        if LimitTimeRange:
            TimeInd = (EndTime >= Time) == (Time >= StartTime)
            Time = Time[TimeInd]
            if RabiFile.startswith('RefRabi'):
                ComplexLowerFreq = ComplexLowerFreq[TimeInd]
                ComplexHigherFreq = ComplexHigherFreq[TimeInd]
            else:
                Complex = Complex[TimeInd]

    if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
        ComplexNormalized = Complex * 10 ** (- ReadoutPower / 20)
        if Calibration:
            RComplex = sbf.FPSweepBackgroundCalibrate(ReadoutFreq, ReadoutPower, Complex, Plus50MHzBackFreq,
                                                      Plus50MHzBackComplex, Plus50MHzBackPower)
        else:
            RComplex = ComplexNormalized

        if RotateComplex:
            RComplex = AutoRotate(RComplex)
        y_data = np.array(RComplex.real, dtype='float64')
        # y_data = np.array(RComplex.imag, dtype='float64')

    else:
        RComplexLowerFreq = sbf.FPSweepBackgroundCalibrate(ReadoutLowerFreq, ReadoutPower, ComplexLowerFreq,
                                                           Plus50MHzBackFreq,
                                                           Plus50MHzBackComplex, Plus50MHzBackPower)
        RComplexHigherFreq = sbf.FPSweepBackgroundCalibrate(ReadoutHigherFreq, ReadoutPower, ComplexHigherFreq,
                                                            Minus50MHzBackFreq,
                                                            Minus50MHzBackComplex, Minus50MHzBackPower)
        ComplexLowerFreqNormalized = ComplexLowerFreq * 10 ** (- ReadoutPower / 20)
        ComplexHigherFreqNormalized = ComplexHigherFreq * 10 ** (- ReadoutPower / 20)

        [amp_cor_re_fit, amp_cor_im_fit, P0_fit, P0_im_fit] = CorrectionParam
        if CircleCorrection:
            RComplexLowerFreq /= (amp_cor_re_fit + amp_cor_im_fit * 1j)
            RComplexLowerFreq = (RComplexLowerFreq - 1) / (P0_fit + P0_im_fit * 1j) * np.abs(
                P0_fit + P0_im_fit * 1j) + 1

        # y_data = np.array(RComplexLowerFreq.real, dtype='float64')
        if FitCorrectedR:
            y_data = np.array(np.real(RComplexLowerFreq / RComplexHigherFreq), dtype='float64')
        else:
            y_data = np.array(RComplexLowerFreq.real, dtype='float64')

    x_data = np.array(Time, dtype='float64')
    if MeasurementType in ('t1', 'transient', 't2echo', 'transient no ref', 't1 no ref'):
        B_guess = y_data[-1]
        A_guess = y_data[0].real - B_guess
        T1_guess = x_data[-1] / 10
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
        else:
            opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
            A_fit, T1_fit, B_fit = opt
            A_std, T1_std, B_std = np.sqrt(cov.diagonal())

            TimeFit = np.linspace(Time.min(), Time.max(), 200)

            FitR = T1_curve(TimeFit, A_fit, T1_fit, B_fit)
            y_pred = T1_curve(Time, A_fit, T1_fit, B_fit)
            y_guess = T1_curve(Time, A_guess, T1_guess, B_guess)
            ParamList = ['A', 'Decay time/ns', 'B']

    elif MeasurementType in ('rabi', 'Ch1 rabi', 'Ch1 pump rabi', 't2', 'rabi no ref'):
        B_guess = y_data.mean()
        A_guess = y_data[0] - B_guess
        T1_guess = x_data[-1]
        MaxInd = y_data.argmax()
        MinInd = y_data.argmin()
        FitTwoExpRabi = False
        if MeasurementType == 't2':
            Tpi_guess = np.abs(x_data[MaxInd] - x_data[MinInd])
        else:
            Tpi_guess = T1_guess / 4
        phi0_guess = 0
        T_out_guess = T1_guess
        if FitTwoExpRabi:
            C_guess = B_guess
            guess = ([A_guess, T1_guess, B_guess * 0, Tpi_guess, phi0_guess, T_out_guess, C_guess])
            # bounds = (
            #     (-2, 1, -1, 1, - np.pi / 2),
            #     (2, np.inf, 1, np.inf, np.pi / 2)
            # )
            bounds = (
                (- np.inf, 1, - 1e-5, 1, - np.pi / 2, - np.inf, - np.inf),
                (np.inf, np.inf, 0, np.inf, np.pi / 2, np.inf, np.inf)
            )
            # print(guess)
            try:
                opt, cov = curve_fit(rabi_two_exp_curve, x_data, y_data, p0=guess, bounds=bounds)
            except RuntimeError:
                print("Error - curve_fit failed")
                opt = guess
                cov = np.zeros([len(opt), len(opt)])
            A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit, T_out_fit, C_fit = opt
            A_std, T1_std, B_std, Tpi_std, phi0_std, T_out_std, C_std = np.sqrt(cov.diagonal())
            TimeFit = np.linspace(Time.min(), Time.max(), 200)
            FitR = rabi_two_exp_curve(TimeFit, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit, T_out_fit, C_fit)
            FitUpper = A_fit * np.exp(- TimeFit / T1_fit) + B_fit * np.exp(-TimeFit / T_out_fit) + C_fit
            FitLower = -A_fit * np.exp(- TimeFit / T1_fit) + B_fit * np.exp(-TimeFit / T_out_fit) + C_fit
            ParamList = ['A', 'Decay time(ns)', 'B', 'Tpi(ns)', 'pho0', 'T_out(ns)', 'C']
        else:
            guess = ([A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess])
            bounds = (
                (-2, 1, -1, 1, - np.pi / 2),
                (2, np.inf, 1, np.inf, np.pi / 2)
            )

            try:
                opt, cov = curve_fit(rabi_curve, x_data, y_data, p0=guess, bounds=bounds)
            except RuntimeError:
                print("Error - curve_fit failed")
                opt = guess
                cov = np.zeros([len(opt), len(opt)])
            A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit = opt
            A_std, T1_std, B_std, Tpi_std, phi0_std = np.sqrt(cov.diagonal())
            TimeFit = np.linspace(Time.min(), Time.max(), 200)
            FitR = rabi_curve(TimeFit, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
            FitUpper = A_fit * np.exp(- TimeFit / T1_fit) + B_fit
            FitLower = -A_fit * np.exp(- TimeFit / T1_fit) + B_fit
            ParamList = ['A', 'Decay time(ns)', 'B', 'Tpi(ns)', 'pho0']

    # print(cov)
    limit = 1.7
    if ShowFig:
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.plot(np.real(RComplex), np.imag(RComplex))
        else:
            plt.plot(np.real(RComplexLowerFreq), np.imag(RComplexLowerFreq))
            plt.plot(np.real(RComplexHigherFreq), np.imag(RComplexHigherFreq))
        if Calibration:
            plt.plot([-2, 2], [0, 0], '--')
            plt.plot([1], [0], 'ro')
        plt.xlabel('Re', fontsize='x-large')
        plt.ylabel('Im', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        if MeasurementType not in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.legend(['%.4GGHz' % ReadoutLowerFreq, '%.4GGHz' % ReadoutHigherFreq])
        plt.tight_layout()
        if Calibration:
            plt.xlim(-limit, limit)
            plt.ylim(-limit, limit)
        ax.set_aspect('equal')

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.plot(Time, np.abs(ComplexNormalized), 'o')
        else:
            plt.plot(Time, np.abs(ComplexLowerFreqNormalized), 'o')
            plt.plot(Time, np.abs(ComplexHigherFreqNormalized), 'o')
        plt.xlabel('Time/ns', fontsize='x-large')
        plt.ylabel('Abs', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        if MeasurementType not in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.legend(['%.4GGHz' % ReadoutLowerFreq, '%.4GGHz' % ReadoutHigherFreq])
        plt.tight_layout()

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.plot(Time, (np.angle(ComplexNormalized) - PhaseSlope * (ReadoutFreq - 4.105)) % (2 * np.pi), 'o')
        else:
            plt.plot(Time,
                     (np.angle(ComplexLowerFreqNormalized) - PhaseSlope * (ReadoutLowerFreq - 4.105)) % (2 * np.pi),
                     'o')
            plt.plot(Time,
                     (np.angle(ComplexHigherFreqNormalized) - PhaseSlope * (ReadoutHigherFreq - 4.105)) % (2 * np.pi),
                     'o')
        plt.xlabel('Time/ns', fontsize='x-large')
        plt.ylabel('Phase', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        if MeasurementType not in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            plt.legend(['%.4GGHz' % ReadoutLowerFreq, '%.4GGHz' % ReadoutHigherFreq])
        plt.tight_layout()

        if MeasurementType not in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plt.plot(Time, np.real(RComplexLowerFreq / RComplexHigherFreq), 'o')
            if FitCorrectedR:
                plt.plot(TimeFit, FitR)
                if MeasurementType == 't1':
                    plt.title('T1=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                        T1_fit / 1000, T1_std / 1000, A_fit, B_fit))
                elif MeasurementType == 'transient':
                    plt.title('T_transient=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                        T1_fit / 1000, T1_std / 1000, A_fit, B_fit))
                elif MeasurementType in ('rabi', 'Ch1 rabi', 'Ch1 pump rabi'):
                    if FitTwoExpRabi:
                        plt.title('Tpi=%.3Gus, T1=%.3Gus, T_out=%.3Gus$\pm$%.2Gus, A=%.3G, B=%.3G, phi0=%.3G' % (
                            Tpi_fit / 1000, T1_fit / 1000, T_out_fit / 1000, T_out_std / 1000, A_fit, B_fit, phi0_fit))
                    else:
                        plt.title('Tpi=%.3Gus, T1=%.3Gus, A=%.3G, B=%.3G, phi0=%.3G' % (
                            Tpi_fit / 1000, T1_fit / 1000, A_fit, B_fit, phi0_fit))
            plt.xlabel('Time/ns', fontsize='x-large')
            plt.ylabel('Re', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        if MeasurementType in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref'):
            if LogScale:
                plt.plot(Time, y_data - B_fit, 'o')
            else:
                plt.plot(Time, y_data, 'o')
        else:
            plt.plot(Time, np.real(RComplexLowerFreq), 'o')
            if FitCorrectedR:
                plt.plot(Time, np.real(RComplexHigherFreq), 'o')
        if not FitCorrectedR:
            if LogScale:
                plt.plot(TimeFit, FitR - B_fit)
            else:
                plt.plot(TimeFit, FitR)
            if MeasurementType in ('t1', 't1 no ref'):
                if FitDoubleExponential:
                    plt.title('TR=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                              'Tqp=%.3G$\pm$%.2Gus, lambda=%.3G$\pm$%.2G' % (
                                  TR_fit / 1000, TR_std / 1000, A_fit, B_fit, Tqp_fit / 1000, Tqp_std / 1000, lamb_fit,
                                  lamb_std))
                else:
                    plt.plot(Time, y_guess, ':')
                    plt.title('T1=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                        T1_fit / 1000, T1_std / 1000, A_fit, B_fit))
            elif MeasurementType in ('transient', 'transient no ref'):
                if FitDoubleExponential:
                    plt.title('TR=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G\n'
                              'Tqp=%.3G$\pm$%.2Gus, lambda=%.3G$\pm$%.2G' % (
                                  TR_fit / 1000, TR_std / 1000, A_fit, B_fit, Tqp_fit / 1000, Tqp_std / 1000, lamb_fit,
                                  lamb_std))
                else:
                    plt.title('T_transient=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                        T1_fit / 1000, T1_std / 1000, A_fit, B_fit))
            elif MeasurementType in ('rabi', 'Ch1 rabi', 'Ch1 pump rabi', 'rabi no ref'):
                plt.plot(TimeFit, FitUpper, ':')
                plt.plot(TimeFit, FitLower, ':')
                if FitTwoExpRabi:
                    plt.title('Tpi=%.3Gus, T1=%.3Gus, T_out=%.3Gus$\pm$%.2Gus\n A=%.3G, B=%.3G, phi0=%.3G, C=%.3G' % (
                        Tpi_fit / 1000, T1_fit / 1000, T_out_fit / 1000, T_out_std / 1000, A_fit, B_fit, phi0_fit, C_fit))
                else:
                    plt.title('Tpi=%.3Gus, T1=%.3Gus \n A=%.3G, B=%.3G, phi0=%.3G' % (
                        Tpi_fit / 1000, T1_fit / 1000, A_fit, B_fit, phi0_fit))
        plt.xlabel('Time/ns', fontsize='x-large')
        plt.ylabel('Re', fontsize='x-large')
        if MeasurementType == 't2':
            plt.title('Tpi=%.3Gus, T2Ramsey=%.3Gus$\pm$%.2Gus\n'
                      'A=%.3G, B=%.3G, phi0=%.3G' % (
                          Tpi_fit / 1000, T1_fit / 1000, T1_std / 1000, A_fit, B_fit, phi0_fit))
        elif MeasurementType == 't2echo':
            plt.title('T2echo=%.3G$\pm$%.2Gus, A=%.3G, B=%.3G' % (
                T1_fit / 1000, T1_std / 1000, A_fit, B_fit))
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()
        if MeasurementType not in ('t2', 't2echo', 'transient no ref', 't1 no ref', 'rabi no ref') and FitCorrectedR:
            plt.legend(['%.4GGHz' % ReadoutLowerFreq, '%.4GGHz' % ReadoutHigherFreq])
        if LogScale:
            ax.set_yscale('log')
        if SaveFig:
            FigPath = DataPath + 'figures/'
            if not os.path.exists(FigPath):
                os.makedirs(FigPath)
            FigName = RabiFile.split('.')[0] + '.PNG'
            plt.savefig(FigPath + FigName)

        if Calibration:
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plt.plot(Time, (PopulationConversionConst[0] - y_data) * PopulationConversionConst[1], 'o')
            if not FitCorrectedR:
                plt.plot(TimeFit, (PopulationConversionConst[0] - FitR) * PopulationConversionConst[1])
            if MeasurementType in ('rabi', 'Ch1 rabi', 'Ch1 pump rabi', 'rabi no ref'):
                plt.plot(TimeFit, (PopulationConversionConst[0] - FitLower) * PopulationConversionConst[1], ':')
                plt.plot(TimeFit, (PopulationConversionConst[0] - FitUpper) * PopulationConversionConst[1], ':')
            plt.xlabel('Time(ns)', fontsize='x-large')
            plt.ylabel('Population', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        if MeasurementType in ('t1', 't1 no ref'):
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            plt.plot(Time, y_data - y_pred, 'o--')
            plt.xlabel('Time/ns', fontsize='x-large')
            plt.ylabel('Residual', fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        if ShowFig:
            plt.show()
        else:
            plt.close('all')

    return {'opt': opt, 'cov': cov, 'ParamList': ParamList}


if __name__ == '__main__':
    DataFolderName = '11112019_back to waveguide'
    DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2020/01\Data_0130\\'
    # DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2019/11\Data_1123\\'
    BackgroundFolder = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
    BackgroundFile = []
    # BackgroundFile = '021219_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-35dBm_0.8_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi300_duty50000readout3us.h5'
    # Plus50MHzBackgroundFile = '012819_rabi_CH2(AWG1Vpp)_no pump_readout_4.146GHz__-20dBm_qubit4.096GHz_-25dBm_4.9_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi100000_duty150000readout3us.h5'
    Plus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    Minus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    BackgroundFile = 'power spectroscopy_116.hdf5'
    RabiFile = 't2_ramsey_Drive1_2.hdf5'
    IQModFreq = 0.05
    CircleCorrection = False
    CorrectionParam = [1.1, 0.044, 0.737, 0.037]
    PhaseSlope = 326.7041108065019
    PhaseReferenceFreq = 4.105
    Calibration = True
    FitCorrectedR = False
    LimitTimeRange = False
    RotateComplex = False
    FitDoubleExponential = False
    LogScale = False
    SaveFig = False
    ShowFig = True
    StartTime = 0
    EndTime = 30e3
    PopulationConversionConst = [1, 1. / 1]
    FitDict = plotReferencedTSweep(DataPath, RabiFile, BackgroundFolder=BackgroundFolder, BackgroundFile=BackgroundFile,
                                   Plus50MHzBackgroundFile=Plus50MHzBackgroundFile,
                                   PopulationConversionConst=PopulationConversionConst,
                                   Minus50MHzBackgroundFile=Minus50MHzBackgroundFile,
                                   IQModFreq=IQModFreq,
                                   CircleCorrection=CircleCorrection,
                                   CorrectionParam=CorrectionParam,
                                   PhaseSlope=PhaseSlope, PhaseReferenceFreq=PhaseReferenceFreq,
                                   Calibration=Calibration,
                                   FitCorrectedR=FitCorrectedR, LimitTimeRange=LimitTimeRange,
                                   RotateComplex=RotateComplex,
                                   StartTime=StartTime, EndTime=EndTime, FitDoubleExponential=FitDoubleExponential,
                                   SaveFig=SaveFig, ShowFig=ShowFig, LogScale=LogScale)
    print(FitDict)
    # print(FitDict['opt'][3])
    # print(FitDict['ParamList'][3])
