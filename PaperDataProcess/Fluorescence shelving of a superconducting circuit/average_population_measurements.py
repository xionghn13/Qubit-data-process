from QubitDataProcessPackages import *
from QubitDecayFunc import T1_curve, rabi_curve, AutoRotate, DoubleExp_curve, TwoExp_curve
import os


def plot_average_population(DataPath, RabiFileList, BackgroundFolder='', BackgroundFile='', CircleCorrection=False,
                            CorrectionParam=[1, 0.027, 0.798, -0.0098], PopulationConversionConstList=[[1, 1.390]],
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

    num_file = len(RabiFileList)
    P0_list = []
    P1_list = []
    P2_list = []
    for i, RabiFile in enumerate(RabiFileList):
        # read data file
        InterleavedMeasurement = 'interleaved' in RabiFile.split('_')
        ContrastMeasurement = 'contrast' in RabiFile.split('_')
        ReadoutPower = edf.readReadoutPowerLabber(DataPath + RabiFile)
        ReadoutFreq = edf.readReadoutFreqLabber(DataPath + RabiFile)
        if RabiFile.startswith('transient'):
            MeasurementType = 'transient no ref'
            if InterleavedMeasurement:
                [Time, Complex] = edf.readMultiRabiInterleavedLabber(DataPath + RabiFile)
            else:
                [Time, Complex] = edf.readMultiRabiLabber(DataPath + RabiFile)
        elif RabiFile.startswith('t1'):
            MeasurementType = 't1 no ref'
            if InterleavedMeasurement:
                [Time, Complex] = edf.readMultiT1InterleavedLabber(DataPath + RabiFile)
            else:
                [Time, Complex] = edf.readMultiT1Labber(DataPath + RabiFile)
        if ContrastMeasurement:
            num_t = Complex.shape[0] // 3
            Complex = np.concatenate((Complex[:, 0].reshape((num_t, 3)), Complex[:, 1].reshape((num_t, 3))), axis=1)
        if LimitTimeRange:
            TimeInd = (EndTime >= Time) == (Time >= StartTime)
            Time = Time[TimeInd]
            if len(Complex.shape) == 1:
                Complex = Complex[TimeInd]
            elif len(Complex.shape) == 2:
                Complex = Complex[TimeInd, :]
            elif len(Complex.shape) == 3:
                Complex = Complex[TimeInd, :, :]

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
        PopulationConversionConst = PopulationConversionConstList[i]
        Population = (PopulationConversionConst[0] - RComplex.real) * PopulationConversionConst[1]
        num_time = len(Time)
        if i == 0:
            population_matrix = np.zeros((num_time, 3, num_file))
        population_matrix[:, :, i] = Population
    plot_population(Time, population_matrix)
    if ShowFig:
        plt.show()
    else:
        plt.close('all')

    return [Time, population_matrix]


def plot_population(Time, population_matrix):
    num_file = population_matrix.shape[2]
    # print(cov)
    if ShowFig:
        for i in range(3):
            fig, ax = plt.subplots()
            ax.grid(linestyle='--')
            for j in range(num_file):
                plt.plot(Time, population_matrix[:, i, j])
            plt.xlabel('Time(ns)', fontsize='x-large')
            plt.ylabel('P' + str(i), fontsize='x-large')
            plt.tick_params(axis='both', which='major', labelsize='x-large')
            plt.tight_layout()

        P_sum_matrix = np.sum(population_matrix, axis=1)
        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        for j in range(num_file):
            plt.plot(Time, P_sum_matrix[:, j])
        plt.xlabel('Time(ns)', fontsize='x-large')
        plt.ylabel('P_sum', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()

        fig, ax = plt.subplots()
        ax.grid(linestyle='--')
        avg_population_matrix = np.average(population_matrix, axis=2)
        std_population_matrix = np.std(population_matrix, axis=2)
        avg_P_sum_matrix = np.average(P_sum_matrix, axis=1)
        std_P_sum_matrix = np.std(P_sum_matrix, axis=1)
        for i in range(3):
            ax.errorbar(Time, avg_population_matrix[:, i], yerr=std_population_matrix[:, i], fmt='o',
                        label='P' + str(i))
        ax.errorbar(Time, avg_P_sum_matrix, yerr=std_P_sum_matrix, fmt='o', label='P_sum')
        plt.legend()
        plt.xlabel('Time(ns)', fontsize='x-large')
        plt.ylabel('Population', fontsize='x-large')
        plt.tick_params(axis='both', which='major', labelsize='x-large')
        plt.tight_layout()


if __name__ == '__main__':
    # DataFolderName = '11112019_back to waveguide'
    DataPath = 'C:\SC Lab\GitHubRepositories\Qubit-data-process\PaperDataProcess\Fluorescence shelving of a superconducting circuit\Fluorescence/'
    # DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2020/01\Data_0130\\'
    BackgroundFolder = 'C:\SC Lab\Projects\Fluxonium\data_process/ziggy4/'
    BackgroundFile = []
    Plus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    Minus50MHzBackgroundFile = 'one_tone_4.05GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_020419.dat'
    BackgroundFile = 'power spectroscopy_116.hdf5'
    TransientFileList = [
        # 'transient_P2_P1_interleaved_2020-02-07-21-39-28.hdf5',
        # 'transient_P2_P1_interleaved_2020-02-08-07-43-27.hdf5',
        # 'transient_P2_P1_interleaved_2020-02-09-02-28-01.hdf5',
        'transient_P2_P1_interleaved_2020-02-09-12-31-53.hdf5',
        'transient_P2_P1_interleaved_2020-02-09-23-06-01.hdf5',
        'transient_P2_P1_interleaved_2020-02-10-09-40-10.hdf5',
        'transient_P2_P1_interleaved_2020-02-10-20-14-23.hdf5',
    ]

    T1FileList = [
        # 't1_P2_P1_interleaved_2020-02-08-02-41-26.hdf5',
        # 't1_P2_P1_interleaved_2020-02-08-12-45-25.hdf5',
        # 't1_P2_P1_interleaved_2020-02-09-07-29-51.hdf5',
        't1_P2_P1_interleaved_2020-02-09-17-33-47.hdf5',
        't1_P2_P1_interleaved_2020-02-10-04-07-55.hdf5',
        't1_P2_P1_interleaved_2020-02-10-14-42-07.hdf5',
        't1_P2_P1_interleaved_2020-02-11-01-16-20.hdf5',
    ]
    PopulationConversionConstList = [
        # [1., 0.9507761023676382],
        # [1., 0.9507761023676382],
        # [1., 0.9684277144489823],
        [1., 0.9684277144489823],
        [1., 0.9684277144489823],
        [1., 0.9684277144489823],
        [1., 0.9684277144489823],
    ]

    # RabiFile = 'transient_P2_P1_interleaved_2020-02-08-17-47-23.hdf5'
    # RabiFile = 'transient_P2_P1_interleaved_2020-02-11-06-48-32.hdf5'

    IQModFreq = 0.05
    CircleCorrection = False
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

    [transient_time, transient_population_matrix] = plot_average_population(DataPath, RabiFileList=TransientFileList,
                                                                            BackgroundFolder=BackgroundFolder,
                                                                            BackgroundFile=BackgroundFile,
                                                                            IQModFreq=IQModFreq,
                                                                            PopulationConversionConstList=PopulationConversionConstList,
                                                                            FitTwoExponential=FitTwoExponential,
                                                                            CircleCorrection=CircleCorrection,
                                                                            CorrectionParam=CorrectionParam,
                                                                            PhaseSlope=PhaseSlope,
                                                                            PhaseReferenceFreq=PhaseReferenceFreq,
                                                                            Calibration=Calibration,
                                                                            FitCorrectedR=FitCorrectedR,
                                                                            LimitTimeRange=LimitTimeRange,
                                                                            RotateComplex=RotateComplex,
                                                                            StartTime=StartTime, EndTime=EndTime,
                                                                            FitDoubleExponential=FitDoubleExponential,
                                                                            SaveFig=SaveFig, ShowFig=ShowFig,
                                                                            LogScale=LogScale)
    [t1_time, t1_population_matrix] = plot_average_population(DataPath, RabiFileList=T1FileList,
                                                              BackgroundFolder=BackgroundFolder,
                                                              BackgroundFile=BackgroundFile,
                                                              IQModFreq=IQModFreq,
                                                              PopulationConversionConstList=PopulationConversionConstList,
                                                              FitTwoExponential=FitTwoExponential,
                                                              CircleCorrection=CircleCorrection,
                                                              CorrectionParam=CorrectionParam,
                                                              PhaseSlope=PhaseSlope,
                                                              PhaseReferenceFreq=PhaseReferenceFreq,
                                                              Calibration=Calibration,
                                                              FitCorrectedR=FitCorrectedR,
                                                              LimitTimeRange=LimitTimeRange,
                                                              RotateComplex=RotateComplex,
                                                              StartTime=StartTime, EndTime=EndTime,
                                                              FitDoubleExponential=FitDoubleExponential,
                                                              SaveFig=SaveFig, ShowFig=ShowFig, LogScale=LogScale)

    f = h5py.File(DataPath + 'average_population_measurements.hdf5', 'w')
    # timeset0 = f.create_dataset('time0', data=TimeList[:1])
    # Yset0 = f.create_dataset('y0', data=y_dataList[:1])
    a = f.create_dataset('transient_time', data=transient_time)
    b = f.create_dataset('transient_population_matrix', data=transient_population_matrix)

    c = f.create_dataset('t1_time', data=t1_time)
    d = f.create_dataset('t1_population_matrix', data=t1_population_matrix)
    f.close()
