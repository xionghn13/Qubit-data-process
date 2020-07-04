import numpy as np
import Labber
from matplotlib import pyplot as plt
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
from scipy.optimize import curve_fit
import DataManipulationFunc as dmf
import FunctionLib as fl


def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


# constants
kB = 1.38e-23
h = 6.626e-34



def histogram_preselected_single_parameter_sweep(sweep_quantity, signal, gg_estimate=[0, 0], num_blob=4):
    rabi_signal = np.zeros(len(sweep_quantity), dtype=complex)
    rabi_signal_preselected = np.zeros((len(sweep_quantity), num_blob), dtype=complex)

    herald_signal = signal[0, 0::2] * 1e6
    select_signal = signal[0, 1::2] * 1e6

    X, Y, H = ssf.complex_array_to_2d_histogram(herald_signal)

    params, params_err = ssf.full_blob_analysis(herald_signal, num_blob=num_blob)
    print('-----')
    print(params)

    heights_fit, centers_fit, sigmas_fit = ssf.unwrap_blob_parameters(params)
    gg_ind = ssf.closest_blob_index(centers_fit, gg_estimate)

    print(gg_ind)
    fit = fg.multi_gaussian(X, Y, params)

    fig, ax = plt.subplots()
    plt.pcolormesh(X, Y, H, cmap='GnBu')
    plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
    plt.contour(X, Y, fit, cmap=plt.cm.copper)

    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    for ind_sweep_quantity in range(len(sweep_quantity)):
        herald_signal = signal[ind_sweep_quantity, 0::2] * 1e6
        select_signal = signal[ind_sweep_quantity, 1::2] * 1e6
        sReal = np.real(herald_signal)
        sImag = np.imag(herald_signal)
        if ind_sweep_quantity == 0:
            X, Y, H = ssf.complex_array_to_2d_histogram(herald_signal)
            plt.pcolormesh(X, Y, H, cmap='GnBu')
            # plt.colorbar()
        rabi_signal[ind_sweep_quantity] = np.average(select_signal)
        for ind_blob in range(num_blob):
            ind = ((sReal - centers_fit[ind_blob, 0]) / sigmas_fit[ind_blob, 0] / width_threshold) ** 2 + (
                    (sImag - centers_fit[ind_blob, 1]) / sigmas_fit[ind_blob, 1] / width_threshold) ** 2 < 1
            rabi_signal_preselected[ind_sweep_quantity, ind_blob] = np.average(select_signal[ind])


    plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
    for ind_blob in range(num_blob):
        plt.plot(np.real(rabi_signal_preselected[:, ind_blob]), np.imag(rabi_signal_preselected[:, ind_blob]),
                 label='preselect blob ' + str(ind_blob))
    plt.legend()
    plt.xlabel('I (uV)')
    plt.ylabel('Q (uV)')

    freq_guess = 4 / sweep_quantity[-1]
    V_complex = dmf.rotate(rabi_signal_preselected[:, gg_ind])
    V_real = V_complex.imag
    fig, ax = plt.subplots()
    plt.plot(sweep_quantity, V_real)
    A_guess = V_real[0] - V_real[-1]
    B_guess = (V_real[0] + V_real[-1]) / 2
    T1_guess = sweep_quantity[-1]
    Tpi_guess = sweep_quantity[-1] / 10
    phi0_guess = 0
    try:
        guess = [A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess]
        opt, cov = curve_fit(fl.rabi_curve, ydata=V_real, xdata=sweep_quantity, p0=guess)
        err = np.sqrt(np.diag(cov))
        axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
        plt.plot(axis_nice, fl.rabi_curve(axis_nice, *opt))
        plt.title('half period: %.5G\u00B1%.4G, decay constant: %.5G\u00B1%.4G' % (opt[3], err[3], opt[1], err[1]))
    except RuntimeError:
        guess = [A_guess, T1_guess / 5, V_real[-1]]
        print(guess)
        try:
            opt, cov = curve_fit(fl.T1_curve, ydata=V_real, xdata=sweep_quantity, p0=guess)
            err = np.sqrt(np.diag(cov))
        except RuntimeError:
            opt = guess
            err = np.zeros_like(opt)
        axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
        plt.plot(axis_nice, fl.T1_curve(axis_nice, *opt))
        plt.title('decay constant: %.5G\u00B1%.4G' % (opt[1], err[1]))
        ax.grid(linestyle='--')

    #
    plt.show()


############################################################

if __name__ == '__main__':
    file_path = 'Z:\Projects\Fluxonium\Data\Augustus 18\\2020\\03\Data_0312\\'
    file_name = 'rabi_A_heralding_2.hdf5'
    f = Labber.LogFile(file_path + file_name)

    gg_estimate = [-100, 200]

    num_blob = 4
    width_threshold = 2
    measurement_type = file_name.split('_')[0]

    sweep_quantity_dict = {
        'rabi': 'Multi-Qubit Pulse Generator - Amplitude #1',
        # 'rabi': 'Multi-Qubit Pulse Generator - Amplitude, 2QB #12',
        # 'Rabi': 'Multi-Qubit Pulse Generator - Width',
        # 'Ramsey': 'Multi-Qubit Pulse Generator - Pulse spacing',
        'Ramsey': 'Multi-Qubit Pulse Generator - Parameter #1',
        'T1': 'Multi-Qubit Pulse Generator - Delay after heralding',

    }
    sweep_quantity_name = sweep_quantity_dict[measurement_type]
    signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')
    sweep_quantity = f.getData(sweep_quantity_name)[0]
    if measurement_type == 'Ramsey':
        sweep_quantity *= 1e6
    elif sweep_quantity_name == 'Multi-Qubit Pulse Generator - Width':
        sweep_quantity *= 1e6
    elif sweep_quantity_name == 'Multi-Qubit Pulse Generator - Delay after heralding':
        sweep_quantity *= 1e6
    # print(sweep_quantity)
    print(signal.shape)
    histogram_preselected_single_parameter_sweep(sweep_quantity, signal, gg_estimate=gg_estimate)
