import numpy as np
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
import Labber
from matplotlib import pyplot as plt
from qutip import *
from scipy.optimize import curve_fit
import FunctionLib as qdf


def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


def randomized_benchmarking_0(x, p, a, b):
    return a * p ** x + b


def randomized_benchmarking_1(x, p, a, b, c):
    return a * p ** x + c * (x - 1) * p ** (x - 2) + b


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\03\Data_0302\RB_heralded_AWG_qubitB_interleaved.hdf5')

num_blob = 4

width_threshold = 2  # sigma

gg_estimate = [100, 150]

p_no_interleaving = 0.9807
p_err_no_interleaving = 0.00145

p1_no_interleaving = 0.9583
p1_err_no_interleaving = 0.000905

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')[:, :]
pulse_num = f.getData('Multi-Qubit Pulse Generator - Number of Cliffords')[0]
pulse_randomize = np.unique(f.getData('Multi-Qubit Pulse Generator - Randomize')[:, 0])
pulse_type = np.unique(f.getData('Multi-Qubit Pulse Generator - Interleaved 1-QB Gate'))
blobs_matrix = np.zeros((len(pulse_randomize), 4, 5))

rb_signal = np.zeros((num_blob, len(pulse_num), len(pulse_randomize), len(pulse_type)), dtype=complex)
rabi_signal = np.zeros(len(pulse_num), dtype=complex)
# rabi_signal_preselected_1 = np.zeros(len(pulse_num), dtype=complex)
# rabi_signal_preselected_2 = np.zeros(len(pulse_num), dtype=complex)
# rabi_signal_preselected_3 = np.zeros(len(pulse_num), dtype=complex)
# rabi_signal_preselected_4 = np.zeros(len(pulse_num), dtype=complex)
rabi_signal_preselected = np.zeros((len(pulse_num), num_blob), dtype=complex)

for ind_pulse_type in range(len(pulse_type)):
    print('analyzing pulse', ind_pulse_type)
    for ind_pulse_randomize in range(len(pulse_randomize)):
        print('analyzing ind_pulse_randomize', ind_pulse_randomize)
        for ind_pulse_num in range(len(pulse_num)):
            # print('1')
            herald_signal = signal[ind_pulse_type * len(pulse_randomize) * len(pulse_num) + ind_pulse_randomize * len(
                pulse_num) + ind_pulse_num, 0::2] * 1e6
            select_signal = signal[ind_pulse_type * len(pulse_randomize) * len(pulse_num) + ind_pulse_randomize * len(
                pulse_num) + ind_pulse_num, 1::2] * 1e6
            sReal = np.real(herald_signal)
            sImag = np.imag(herald_signal)
            if ind_pulse_num == 0:
                H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
                X, Y = np.meshgrid(xedges[1:], yedges[1:])
                H = H.T

                centers, sigmas = ssf.getBlobCenters(sReal, sImag, num_blob)
                heights = ssf.getCenterHeights(X, Y, H, centers)
                param_mat = np.concatenate((heights.reshape(num_blob, 1), centers, sigmas), axis=1)

                params, params_err = fg.fitgaussian(X, Y, H, param_mat)
                fit = fg.multi_gaussian(X, Y, params)

                if ind_pulse_randomize + ind_pulse_type == 0:
                    print('-----')
                    print(params)
                    centers_fit = params[:, 1:3]
                    sigmas_fit = params[:, 3:]
                    heights_fit = params[:, 0]
                    most_pts_ind = np.argmax(heights_fit)
                    gg_ind = np.argmin(np.sum((centers_fit - gg_estimate) ** 2, axis=1))
                    print(gg_ind)
                    # most_pts_ind = 0
                    # print(most_pts_ind)
                    # print(centers_fit)
                    fig, ax = plt.subplots()
                    plt.pcolormesh(X, Y, H, cmap='GnBu')
                    plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
                    plt.contour(X, Y, fit, cmap=plt.cm.copper)
                elif ind_pulse_type == len(pulse_type) - 1:
                    blobs_matrix[ind_pulse_randomize, :, :] = params

            rabi_signal[ind_pulse_num] = np.average(select_signal)

            for ind_blob in range(num_blob):
                ind = ((sReal - centers_fit[ind_blob, 0]) / sigmas_fit[ind_blob, 0] / width_threshold) ** 2 + (
                        (sImag - centers_fit[ind_blob, 1]) / sigmas_fit[ind_blob, 1] / width_threshold) ** 2 < 1
                rabi_signal_preselected[ind_pulse_num, ind_blob] = np.average(select_signal[ind])
                rb_signal[ind_blob, ind_pulse_num, ind_pulse_randomize, ind_pulse_type] = rabi_signal_preselected[
                    ind_pulse_num, ind_blob]

plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
for ind_blob in range(num_blob):
    plt.plot(np.real(rb_signal[ind_blob, :, -1, -1]), np.imag(rb_signal[ind_blob, :, -1, -1]),
             label='preselect blob' + str(ind_blob))

plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

avg_signal = np.average(rb_signal, axis=2)
std_signal = np.std(rb_signal, axis=2)

p_array = np.zeros(num_blob)
p_err_array = np.zeros_like(p_array)
p1_array = np.zeros_like(p_array)
p1_err_array = np.zeros_like(p_array)


def p_to_fidelity_interleaved(p_arr, p_err_arr, dim, p, p_err):
    parameter_0 = p_arr[gg_ind]
    parameter_0_err = p_err_arr[gg_ind]
    parameter = parameter_0 / p
    parameter_err = parameter * np.sqrt((parameter_0_err / parameter_0) ** 2 + (p_err / p) ** 2)
    error = abs((dim - 1) * (1 - parameter) / dim)
    error_err = (dim - 1) * parameter_err / dim
    return [parameter, parameter_err, error, error_err]


gate_list = ['Xp', 'Xm', 'X2p', 'X2m', 'Yp', 'Ym', 'Y2p', 'Y2m', 'I']
n = 1
d = 2 ** n
for ind_pulse_type in range(len(pulse_type)):
    for ind_blob in range(num_blob):
        V_complex = qdf.AutoRotate(avg_signal[ind_blob, :, ind_pulse_type])
        toFit = np.real(V_complex)
        # 0th order
        guess = ([0.9, np.max(toFit) - np.min(toFit), np.min(toFit)])
        opt, cov = curve_fit(randomized_benchmarking_0, ydata=toFit, xdata=pulse_num, p0=guess)
        err = (np.sqrt(np.diag(cov)))
        p_array[ind_blob] = opt[0]
        p_err_array[ind_blob] = err[0]

        if ind_pulse_type == 0:
            fig = plt.figure(2)
            ax = fig.add_subplot(111)
            ax.grid(linestyle='--')
            ax.errorbar(pulse_num, np.real(V_complex), yerr=std_signal[ind_blob, :, ind_pulse_type], fmt='o')
            plt.plot(pulse_num, randomized_benchmarking_0(pulse_num, *opt), label='Zeroth order fit')

    [parameter, parameter_err, error, error_err] = p_to_fidelity_interleaved(p_array, p_err_array, d, p_no_interleaving,
                                                                             p_err_no_interleaving)
    print('For gate %s, p_C/p=%.4G\u00B1%.4G, 0-order model fidelity %.4G\u00B1%.4G' % (gate_list[ind_pulse_type],
                                                                                        parameter, parameter_err,
                                                                                        1 - error, error_err))


for ind_pulse_type in range(len(pulse_type)):
    for ind_blob in range(num_blob):
        V_complex = qdf.AutoRotate(avg_signal[ind_blob, :, ind_pulse_type])
        toFit = np.real(V_complex)

        # 1st order
        guess = ([0.9, np.max(toFit) - np.min(toFit), np.min(toFit), np.min(toFit)])
        opt, cov = curve_fit(randomized_benchmarking_1, ydata=toFit, xdata=pulse_num, p0=guess)
        err = (np.sqrt(np.diag(cov)))
        p1_array[ind_blob] = opt[0]
        p1_err_array[ind_blob] = err[0]


        if ind_pulse_type == 0:
            fig = plt.figure(3)
            ax = fig.add_subplot(111)
            ax.grid(linestyle='--')
            ax.errorbar(pulse_num, np.real(V_complex), yerr=std_signal[ind_blob, :, ind_pulse_type], fmt='o')
            plt.plot(pulse_num, randomized_benchmarking_1(pulse_num, *opt), label='Zeroth order fit')

    [parameter, parameter_err, error, error_err] = p_to_fidelity_interleaved(p1_array, p1_err_array, d, p1_no_interleaving,
                                                                             p1_err_no_interleaving)
    print('For gate %s, p_C/p=%.4G\u00B1%.4G, 1-order model fidelity %.4G\u00B1%.4G' % (gate_list[ind_pulse_type],
                                                                                        parameter, parameter_err,
                                                                                        1 - error, error_err))

fig, ax = plt.subplots()
x_blob = []
y_blob = []
for ind_pulse_randomize in range(len(pulse_randomize)):
    most_pts_ind = np.argmin(blobs_matrix[ind_pulse_randomize, :, 0])
    x_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 1]]
    y_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 2]]
plt.pcolormesh(X, Y, H, cmap='GnBu')
# plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
plt.contour(X, Y, fit, cmap=plt.cm.copper)
plt.plot(x_blob, y_blob, c='r')

plt.show()
