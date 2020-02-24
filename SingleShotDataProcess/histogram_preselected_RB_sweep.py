import numpy as np
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
import Labber
from matplotlib import pyplot as plt
from qutip import *
from scipy.optimize import curve_fit
import QubitDecayFunc as qdf


def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


def randomized_benchmarking_0(x, p, a, b):
    return a * p ** x + b


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0222\RB_heralded_AWG_qubitB_interleaved.hdf5')
num_blob = 4

width_threshold = 2  # sigma

p_no_interleaving = 0.9748
p_err_no_interleaving = 0.003257

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

                params = fg.fitgaussian(X, Y, H, param_mat)
                fit = fg.multi_gaussian(X, Y, params)

                if ind_pulse_randomize + ind_pulse_type == 0:
                    print('-----')
                    print(params)
                    centers_fit = params[:, 1:3]
                    sigmas_fit = params[:, 3:]
                    heights_fit = params[:, 0]
                    most_pts_ind = np.argmax(heights_fit)
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

p_array = np.zeros(num_blob)
p_err_array = np.zeros_like(p_array)
fig, ax = plt.subplots()
ax.grid(linestyle='--')


gate_list = ['Xp', 'Xm', 'X2p', 'X2m', 'Yp', 'Ym', 'Y2p', 'Y2m', 'I']
n = 1
d = 2 ** n
for ind_pulse_type in range(len(pulse_type)):
    for ind_blob in range(num_blob):
        V_complex = qdf.AutoRotate(avg_signal[ind_blob, :, ind_pulse_type])
        toFit = np.real(V_complex)
        guess = ([0.9, np.max(toFit) - np.min(toFit), np.min(toFit)])
        opt, cov = curve_fit(randomized_benchmarking_0, ydata=toFit, xdata=pulse_num, p0=guess)
        err = (np.sqrt(np.diag(cov)))
        p_array[ind_blob] = opt[0]
        p_err_array[ind_blob] = err[0]
        if ind_pulse_type == 0:
            plt.plot(pulse_num, np.real(V_complex))
            plt.plot(pulse_num, randomized_benchmarking_0(pulse_num, *opt), label='Zeroth order fit')

    parameter = np.average(p_array)
    parameter_err = parameter * np.sqrt(np.sum((p_err_array / p_array) ** 2))
    parameter /= p_no_interleaving
    parameter_err = parameter * np.sqrt((err[0] / opt[0]) ** 2 + (p_err_no_interleaving / p_no_interleaving) ** 2)
    error = abs((d - 1) * (1 - parameter) / d)
    error_err = (d - 1) * parameter_err / d
    error = error / 1.875
    error_err = error_err / 1.875
    print('For gate %s, p_C/p=%.4G\u00B1%.4G, 0-order model fidelity %.4G\u00B1%.4G' % (gate_list[ind_pulse_type],
    parameter, parameter_err, 1 - error, error_err))

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
