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
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0222\RB_heralded_AWG_qubitB.hdf5')
num_blob = 4

width_threshold = 2  # sigma

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')[:, :]
pulse_num = f.getData('Multi-Qubit Pulse Generator - Number of Cliffords')[0]
pulse_randomize = f.getData('Multi-Qubit Pulse Generator - Randomize')[:, 0]
rb_signal = np.zeros((num_blob, len(pulse_num), len(pulse_randomize)), dtype=complex)
rabi_signal = np.zeros(len(pulse_num), dtype=complex)
rabi_signal_preselected = np.zeros((len(pulse_num), num_blob), dtype=complex)
blobs_matrix = np.zeros((len(pulse_randomize), 4, 5))

for ind_pulse_randomize in range(len(pulse_randomize)):
    print('analyzing ind_pulse_randomize', ind_pulse_randomize)
    for ind_pulse_num in range(len(pulse_num)):
        herald_signal = signal[ind_pulse_randomize * len(pulse_num) + ind_pulse_num, 0::2] * 1e6
        select_signal = signal[ind_pulse_randomize * len(pulse_num) + ind_pulse_num, 1::2] * 1e6
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
            blobs_matrix[ind_pulse_randomize, :, :] = params
            if ind_pulse_randomize == 0:
                print('-----')
                print(params)
                centers_fit = params[:, 1:3]
                sigmas_fit = params[:, 3:]
                heights_fit = params[:, 0]
                most_pts_ind = np.argmax(heights_fit)
                # most_pts_ind = 0
                print(most_pts_ind)
                # print(centers_fit)

                fig, ax = plt.subplots()
                plt.pcolormesh(X, Y, H, cmap='GnBu')
                plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
                plt.contour(X, Y, fit, cmap=plt.cm.copper)

        rabi_signal[ind_pulse_num] = np.average(select_signal)
        preselected_signal1 = []
        preselected_signal2 = []
        preselected_signal3 = []
        preselected_signal4 = []
        # for idy in range(len(herald_signal)):
        #     if ((sReal[idy] - centers_fit[0, 0]) / sigmas_fit[0, 0] / width_threshold) ** 2 + (
        #             (sImag[idy] - centers_fit[0, 1]) / sigmas_fit[0, 1] / width_threshold) ** 2 < 1:
        #         preselected_signal1 = np.append(preselected_signal1, select_signal[idy])
        #     elif ((sReal[idy] - centers_fit[1, 0]) / sigmas_fit[1, 0] / width_threshold) ** 2 + (
        #             (sImag[idy] - centers_fit[1, 1]) / sigmas_fit[1, 1] / width_threshold) ** 2 < 1:
        #         preselected_signal2 = np.append(preselected_signal2, select_signal[idy])
        #     elif ((sReal[idy] - centers_fit[2, 0]) / sigmas_fit[2, 0] / width_threshold) ** 2 + (
        #             (sImag[idy] - centers_fit[2, 1]) / sigmas_fit[2, 1] / width_threshold) ** 2 < 1:
        #         preselected_signal3 = np.append(preselected_signal3, select_signal[idy])
        #     elif ((sReal[idy] - centers_fit[3, 0]) / sigmas_fit[3, 0] / width_threshold) ** 2 + (
        #             (sImag[idy] - centers_fit[3, 1]) / sigmas_fit[3, 1] / width_threshold) ** 2 < 1:
        #         preselected_signal4 = np.append(preselected_signal4, select_signal[idy])

        for ind_blob in range(num_blob):
            ind = ((sReal - centers_fit[ind_blob, 0]) / sigmas_fit[ind_blob, 0] / width_threshold) ** 2 + (
                    (sImag - centers_fit[ind_blob, 1]) / sigmas_fit[ind_blob, 1] / width_threshold) ** 2 < 1
            rabi_signal_preselected[ind_pulse_num, ind_blob] = np.average(select_signal[ind])
            rb_signal[ind_blob, ind_pulse_num, ind_pulse_randomize] = np.average(select_signal[ind])

plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
for ind_blob in range(num_blob):
    plt.plot(np.real(rb_signal[ind_blob, :, -1]), np.imag(rb_signal[ind_blob, :, -1]),
             label='preselect blob ' + str(ind_blob))
plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

avg_signal = np.average(rb_signal, axis=2)


n = 1
d = 2 ** n
p_array = np.zeros(num_blob)
p_err_array = np.zeros_like(p_array)
fig, ax = plt.subplots()
ax.grid(linestyle='--')
for ind_blob in range(num_blob):
    V_complex = qdf.AutoRotate(avg_signal[ind_blob, :])
    plt.plot(pulse_num, np.real(V_complex))
    toFit = np.real(V_complex)
    guess = ([0.9, np.max(toFit) - np.min(toFit), np.min(toFit)])
    opt, cov = curve_fit(randomized_benchmarking_0, ydata=toFit, xdata=pulse_num, p0=guess)
    err = (np.sqrt(np.diag(cov)))
    p_array[ind_blob] = opt[0]
    p_err_array[ind_blob] = err[0]
    plt.plot(pulse_num, randomized_benchmarking_0(pulse_num, *opt), label='Zeroth order fit')

print('p and p_err for different blobs')
print(p_array)
print(p_err_array)
parameter = np.average(p_array)
parameter_err = parameter * np.sqrt(np.sum((p_err_array / p_array) ** 2))
error = abs((d - 1) * (1 - parameter) / d)
error_err = (d - 1) * parameter_err / d
error = error / 1.875
error_err = error_err / 1.875
print('p=%.4G\u00B1%.4G, 0-order model fidelity %.4G\u00B1%.4G' % (parameter, parameter_err, 1 - error, error_err))

fig, ax = plt.subplots()
x_blob = []
y_blob = []
for ind_pulse_randomize in range(len(pulse_randomize)):
    most_pts_ind = np.argmax(blobs_matrix[ind_pulse_randomize, :, 0])
    x_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 1]]
    y_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 2]]
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
plt.contour(X, Y, fit, cmap=plt.cm.copper)
plt.plot(x_blob, y_blob, c='r')
# print(x_blob)
# print(y_blob)
plt.show()
