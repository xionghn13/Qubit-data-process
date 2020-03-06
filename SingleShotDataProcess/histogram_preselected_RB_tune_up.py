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


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\03\Data_0302\RB_heralded_AWG_qubitB_sweep_11.hdf5')
num_blob = 4

width_threshold = 2  # sigma
gg_estimate = [100, 150]

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')[:, :]
pulse_randomize = np.unique(f.getData('Multi-Qubit Pulse Generator - Randomize'))
sweep_quantity = np.unique(f.getData('IQ 1 - Frequency'))
# if len(sweep_quantity) == 1:
#         sweep_quantity = np.unique(f.getData('Multi-Qubit Pulse Generator - Amplitude, 2QB #12'))
if len(sweep_quantity) == 1:
    sweep_quantity = np.unique(f.getData('Multi-Qubit Pulse Generator - Frequency #1'))
if len(sweep_quantity) == 1:
    sweep_quantity = np.unique(f.getData('Multi-Qubit Pulse Generator - Amplitude #1'))


blobs_matrix = np.zeros((len(pulse_randomize), 4, 5))

rb_signal = np.zeros((num_blob, len(pulse_randomize), len(sweep_quantity)), dtype=complex)

rabi_signal_preselected = np.zeros(num_blob, dtype=complex)

for ind_sweep_quantity in range(len(sweep_quantity)):
    print('analyzing pulse', ind_sweep_quantity)
    for ind_pulse_randomize in range(len(pulse_randomize)):
        print('analyzing ind_pulse_randomize', ind_pulse_randomize)
        herald_signal = signal[ind_sweep_quantity * len(pulse_randomize) + ind_pulse_randomize, 0::2] * 1e6
        select_signal = signal[ind_sweep_quantity * len(pulse_randomize) + ind_pulse_randomize, 1::2] * 1e6
        sReal = np.real(herald_signal)
        sImag = np.imag(herald_signal)
        H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
        X, Y = np.meshgrid(xedges[1:], yedges[1:])
        H = H.T

        centers, sigmas = ssf.get_blob_centers(sReal, sImag, num_blob)
        heights = ssf.get_center_heights(X, Y, H, centers)
        param_mat = np.concatenate((heights.reshape(num_blob, 1), centers, sigmas), axis=1)

        params, params_err = fg.fit_gaussian(X, Y, H, param_mat)
        fit = fg.multi_gaussian(X, Y, params)

        if ind_pulse_randomize + ind_sweep_quantity == 0:
            print('-----')
            print(params)
            centers_fit = params[:, 1:3]
            sigmas_fit = params[:, 3:]
            heights_fit = params[:, 0]
            most_pts_ind = np.argmax(heights_fit)
            gg_ind = np.argmin(np.sum((centers_fit - gg_estimate) ** 2, axis=1))
            print(gg_ind)
            fig, ax = plt.subplots()
            plt.pcolormesh(X, Y, H, cmap='GnBu')
            plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
            plt.contour(X, Y, fit, cmap=plt.cm.copper)
        elif ind_sweep_quantity == len(sweep_quantity) - 1:
            blobs_matrix[ind_pulse_randomize, :, :] = params

        for ind_blob in range(num_blob):
            ind = ((sReal - centers_fit[ind_blob, 0]) / sigmas_fit[ind_blob, 0] / width_threshold) ** 2 + (
                    (sImag - centers_fit[ind_blob, 1]) / sigmas_fit[ind_blob, 1] / width_threshold) ** 2 < 1
            rabi_signal_preselected[ind_blob] = np.average(select_signal[ind])
            rb_signal[ind_blob, ind_pulse_randomize, ind_sweep_quantity] = rabi_signal_preselected[ind_blob]

avg_signal = np.average(rb_signal, axis=1)


for ind_blob in range(num_blob):
    plt.plot(np.real(avg_signal[ind_blob, :]), np.imag(avg_signal[ind_blob, :]),
             label='preselect blob' + str(ind_blob))

plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

# print(avg_signal[gg_ind, :])
# print((centers_fit[gg_ind, 0] + 1j * centers_fit[gg_ind, 1]))
distance = np.abs(avg_signal[gg_ind, :] - (centers_fit[gg_ind, 0] + 1j * centers_fit[gg_ind, 1]))


fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(sweep_quantity, distance)

fig, ax = plt.subplots()
x_blob = []
y_blob = []
for ind_pulse_randomize in range(len(pulse_randomize)):
    most_pts_ind = np.argmin(blobs_matrix[ind_pulse_randomize, :, 0])
    x_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 1]]
    y_blob += [blobs_matrix[ind_pulse_randomize, most_pts_ind, 2]]
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.contour(X, Y, fit, cmap=plt.cm.copper)
plt.plot(x_blob, y_blob, c='r')

plt.show()
