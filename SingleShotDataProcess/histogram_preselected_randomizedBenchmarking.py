import numpy as np
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
import Labber
from matplotlib import pyplot as plt
from qutip import *
from scipy.optimize import curve_fit


def gaussian4(x, a1, x1, a2, x2, a3, x3, a4, x4, sigma):
    return a1 * np.exp(-(x - x1) ** 2 / sigma ** 2) + a2 * np.exp(-(x - x2) ** 2 / sigma ** 2) + \
           a3 * np.exp(-(x - x3) ** 2 / sigma ** 2) + a4 * np.exp(-(x - x4) ** 2 / sigma ** 2)


def gaussian2d(coord, a, x0, y0, sigma):
    return a * np.exp((-(coord[0] - x0) ** 2 - (coord[1] - y0) ** 2) / sigma ** 2)


def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0217\RB_heralded_AWG_qubitB_2.hdf5')
num_blob = 4

width_threshold = 2  # sigma

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')
pulseNum = f.getData('Multi-Qubit Pulse Generator - Number of Cliffords')[0]
pulseAmp = f.getData('Multi-Qubit Pulse Generator - Randomize')[:, 0]
rb_signal = np.zeros((len(pulseNum), len(pulseAmp)), dtype=complex)
rabi_signal = np.zeros(len(pulseNum), dtype=complex)
rabi_signal_preselected_1 = np.zeros(len(pulseNum), dtype=complex)
rabi_signal_preselected_2 = np.zeros(len(pulseNum), dtype=complex)
rabi_signal_preselected_3 = np.zeros(len(pulseNum), dtype=complex)
rabi_signal_preselected_4 = np.zeros(len(pulseNum), dtype=complex)

for pulseAmp_idx in range(len(pulseAmp)):
    for idx in range(len(pulseNum)):
        herald_signal = signal[pulseAmp_idx * len(pulseNum) + idx, 0::2] * 1e6
        select_signal = signal[pulseAmp_idx * len(pulseNum) + idx, 1::2] * 1e6
        sReal = np.real(herald_signal)
        sImag = np.imag(herald_signal)
        if idx == 0 and pulseAmp_idx == 0:
            H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
            X, Y = np.meshgrid(xedges[1:], yedges[1:])
            H = H.T

            centers, sigmas = ssf.getBlobCenters(sReal, sImag, num_blob)
            heights = ssf.getCenterHeights(X, Y, H, centers)
            param_mat = np.concatenate((heights.reshape(num_blob, 1), centers, sigmas), axis=1)

            params = fg.fitgaussian(X, Y, H, param_mat)
            print('-----')
            print(params)
            centers_fit = params[:, 1:3]
            sigmas_fit = params[:, 3:]

            # print(centers_fit)
            fit = fg.multi_gaussian(X, Y, params)

            fig, ax = plt.subplots()
            plt.pcolormesh(X, Y, H, cmap='GnBu')
            plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
            plt.contour(X, Y, fit, cmap=plt.cm.copper)

        rabi_signal[idx] = np.average(select_signal)
        preselected_signal1 = []
        preselected_signal2 = []
        preselected_signal3 = []
        preselected_signal4 = []
        for idy in range(len(herald_signal)):
            if ((sReal[idy] - centers_fit[0, 0]) / sigmas_fit[0, 0] / width_threshold) ** 2 + (
                    (sImag[idy] - centers_fit[0, 1]) / sigmas_fit[0, 1] / width_threshold) ** 2 < 1:
                preselected_signal1 = np.append(preselected_signal1, select_signal[idy])
            elif ((sReal[idy] - centers_fit[1, 0]) / sigmas_fit[1, 0] / width_threshold) ** 2 + (
                    (sImag[idy] - centers_fit[1, 1]) / sigmas_fit[1, 1] / width_threshold) ** 2 < 1:
                preselected_signal2 = np.append(preselected_signal2, select_signal[idy])
            elif ((sReal[idy] - centers_fit[2, 0]) / sigmas_fit[2, 0] / width_threshold) ** 2 + (
                    (sImag[idy] - centers_fit[2, 1]) / sigmas_fit[2, 1] / width_threshold) ** 2 < 1:
                preselected_signal3 = np.append(preselected_signal3, select_signal[idy])
            elif ((sReal[idy] - centers_fit[3, 0]) / sigmas_fit[3, 0] / width_threshold) ** 2 + (
                    (sImag[idy] - centers_fit[3, 1]) / sigmas_fit[3, 1] / width_threshold) ** 2 < 1:
                preselected_signal4 = np.append(preselected_signal4, select_signal[idy])

        rabi_signal_preselected_1[idx] = np.average(preselected_signal1)
        rabi_signal_preselected_2[idx] = np.average(preselected_signal2)
        rabi_signal_preselected_3[idx] = np.average(preselected_signal3)
        rabi_signal_preselected_4[idx] = np.average(preselected_signal4)
        rb_signal[idx, pulseAmp_idx] = np.average(preselected_signal3)

plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
plt.plot(np.real(rabi_signal_preselected_1), np.imag(rabi_signal_preselected_1), label='preselect gg')
plt.plot(np.real(rabi_signal_preselected_2), np.imag(rabi_signal_preselected_2), label='preselect eg')
plt.plot(np.real(rabi_signal_preselected_3), np.imag(rabi_signal_preselected_3), label='preselect ge')
plt.plot(np.real(rabi_signal_preselected_4), np.imag(rabi_signal_preselected_4), label='preselect ee')
plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

avg_signal = np.average(rb_signal, axis=1)

plt.figure(2)
plt.plot(pulseNum, np.real(avg_signal))
plt.figure(3)
plt.plot(pulseNum, np.imag(avg_signal))
plt.figure(4)
plt.plot(pulseNum, abs(avg_signal))


def randomized_benchmarking_0(x, p, a, b):
    return a * p ** x + b


n = 1
d = 2 ** n
toFit = np.imag(avg_signal)
guess = ([0.9, np.max(toFit) - np.min(toFit), np.min(toFit)])
opt, cov = curve_fit(randomized_benchmarking_0, ydata=toFit, xdata=pulseNum, p0=guess)
err = (np.sqrt(np.diag(cov)))
parameter = (opt[0])
parameter_err = err[0]
error = abs((d - 1) * (1 - parameter) / d)
error_err = (d - 1) * parameter_err / d
error = error / 1.875
error_err = error_err / 1.875
print('0-order model fidelity', (1 - error))
print('0-order model fidelity error', (error_err))
plt.figure(3)
plt.plot(pulseNum, randomized_benchmarking_0(pulseNum, *opt), label='Zeroth order fit')

plt.show()
