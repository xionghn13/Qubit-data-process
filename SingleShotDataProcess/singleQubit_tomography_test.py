import numpy as np
from qutip import *
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
import Labber

##############single qubit tomography##############
# beta calibration

# Gate sequence in Labber is I, X2p, Y2m
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0217\SingleQubit_tomo_qubitB_test.hdf5')

num_blob = 4
gg_estimate = [0, -2000]

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')

preselected_data = np.zeros(len(signal[:, 0]), dtype=complex)

herald_signal = signal[0, 0::2] * 1e6
select_signal = signal[0, 1::2] * 1e6

sReal = np.real(herald_signal)
sImag = np.imag(herald_signal)
H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
X, Y = np.meshgrid(xedges[1:], yedges[1:])
H = H.T

centers, sigmas = ssf.getBlobCenters(sReal, sImag, num_blob)
heights = ssf.getCenterHeights(X, Y, H, centers)
param_mat = np.concatenate((heights.reshape(num_blob, 1), centers, sigmas), axis=1)

params, params_err = fg.fitgaussian(X, Y, H, param_mat)
print('-----')
print(params)
centers_fit = params[:, 1:3]
sigmas_fit = params[:, 3:]
gg_ind = np.argmin(np.sum((centers_fit - gg_estimate) ** 2, axis=1))
print(gg_ind)
fit = fg.multi_gaussian(X, Y, params)

fig, ax = plt.subplots()
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
plt.contour(X, Y, fit, cmap=plt.cm.copper)
preselected_signal = [[], []]

for pulse_idx in range(len(signal[:, 0])):
    herald_signal = signal[pulse_idx, 0::2] * 1e6
    sReal = np.real(herald_signal)
    sImag = np.imag(herald_signal)
    select_signal = signal[pulse_idx, 1::2] * 1e6
    for record_idx in range(len(herald_signal)):
        if ((sReal[record_idx] - centers_fit[gg_ind, 0]) / sigmas_fit[gg_ind, 0] / 2) ** 2 + (
                (sImag[record_idx] - centers_fit[gg_ind, 1]) / sigmas_fit[gg_ind, 1] / 2) ** 2 < 1:
            preselected_signal[pulse_idx] = np.append(preselected_signal[pulse_idx], select_signal[record_idx])
    preselected_data[pulse_idx] = np.average(preselected_signal[pulse_idx])

g_signal = preselected_signal[0]
print(len(g_signal))
sReal = np.real(g_signal)
sImag = np.imag(g_signal)
H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
X, Y = np.meshgrid(xedges[1:], yedges[1:])
H = H.T

fig, ax = plt.subplots()
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')


pi_signal = preselected_signal[1]
print(len(pi_signal))
sReal = np.real(pi_signal)
sImag = np.imag(pi_signal)
H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
X, Y = np.meshgrid(xedges[1:], yedges[1:])
H = H.T

fig, ax = plt.subplots()
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')


plt.show()
