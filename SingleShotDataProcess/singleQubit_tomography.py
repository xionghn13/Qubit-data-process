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
f = Labber.LogFile('C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0217\SingleQubit_tomo_qubitB_2.hdf5')

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

params = fg.fitgaussian(X, Y, H, param_mat)
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

for pulse_idx in range(len(signal[:, 0])):
    preselected_signal = []
    herald_signal = signal[pulse_idx, 0::2] * 1e6
    select_signal = signal[pulse_idx, 1::2] * 1e6
    sReal = np.real(herald_signal)
    sImag = np.imag(herald_signal)
    for record_idx in range(len(herald_signal)):
        if ((sReal[record_idx] - centers_fit[gg_ind, 0]) / sigmas_fit[gg_ind, 0] / 2) ** 2 + (
                (sImag[record_idx] - centers_fit[gg_ind, 1]) / sigmas_fit[gg_ind, 1] / 2) ** 2 < 1:
            preselected_signal = np.append(preselected_signal, select_signal[record_idx])
    preselected_data[pulse_idx] = np.average(preselected_signal)

s0 = preselected_data[0]
sz = preselected_data[6]
print(s0, sz)
betaI = s0 + sz
betaZ = s0 - sz
# Start in ground state
m = preselected_data[0:3]

measurement_matrix = 0.5 * np.array([[0, 0, betaZ], [0, betaZ, 0], [betaZ, 0, 0]])
avgX, avgY, avgZ = np.linalg.inv(measurement_matrix).dot(m.transpose() - 0.5 * betaI).transpose()
rho_reconstructed = 0.5 * (qeye(2) + avgX * sigmax() + avgY * sigmay() + avgZ * sigmaz())
matrix_histogram_complex(rho_reconstructed)
# rho_ideal = ket2dm(basis(2,0))
# plt.title(fidelity(rho_ideal, rho_reconstructed))
# matrix_histogram_complex(rho_ideal)
# Start in excited state
print('ground state')
print(m)
m = preselected_data[6:9]
avgX, avgY, avgZ = np.linalg.inv(measurement_matrix).dot(m.transpose() - 0.5 * betaI).transpose()
rho_reconstructed = 0.5 * (qeye(2) + avgX * sigmax() + avgY * sigmay() + avgZ * sigmaz())
matrix_histogram_complex(rho_reconstructed)
# rho_ideal = ket2dm(basis(2,1))
# plt.title(fidelity(rho_ideal, rho_reconstructed))
# matrix_histogram_complex(rho_ideal)
# Start in superposition state
print('pi pulse')
print(m)
m = preselected_data[9:]
avgX, avgY, avgZ = np.linalg.inv(measurement_matrix).dot(m.transpose() - 0.5 * betaI).transpose()
rho_reconstructed = 0.5 * (qeye(2) + avgX * sigmax() + avgY * sigmay() + avgZ * sigmaz())
matrix_histogram_complex(rho_reconstructed)
# rho_ideal = ket2dm(rx(phi=np.pi/2)*basis(2,0))
# plt.title(fidelity(rho_ideal, rho_reconstructed))
# matrix_histogram_complex(rho_ideal)
print('pi/2 pulse')
print(m)
#####################################################################################################
# sI = np.array([[1,0],[0,1]])
# sX = np.array([[0,1],[1,0]])
# sY = np.array([[0,-1j],[1j,0]])
# sZ = np.array([[1,0],[0,-1]])
#
# M = np.zeros((3, 2, 2), dtype = complex)
# M[0,:,:] = betaI +betaZ*sZ
# M[1,:,:] = betaI +betaZ*sY
# M[2,:,:] = betaI +betaZ*sX

# def density_matrix(t1,t2,t3):
#     tau = np.array([[t1, t2+1j*t3], [0, 1-t1]])
#     rho = np.conj(tau.transpose())*tau / np.trace(np.conj(tau.transpose())*tau)
#     return rho
#
# def likelihood(x):
#     dist = 0
#     for idx in range(3):
#         dist = dist + np.real((m[idx] - np.trace(M[idx, :, :].dot(density_matrix(x[0],x[1],x[2])))))**2 + np.imag((m[idx] - np.trace(M[idx, :, :].dot(density_matrix(x[0],x[1],x[2])))))**2
#     return dist
#
# guess = np.ones(3)*0.5
# # guess[0] = 1
# res = minimize(likelihood, guess, method='nelder-mead')
# t = res.x
# rho_reconstructed_mle = density_matrix(*t)
# matrix_histogram_complex(rho_reconstructed_mle)
plt.show()
