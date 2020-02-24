import numpy as np
import Labber
from matplotlib import pyplot as plt
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
from scipy.optimize import curve_fit
import QubitDecayFunc as qdf

def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
file_path = 'C:\SC Lab\Labber\\11112019_back to waveguide\\2020\\02\Data_0224\\'
file_name = 'rabi_histogram_2.hdf5'
f = Labber.LogFile(file_path + file_name)
num_blob = 1
width_threshold = 2
measurement_type = file_name.split('_')[0]

sweep_quantity_dict = {
    # 'Rabi': 'Multi-Qubit Pulse Generator - Amplitude #1',
    'rabi': 'Pulse Generator - Width #1',
    # 'Ramsey': 'Multi-Qubit Pulse Generator - Pulse spacing',
}
signal = f.getData('Alazar - Channel A - Demodulated values')
# signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')

sweep_quantity = f.getData(sweep_quantity_dict[measurement_type])[0]
num_sweep_quantity = len(sweep_quantity)
rabi_signal = np.zeros(num_sweep_quantity, dtype=complex)


for ind in range(num_sweep_quantity):
    signal_array = signal[ind, :] * 1e6
    rabi_signal[ind] = np.average(signal_array)
    sReal = np.real(signal_array)
    sImag = np.imag(signal_array)
    H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
    X, Y = np.meshgrid(xedges[1:], yedges[1:])
    H = H.T

    centers, sigmas = ssf.getBlobCenters(sReal, sImag, num_blob)
    heights = ssf.getCenterHeights(X, Y, H, centers)
    param_mat = np.concatenate((heights.reshape(num_blob, 1), centers, sigmas), axis=1)

    params, params_err = fg.fitgaussian(X, Y, H, param_mat)

    centers_fit = params[:, 1:3]
    sigmas_fit = params[:, 3:]
    heights_fit = params[:, 0]

    print('-----')
    print(params)
    print(params_err)
    print('-----')

    most_pts_ind = np.argmax(heights_fit)
    # most_pts_ind = 0
    # print(centers_fit)
    fit = fg.multi_gaussian(X, Y, params)

    fig, ax = plt.subplots()
    plt.pcolormesh(X, Y, H, cmap='GnBu')
    plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
    plt.contour(X, Y, fit, cmap=plt.cm.copper)

fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
plt.legend()
plt.xlabel('real (uV)')
plt.ylabel('imag (uV)')

#
plt.show()
