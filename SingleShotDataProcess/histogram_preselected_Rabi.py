import numpy as np
import Labber
from matplotlib import pyplot as plt
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
from scipy.optimize import curve_fit


def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time
file_path = 'C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0217\\'
file_name = 'Ramsey_heralded_AWG_qubitB_3.hdf5'
f = Labber.LogFile(file_path + file_name)
num_blob = 4

measurement_type = file_name.split('_')[0]

sweep_quantity_dict = {
    'Rabi': 'Multi-Qubit Pulse Generator - Amplitude #1',
    'Ramsey': 'Multi-Qubit Pulse Generator - Pulse spacing',
}
signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')
sweep_quantity = f.getData(sweep_quantity_dict[measurement_type])[0]
rabi_signal = np.zeros(len(sweep_quantity), dtype=complex)
rabi_signal_preselected_1 = np.zeros(len(sweep_quantity), dtype=complex)
rabi_signal_preselected_2 = np.zeros(len(sweep_quantity), dtype=complex)
rabi_signal_preselected_3 = np.zeros(len(sweep_quantity), dtype=complex)
rabi_signal_preselected_4 = np.zeros(len(sweep_quantity), dtype=complex)


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

# print(centers_fit)
fit = fg.multi_gaussian(X, Y, params)

fig, ax = plt.subplots()
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
plt.contour(X, Y, fit, cmap=plt.cm.copper)

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for idx in range(len(sweep_quantity)):
    herald_signal = signal[idx, 0::2] * 1e6
    select_signal = signal[idx, 1::2] * 1e6
    sReal = np.real(herald_signal)
    sImag = np.imag(herald_signal)
    if idx == 0:
        H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(X, Y, H, cmap='GnBu')
        # plt.colorbar()
    rabi_signal[idx] = np.average(select_signal)
    preselected_signal1 = []
    preselected_signal2 = []
    preselected_signal3 = []
    preselected_signal4 = []
    for idy in range(len(herald_signal)):
        if ((sReal[idy] - centers_fit[0, 0]) / sigmas_fit[0, 0] / 2) ** 2 + (
                (sImag[idy] - centers_fit[0, 1]) / sigmas_fit[0, 1] / 2) ** 2 < 1:
            preselected_signal1 = np.append(preselected_signal1, select_signal[idy])
        elif ((sReal[idy] - centers_fit[1, 0]) / sigmas_fit[1, 0] / 2) ** 2 + (
                (sImag[idy] - centers_fit[1, 1]) / sigmas_fit[1, 1] / 2) ** 2 < 1:
            preselected_signal2 = np.append(preselected_signal2, select_signal[idy])
        elif ((sReal[idy] - centers_fit[2, 0]) / sigmas_fit[2, 0] / 2) ** 2 + (
                (sImag[idy] - centers_fit[2, 1]) / sigmas_fit[2, 1] / 2) ** 2 < 1:
            preselected_signal3 = np.append(preselected_signal3, select_signal[idy])
        elif ((sReal[idy] - centers_fit[3, 0]) / sigmas_fit[3, 0] / 2) ** 2 + (
                (sImag[idy] - centers_fit[3, 1]) / sigmas_fit[3, 1] / 2) ** 2 < 1:
            preselected_signal4 = np.append(preselected_signal4, select_signal[idy])

    rabi_signal_preselected_1[idx] = np.average(preselected_signal1)
    rabi_signal_preselected_2[idx] = np.average(preselected_signal2)
    rabi_signal_preselected_3[idx] = np.average(preselected_signal3)
    rabi_signal_preselected_4[idx] = np.average(preselected_signal4)

plt.plot(np.real(rabi_signal), np.imag(rabi_signal), label='raw Rabi')
plt.plot(np.real(rabi_signal_preselected_1), np.imag(rabi_signal_preselected_1), label='preselect gg')
plt.plot(np.real(rabi_signal_preselected_2), np.imag(rabi_signal_preselected_2), label='preselect eg')
plt.plot(np.real(rabi_signal_preselected_3), np.imag(rabi_signal_preselected_3), label='preselect ge')
plt.plot(np.real(rabi_signal_preselected_4), np.imag(rabi_signal_preselected_4), label='preselect ee')
plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

freq_guess = 4 / sweep_quantity[-1]
fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(sweep_quantity, np.real(rabi_signal_preselected_1))
guess = ([np.max(np.real(rabi_signal_preselected_1)) - np.min(np.real(rabi_signal_preselected_1)), freq_guess, 0,
          (np.max(np.real(rabi_signal_preselected_1)) + np.min(np.real(rabi_signal_preselected_1))) / 2])
print(guess)
opt, cov = curve_fit(osc_func, ydata=np.real(rabi_signal_preselected_1), xdata=sweep_quantity, p0=guess)
axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
plt.plot(axis_nice, osc_func(axis_nice, *opt))
plt.title('Pi pulse: %.3G' % (0.5 / opt[1]))

fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(sweep_quantity, np.real(rabi_signal_preselected_2))
guess = ([np.max(np.real(rabi_signal_preselected_2)) - np.min(np.real(rabi_signal_preselected_2)), freq_guess, 0,
          np.real(rabi_signal_preselected_2)[0]])
opt, cov = curve_fit(osc_func, ydata=np.real(rabi_signal_preselected_2), xdata=sweep_quantity, p0=guess)
axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
plt.plot(axis_nice, osc_func(axis_nice, *opt))
plt.title('Pi pulse: %.3G' % (0.5 / opt[1]))
#
plt.show()
