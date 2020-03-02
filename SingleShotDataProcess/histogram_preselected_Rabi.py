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
file_path = 'C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0228\\'
file_name = 'Rabi_heralded_AWG_qubitB_4.hdf5'
f = Labber.LogFile(file_path + file_name)

num_blob = 4
width_threshold = 2
measurement_type = file_name.split('_')[0]

sweep_quantity_dict = {
    'Rabi': 'Multi-Qubit Pulse Generator - Amplitude #1',
    # 'Rabi': 'Multi-Qubit Pulse Generator - Width',
    # 'Ramsey': 'Multi-Qubit Pulse Generator - Pulse spacing',
    'Ramsey': 'Multi-Qubit Pulse Generator - Parameter #1',

}
sweep_quantity_name = sweep_quantity_dict[measurement_type]
signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')
sweep_quantity = f.getData(sweep_quantity_name)[0]
if measurement_type == 'Ramsey':
    sweep_quantity *= 1e6
if sweep_quantity_name == 'Multi-Qubit Pulse Generator - Width':
    sweep_quantity *= 1e6
# print(sweep_quantity)
rabi_signal = np.zeros(len(sweep_quantity), dtype=complex)
rabi_signal_preselected = np.zeros((len(sweep_quantity), num_blob), dtype=complex)

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
heights_fit = params[:, 0]
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
for ind_sweep_quantity in range(len(sweep_quantity)):
    herald_signal = signal[ind_sweep_quantity, 0::2] * 1e6
    select_signal = signal[ind_sweep_quantity, 1::2] * 1e6
    sReal = np.real(herald_signal)
    sImag = np.imag(herald_signal)
    if ind_sweep_quantity == 0:
        H, xedges, yedges = np.histogram2d(sReal, sImag, bins=[100, 100])
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(X, Y, H, cmap='GnBu')
        # plt.colorbar()
    rabi_signal[ind_sweep_quantity] = np.average(select_signal)
    preselected_signal1 = []
    preselected_signal2 = []
    preselected_signal3 = []
    preselected_signal4 = []
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
V_complex = qdf.AutoRotate(rabi_signal_preselected[:, most_pts_ind])
V_real = V_complex.real
fig, ax = plt.subplots()
ax.grid(linestyle='--')
plt.plot(sweep_quantity, V_real)
A_guess = np.max(V_real) - np.min(V_real)
B_guess = (np.max(V_real) + np.min(V_real)) / 2
T1_guess = sweep_quantity[-1]
Tpi_guess = sweep_quantity[-1] / 10
phi0_guess = 0
try:
    guess = [A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess]
    opt, cov = curve_fit(qdf.rabi_curve, ydata=V_real, xdata=sweep_quantity, p0=guess)
    err = np.sqrt(np.diag(cov))
    axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
    plt.plot(axis_nice, qdf.rabi_curve(axis_nice, *opt))
    plt.title('half period: %.5G\u00B1%.4G, decay constant: %.5G\u00B1%.4G' % (opt[3], err[3], opt[1], err[1]))
except RuntimeError:
    guess = [A_guess, T1_guess, B_guess]
    opt, cov = curve_fit(qdf.T1_curve, ydata=V_real, xdata=sweep_quantity, p0=guess)
    err = np.sqrt(np.diag(cov))
    axis_nice = np.linspace(sweep_quantity[0], sweep_quantity[-1], 1001)
    plt.plot(axis_nice, qdf.T1_curve(axis_nice, *opt))
    plt.title('decay constant: %.5G\u00B1%.4G' % (opt[1], err[1]))

#
plt.show()
