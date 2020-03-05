import numpy as np
import Labber
from matplotlib import pyplot as plt
import SingleShotDataProcess.SingleShotFunc as ssf
import SingleShotDataProcess.FitGaussians as fg
from scipy.optimize import curve_fit
import FunctionLib as qdf

def osc_func(x, amp, freq, offset1, offset2):
    return amp * np.cos(2 * np.pi * freq * (x - offset1)) - offset2


# constants
kB = 1.38e-23
h = 6.626e-34
############################################################
# Vary heralding wait time

file_path = 'C:\SC Lab\Labber\data\Augustus 18\\2020\\02\Data_0229\\'
file_name = 'TuneUp_heralded_AWG_qubitB_18.hdf5'

f = Labber.LogFile(file_path + file_name)

num_blob = 4
width_threshold = 2
measurement_type = file_name.split('_')[0]

# sweep_quantity_name = 'Multi-Qubit Pulse Generator - Amplitude #1'
sweep_quantity_name = 'Multi-Qubit Pulse Generator - Frequency #1'

signal = f.getData('AlazarTech Signal Demodulator - Channel A - Demodulated values')
pulse_num = np.unique(f.getData('Multi-Qubit Pulse Generator - # of pulses'))
sweep_quantity = np.unique( f.getData(sweep_quantity_name))

print(pulse_num)

# print(sweep_quantity)
rabi_signal = np.zeros((len(pulse_num), len(sweep_quantity)), dtype=complex)
rabi_signal_preselected = np.zeros((len(pulse_num), len(sweep_quantity), num_blob), dtype=complex)

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

for ind_sweep_quantity in range(len(sweep_quantity)):
    print('ind_sweep_quantity', ind_sweep_quantity)
    for ind_pulse_num in range(len(pulse_num)):
        print('ind_pulse_num', ind_pulse_num)
        herald_signal = signal[ind_sweep_quantity * len(pulse_num) + ind_pulse_num, 0::2] * 1e6
        select_signal = signal[ind_sweep_quantity * len(pulse_num) + ind_pulse_num, 1::2] * 1e6
        sReal = np.real(herald_signal)
        sImag = np.imag(herald_signal)
        rabi_signal[ind_pulse_num, ind_sweep_quantity] = np.average(select_signal)
        for ind_blob in range(num_blob):
            ind = ((sReal - centers_fit[ind_blob, 0]) / sigmas_fit[ind_blob, 0] / width_threshold) ** 2 + (
                    (sImag - centers_fit[ind_blob, 1]) / sigmas_fit[ind_blob, 1] / width_threshold) ** 2 < 1
            rabi_signal_preselected[ind_pulse_num, ind_sweep_quantity, ind_blob] = np.average(select_signal[ind])



fig, ax = plt.subplots()
plt.pcolormesh(X, Y, H, cmap='GnBu')
plt.scatter(centers_fit[:, 0], centers_fit[:, 1], c='r')
plt.contour(X, Y, fit, cmap=plt.cm.copper)


plt.plot(np.real(rabi_signal[0, :]), np.imag(rabi_signal[0, :]), label='no selection data')
for ind_blob in range(num_blob):
    plt.plot(np.real(rabi_signal_preselected[0, :, ind_blob]), np.imag(rabi_signal_preselected[0, :, ind_blob]),
             label='preselect blob ' + str(ind_blob))
plt.legend()
plt.xlabel('I (uV)')
plt.ylabel('Q (uV)')

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for ind_pulse_num in range(len(pulse_num)):
    V_complex = rabi_signal_preselected[ind_pulse_num, :, most_pts_ind]
    V_real = V_complex.real
    plt.plot(sweep_quantity, V_real, label='pulse number %.5G' % pulse_num[ind_pulse_num])
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
plt.ylabel('real')
plt.tight_layout()

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for ind_pulse_num in range(len(pulse_num)):
    V_complex = rabi_signal_preselected[ind_pulse_num, :, most_pts_ind]
    V_real = V_complex.imag
    plt.plot(sweep_quantity, V_real, label='pulse number %.5G' % pulse_num[ind_pulse_num])
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
plt.ylabel('imag')
plt.tight_layout()

fig, ax = plt.subplots()
ax.grid(linestyle='--')
for ind_sweep_quantity in range(len(sweep_quantity)):
    V_complex = rabi_signal_preselected[:, ind_sweep_quantity, most_pts_ind]
    V_real = V_complex.real
    plt.plot(pulse_num, V_real, label='%.5G' % sweep_quantity[ind_sweep_quantity])
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
plt.tight_layout()
#
plt.show()
