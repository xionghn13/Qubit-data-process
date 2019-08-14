import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from QubitDecayFunc import T1_curve, rabi_curve, AutoRotate, DoubleExp_curve
import ExtractDataFunc as edf
import os

DataPath = 'C:\\Users/admin\Labber\Data/2019/07/Data_0718/'
T1File = 't1_2019-07-18-11-00-56.hdf5'

[Time, Complex] = edf.readT1Labber(DataPath + T1File)
ReadoutPower = edf.readQubitPowerLabber(DataPath + T1File)
RComplex = Complex * 10 ** (- ReadoutPower / 20)
RComplex = - AutoRotate(RComplex)
y_data = np.array(RComplex.real, dtype='float64')
x_data = np.array(Time, dtype='float64')

B_guess = y_data[-1]
A_guess = y_data[0].real - B_guess
T1_guess = x_data[-1] / 10

# Double exp fit
try:
    opt, cov = curve_fit(DoubleExp_curve, x_data, y_data,
                         p0=[A_guess, T1_guess, B_guess, T1_guess * 0.1, 0.5],
                         maxfev=300000)
    # print('guess = %s' % str([A_guess, T1_guess, B_guess, T1_guess * 0.1, 1]))
    # print('Double exp fit opt = %s' % str(opt))
except RuntimeError:
    print("Error - curve_fit failed")
    opt = np.array([A_guess, T1_guess, B_guess, T1_guess, 1])
    cov = np.zeros([len(opt), len(opt)])

A_double_fit, TR_fit, B_double_fit, Tqp_fit, lamb_fit = opt
A_double_std, TR_std, B_double_std, Tqp_std, lamb_std = np.sqrt(cov.diagonal())

TimeFit = np.linspace(Time.min(), Time.max(), 200)

FitR_double = DoubleExp_curve(TimeFit, A_double_fit, TR_fit, B_double_fit, Tqp_fit, lamb_fit)
# y_pred = DoubleExp_curve(Time, A_fit, TR_fit, B_fit, Tqp_fit, lamb_fit)
# ParamList = ['A', 'TR/ns', 'B', 'Tqp/ns', 'lambda']
# Single exp fit
opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
A_single_fit, T1_fit, B_single_fit = opt
A_single_std, T1_std, B_single_std = np.sqrt(cov.diagonal())

# TimeFit = np.linspace(Time.min(), Time.max(), 200)

FitR_single = T1_curve(TimeFit, A_single_fit, T1_fit, B_single_fit)
# y_pred = T1_curve(Time, A_fit, T1_fit, B_fit)
# y_guess = T1_curve(Time, A_guess, T1_guess, B_guess)
# ParamList = ['A', 'Decay time/ns', 'B']
print('double offset: %.3G, single offset: %.3G' % (B_double_fit, B_single_fit))
B_fit = (B_double_fit + B_single_fit) / 2
fig, ax = plt.subplots()
plt.plot(Time / 1000, y_data - B_double_fit, 'o')
plt.plot(TimeFit / 1000, FitR_double - B_double_fit)
plt.plot(TimeFit / 1000, FitR_single - B_single_fit)

plt.title('Double Exp Fit: TR=%.3G$\pm$%.2Gus, Tqp=%.3G$\pm$%.2Gus, nqp=%.3G$\pm$%.2G\n'
          'Single Exp Fit: T1=%.3G$\pm$%.2Gus' % (
          TR_fit / 1000, TR_std / 1000, Tqp_fit / 1000, Tqp_std / 1000, lamb_fit, lamb_std, T1_fit / 1000,
          T1_std / 1000))

plt.xlabel('Time/us', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.ylim([1e-6, 1e-3])
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_yscale('log')
plt.show()
