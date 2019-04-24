import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit


def reflection(fPower, f0, gamma_f, P0, A):
    f = fPower[0]
    num_f = len(f)
    Power = fPower[1]
    num_Power = len(Power)
    f = np.tensordot(f, np.ones(num_Power), axes=0)
    Power = np.tensordot(np.ones(num_f), Power, axes=0)
    r = 1 - 2 * P0 * (1 + 2j * (f - f0) / gamma_f) / (1 + 4 * (f - f0) ** 2 / gamma_f ** 2 + A * Power)
    return r


NumPower = 5
BackPowerdBm = 5
OneToneFreq = np.linspace(8.08, 8.1, 500)
OneToneFreq = np.concatenate(([0], OneToneFreq, [0]))
OneTonePowerdBm = np.linspace(-25, -5, NumPower)
OneTonePower = 1e-3 * 10 ** (OneTonePowerdBm / 10)
BackPower = 1e-3 * 10 ** (BackPowerdBm / 10)
f0 = 8.09
gamma_f = 1e-3
P0 = 0.7
A = 5e4
amp_cor = 1
# df = 0.15
df = 0.00
guess = [f0, gamma_f, P0, A, amp_cor]
bounds = ([1, 0, 0, 0, 0.5], [20, 20, 1, np.inf, 1.5])

OneToneR = reflection([OneToneFreq, OneTonePower], f0, gamma_f, P0, A)
BackR = reflection([OneToneFreq, [BackPower]], f0 - df, gamma_f, P0, A)
OneToneR2 = OneToneR / BackR
opt, cov, LargerFreqRange, FittedComplex = qsf.fitReflectionCircles(OneToneFreq, OneTonePowerdBm, OneToneR2,
                                                                    guess, bounds)
f0_fit, gamma_f_fit, P0_fit, A_fit, amp_cor = opt

fig, ax = plt.subplots()
for i in range(NumPower):
    plt.plot(OneToneFreq[1: -1], np.abs(OneToneR[1: -1, i]))
plt.plot(OneToneFreq[1: -1], np.abs(BackR[1: -1]), '--')
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Abs', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
for i in range(NumPower):
    plt.plot(np.real(OneToneR2[100: -100, i]), np.imag(OneToneR2[100: -100, i]))
    plt.plot(np.real(OneToneR[:, i]), np.imag(OneToneR[:, i]), 'r')
plt.plot([-2, 2], [0, 0], '--')
plt.plot([1], [0], 'ro')
plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.xlabel('Re', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_aspect('equal')
plt.show()