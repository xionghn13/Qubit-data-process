import numpy as np

h = 6.63e-34
f = 0.5e9
kB = 1.38e-23
T0 = 0.07

P0 = 0.2
T = h * f / kB / np.log((1 - P0) / P0)
print('T = %.3GmK' % (T * 1000))
print(np.exp(- h * f / kB / T0))
print(1 / np.tanh(h * f / kB / T0))

print('---------------A B to P0 P1----------------------')

A = 0.131
B = 0.396
f = 1.2045e9
r01 = (1 - (B - A / 2)) / (1 - (B + A / 2))
p0 = r01 / (r01 + 1)
T = h * f / kB / np.log(r01)
print(p0)
print('T = %.3GmK' % (T * 1000))

print('---------------Gamma to temp----------------------')
Gamma = 2 * np.pi * 4.69
Gamma_up = 1 / 11.6
Time_up_std = 0.0
Gamma_up_high = 1 / (1 / Gamma_up - Time_up_std)
Gamma_up_low = 1 / (1 / Gamma_up + Time_up_std)
Gamma_down = Gamma - Gamma_up
Gamma_down = 2.95
Time_down_std = 0.0
Gamma_down_high = 1 / (1 / Gamma_down - Time_down_std)
Gamma_down_low = 1 / (1 / Gamma_down + Time_down_std)
f = 4.104e9
T = h * f / kB / np.log(Gamma_down / Gamma_up)
T_high = h * f / kB / np.log(Gamma_down_high / Gamma_up_low)
T_low = h * f / kB / np.log(Gamma_down_low / Gamma_up_high)
print('Gamma_up = %.3GMHz, Gamma_down = %.3GMHz' % (Gamma_up, Gamma_down))
print('T = %.3GmK' % (T * 1000))
print('T_high = %.3GmK' % (T_high * 1000))
print('T_low = %.3GmK' % (T_low * 1000))
print('---------------temp to Gamma----------------------')
Gamma_down = 2 * np.pi * 5
f = 11.5e9
T = 60e-3
Gamma_up = Gamma_down * np.exp(- h * f / kB / T)
print('Gamma_up = %.3GMHz, Gamma_down = %.3GMHz' % (Gamma_up, Gamma_down))
print('T_2 = %.3Gus' % (1 / Gamma_up))
print('---------------temp to freq----------------------')
Gamma_down = 2 * np.pi * 5
Gamma_up = 1 / 200
T = 70e-3
freq = T * kB / h * np.log(Gamma_down / Gamma_up)
print('f = %.3GGHz' % (freq / 1e9))
print('--------------------freq to wavelength----------------')
f = 0.02e9
c = 3e8
epsilon = 2
T = 1 / f
lamb = c / np.sqrt(epsilon) * T * 1e3
print('lambda = %.3Gmm' % lamb)
print('--------------------wavelength to freq----------------')
lamb = 4e-3 * 2
c = 3e8
epsilon = (1 + 11.9) / 2
T = lamb * np.sqrt(epsilon) / c
f = 1 / T * 1e-9
print('freq = %.3GGHz' % f)
print('---------------temp to S_ratio----------------------')
f1 = 0.0001e9
f2 = 1e9
T = 70e-3
S_ratio = (1 / np.tanh(h * f1 / kB / T / 2) + 1) / (1 / np.tanh(h * f2 / kB / T / 2) + 1)
print('S_ratio = %.3G' % (S_ratio))
print('--------------------n to temp------------------')
n = 0.016
f = 7.79e9
T = h * f / kB / np.log(1 + 1 / n)
print('T = %.3GmK' % (T * 1e3))
