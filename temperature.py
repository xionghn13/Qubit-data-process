import numpy as np

h = 6.63e-34
f = 7.96e9
kB = 1.38e-23
T0 = 46e-3

P0 = 0.54
P1 = 0.54 * 0.68
T = h * f / kB / np.log(P0 / P1)
print('T = %.3GmK' % (T * 1000))
print(np.exp( - h * f / kB / T0))
print('n_eff = %.3G' % (1 / (np.exp(h * f / kB / T0) - 1)))
print(1 / np.tanh(h * f / kB / T0))

print('---------------A B to P0 P1----------------------')

A = 0.243
B = 0.525
f = 1.1558e9
r01 = (1 - (B - A)) / (1 - (B + A))
p0 = r01 / (r01 + 1)
T = h * f / kB / np.log(r01)
print(r01)
print(p0)
print('T = %.3GmK' % (T * 1000))

print('---------------A1 B1 A2 B2 to P0 P1 P2----------------------')

A1, B1, f1 = [0.215, 0.498, 1.1558e9]
A2, B2, f2 = [0.354, 0.653, 3.8816e9]
end_point = 1.
r01 = (end_point - (B1 - A1)) / (end_point - (B1 + A1))
r02 = abs((end_point - (B2 - A2)) / (end_point - (B2 + A2)))
T1 = h * f1 / kB / np.log(r01)
T2 = h * f2 / kB / np.log(r02)
p0 = 1 / (1 / r01 + 1 / r02 + 1)
R2P = p0 / np.mean((1 - (B1 - A1), 1 - (B2 - A2)))
print(r01)
print('[P0, P1, P2] = [', p0, ',', p0 / r01, ',', p0 / r02, ']')
print('T01 = %.3GmK' % (T1 * 1000))
print('T02 = %.3GmK' % (T2 * 1000))
print('end point is', end_point)
print('ratio between P0 and end_point - r is', R2P)


print('---------------Gamma to temp----------------------')
Gamma = 2 * np.pi * 4.69
Gamma_up = 0.0067
Time_up_std = 0.0
Gamma_up_high = 1 / (1 / Gamma_up - Time_up_std)
Gamma_up_low = 1 / (1 / Gamma_up + Time_up_std)
Gamma_down = Gamma - Gamma_up
Gamma_down = 1 / 60
Time_down_std = 0.0
Gamma_down_high = 1 / (1 / Gamma_down - Time_down_std)
Gamma_down_low = 1 / (1 / Gamma_down + Time_down_std)
f = 1.1518e9
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
T = 63e-3
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
lamb = 10e-3 * 4
c = 3e8
epsilon = (1 + 1) / 2
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
print('--------------------n to temp------------------')
g = 511
fc = 7.45e3
fq = 109
kappa = 5.47 * 2 * np.pi
Delta = fc - fq
Gamma_purcell = kappa * (g / Delta) ** 2
T1_purcell = 1 / Gamma_purcell
print('T1_purcell=%.3Gus' % T1_purcell)
print('--------------------T1 ratio vs Temp------------------')
T1_lowT = 120
dT1_lowT = 18
temp_low = 13.7e-3
T1_highT = 84.3
dT1_highT = 9.5
temp_high = 15.6e-3
freq = 528.7e6
T1_low_high_ratio = T1_lowT / T1_highT
dT1_low_high_ratio = np.sqrt((dT1_lowT / T1_highT) ** 2 + (dT1_highT * T1_lowT / T1_highT ** 2) ** 2)
T1_ratio_from_temp = 1 / (1 + np.exp(- h * freq / kB / temp_low)) * (1 + np.exp(- h * freq / kB / temp_high))
print(u"T1_low_high_ratio=%.3G\u00B1%.3G, T1_ratio_from_temp=%.3G" % (T1_low_high_ratio, dT1_low_high_ratio, T1_ratio_from_temp))
