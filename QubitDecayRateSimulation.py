from QubitDataProcessPackages import *
import Single_small_junction as ssj


# for dielectric loss estimation
N = 50
EL = 0.4187
EC = 2.5225
EJ = 6.5852

EC = 1.0645
EJ = 5.2024
EL = 0.7235

# EC = 1.0690
# EJ = 3.8684
# EL = 0.7701

# loss_tan = 1/e-6
# Q_cap = 1 / loss_tan
Q_cap = 3.0e5 / 8.54 * 6.23
T = 1e-9
epsilon = 0

flux = 0.50
[pem01, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
T1diel01 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap / (freq / 1.1558) ** epsilon, freq, pem01, T)
print(freq)
[pem23, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 2, 3))
print(freq)
T1diel23 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap / (freq / 1.1558) ** epsilon, freq, pem23, T)
[pem12, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 1, 2))
print(freq)
T1diel12 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap / (freq / 1.1558) ** epsilon, freq, pem12, T)
[pem02, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 2))
print(freq)
T1diel02 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap / (freq / 1.1558) ** epsilon, freq, pem02, T)

[qpem02, freq] = np.abs(ssj.qp_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 2))
T1qp02 = 1e6 / ssj.relaxation_rate_qp_finiteTemp(EL, EC, EJ, freq, qpem02, T)


# [pem01, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
# T1diel01 = 1e6 / ssj.relaxation_rate_qp_finiteTemp(EL, EC, EJ, Q_cap, freq, pem01, T)

print('10 time', T1diel01, 'us, ', '32 time', T1diel23, 'us. 21 time', T1diel12, 'us. T102 is', T1diel02, 'us.')