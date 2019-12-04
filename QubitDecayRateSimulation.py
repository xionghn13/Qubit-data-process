from QubitDataProcessPackages import *
import Single_small_junction as ssj


# for dielectric loss estimation
N = 50
EL = 0.619
EC = 1.184
EJ = 1.967
# loss_tan = 1/e-6
# Q_cap = 1 / loss_tan
Q_cap = 3.9e5
T = 20e-3
flux = 0.5
[pem01, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
T1diel01 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap, freq, pem01, T)
[pem23, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 2, 3))
T1diel23 = 1e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap, freq, pem23, T)

[qpem02, freq] = np.abs(ssj.qp_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 2))
T1qp02 = 1e6 / ssj.relaxation_rate_qp_finiteTemp(EL, EC, EJ, freq, qpem02, T)


# [pem01, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
# T1diel01 = 1e6 / ssj.relaxation_rate_qp_finiteTemp(EL, EC, EJ, Q_cap, freq, pem01, T)

print('10 time ', T1diel01, 'us, ', '32 time ', T1diel23, 'us. T1qp02 is', T1qp02, 'us.')