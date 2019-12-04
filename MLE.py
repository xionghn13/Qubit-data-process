import numpy as np
from qutip import *
from matplotlib import pyplot as plt


def vec2dm(vec):
    return 0.5 * (qeye(2) + vec[0] * sigmax() + vec[1] * sigmay() + vec[2] * sigmaz())
##############single qubit tomography##############
n_rho1 = np.array([1, 0, 0])
n_noise = np.random.randn(3) / 20
n_rho1_noise = n_rho1 + n_noise
# ideal measurement
rho1_measure = vec2dm(n_rho1)
# real measurement with voltages
# beta coeffcients are experimentally determined from either single shot readout or from Rabi calibration + temperature measurement
betaI = (np.random.randint(1, 1000) + 1j * np.random.randint(1, 1000)) * 1e-6
betaZ = (np.random.randint(1, 1000) - 1j * np.random.randint(1, 1000)) * 1e-6
mI = 0.5 * expect(betaI * qeye(2) + betaZ * sigmaz(), rho1)  # Measurement in ground state
mZ = 0.5 * expect(betaI * qeye(2) - betaZ * sigmaz(), rho1)  # measurement after flipping the state
mX = 0.5 * expect(betaI * qeye(2) + betaZ * sigmay(), rho1)  # measurement after doing X pi/2 pulse
mY = 0.5 * expect(betaI * qeye(2) + betaZ * sigmax(), rho1)  # measurement after doing Y -pi/2 pulse
avgI = (mI + mZ) / betaI  # this is in principle equal to 1 and finding it is not necessary
avgX = (2 * mY - avgI * betaI) / betaZ
avgY = (2 * mX - avgI * betaI) / betaZ
avgZ = (mI - mZ) / betaZ
# smarter math
measurement_matrix = 0.5 * np.array(
    [[betaI, 0, 0, betaZ], [betaI, 0, betaZ, 0], [betaI, betaZ, 0, 0], [betaI, 0, 0, -betaZ]])
avgI, avgX, avgY, avgZ = np.linalg.inv(measurement_matrix).dot(np.array([mI, mX, mY, mZ]).transpose()).transpose()
rho1_recostructed = 0.5 * (avgI * qeye(2) + avgX * sigmax() + avgY * sigmay() + avgZ * sigmaz())
# plotting
# matrix_histogram_complex(rho1) # direct plotting
# matrix_histogram_complex(rho1_measure) #ideal measurement plotting
# matrix_histogram_complex(rho1_recostructed) #measurement with voltages and rotations plotting