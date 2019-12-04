# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
from qutip import *
from matplotlib import pyplot as plt

rho1 = rand_dm(2)
rho2 = rand_dm(2)

##############single qubit tomography##############
# ideal measurement
rho1_measure = 0.5 * (qeye(2) + expect(sigmax(), rho1) * sigmax() + expect(sigmay(), rho1) * sigmay() + expect(sigmaz(),
                                                                                                               rho1) * sigmaz())
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

##############two qubits tomography##############
rho = tensor(rho1, rho2)
# ideal measurement
rho_measure = 0.25 * (tensor(qeye(2), qeye(2))
                      + expect(tensor(qeye(2), sigmax()), rho) * tensor(qeye(2), sigmax())
                      + expect(tensor(qeye(2), sigmay()), rho) * tensor(qeye(2), sigmay())
                      + expect(tensor(qeye(2), sigmaz()), rho) * tensor(qeye(2), sigmaz())
                      + expect(tensor(sigmax(), qeye(2)), rho) * tensor(sigmax(), qeye(2))
                      + expect(tensor(sigmax(), sigmax()), rho) * tensor(sigmax(), sigmax())
                      + expect(tensor(sigmax(), sigmay()), rho) * tensor(sigmax(), sigmay())
                      + expect(tensor(sigmax(), sigmaz()), rho) * tensor(sigmax(), sigmaz())
                      + expect(tensor(sigmay(), qeye(2)), rho) * tensor(sigmay(), qeye(2))
                      + expect(tensor(sigmay(), sigmax()), rho) * tensor(sigmay(), sigmax())
                      + expect(tensor(sigmay(), sigmay()), rho) * tensor(sigmay(), sigmay())
                      + expect(tensor(sigmay(), sigmaz()), rho) * tensor(sigmay(), sigmaz())
                      + expect(tensor(sigmaz(), qeye(2)), rho) * tensor(sigmaz(), qeye(2))
                      + expect(tensor(sigmaz(), sigmax()), rho) * tensor(sigmaz(), sigmax())
                      + expect(tensor(sigmaz(), sigmay()), rho) * tensor(sigmaz(), sigmay())
                      + expect(tensor(sigmaz(), sigmaz()), rho) * tensor(sigmaz(), sigmaz())
                      )
# real measurement with voltages
# beta coeffcients are experimentally determined from either single shot readout or from Rabi calibration + temperature measurement
betaII = (np.random.randint(1, 1000) + 1j * np.random.randint(1, 1000)) * 1e-6
betaIZ = (np.random.randint(1, 1000) + 1j * np.random.randint(1, 1000)) * 1e-6
betaZI = (np.random.randint(1, 1000) + 1j * np.random.randint(1, 1000)) * 1e-6
betaZZ = (np.random.randint(1, 1000) + 1j * np.random.randint(1, 1000)) * 1e-6
MII = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmaz(), qeye(
    2)) + betaZZ * tensor(sigmaz(), sigmaz()))
MIX = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmay()) + betaZI * tensor(sigmaz(), qeye(
    2)) + betaZZ * tensor(sigmaz(), sigmay()))
MIY = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmax()) + betaZI * tensor(sigmaz(), qeye(
    2)) + betaZZ * tensor(sigmaz(), sigmax()))
MIZ = 0.25 * (betaII * tensor(qeye(2), qeye(2)) - betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmaz(), qeye(
    2)) - betaZZ * tensor(sigmaz(), sigmaz()))

MXI = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmay(), qeye(
    2)) + betaZZ * tensor(sigmay(), sigmaz()))
MXX = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmay()) + betaZI * tensor(sigmay(), qeye(
    2)) + betaZZ * tensor(sigmay(), sigmay()))
MXY = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmax()) + betaZI * tensor(sigmay(), qeye(
    2)) + betaZZ * tensor(sigmay(), sigmax()))
MXZ = 0.25 * (betaII * tensor(qeye(2), qeye(2)) - betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmay(), qeye(
    2)) - betaZZ * tensor(sigmay(), sigmaz()))

MYI = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmax(), qeye(
    2)) + betaZZ * tensor(sigmax(), sigmaz()))
MYX = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmay()) + betaZI * tensor(sigmax(), qeye(
    2)) + betaZZ * tensor(sigmax(), sigmay()))
MYY = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmax()) + betaZI * tensor(sigmax(), qeye(
    2)) + betaZZ * tensor(sigmax(), sigmax()))
MYZ = 0.25 * (betaII * tensor(qeye(2), qeye(2)) - betaIZ * tensor(qeye(2), sigmaz()) + betaZI * tensor(sigmax(), qeye(
    2)) - betaZZ * tensor(sigmax(), sigmaz()))

MZI = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmaz()) - betaZI * tensor(sigmaz(), qeye(
    2)) - betaZZ * tensor(sigmaz(), sigmaz()))
MZX = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmay()) - betaZI * tensor(sigmaz(), qeye(
    2)) - betaZZ * tensor(sigmaz(), sigmay()))
MZY = 0.25 * (betaII * tensor(qeye(2), qeye(2)) + betaIZ * tensor(qeye(2), sigmax()) - betaZI * tensor(sigmaz(), qeye(
    2)) - betaZZ * tensor(sigmaz(), sigmax()))
MZZ = 0.25 * (betaII * tensor(qeye(2), qeye(2)) - betaIZ * tensor(qeye(2), sigmaz()) - betaZI * tensor(sigmaz(), qeye(
    2)) + betaZZ * tensor(sigmaz(), sigmaz()))
# In experiments, mij is determined via dispersive measurements.
mII = expect(MII, rho)  # Measurement in ground state
mIX = expect(MIX, rho)  # measurement after doing nothing on qubit A and X pi/2 pulse on qubit B
mIY = expect(MIY, rho)  # measurement after doing nothing on qubit A and Y -pi/2 pulse on qubit B
mIZ = expect(MIZ, rho)  # measurement after flipping the state of qubit B
mXI = expect(MXI, rho)  # measurement after doing X pi/2 pulse on qubit A and nothing on qubit B
mXX = expect(MXX, rho)  # measurement after doing X pi/2 pulse on qubit A and X pi/2 pulse on qubit B
mXY = expect(MXY, rho)  # measurement after doing X pi/2 pulse on qubit A and Y -pi/2 pulse on qubit B
mXZ = expect(MXZ, rho)  # measurement after doing X pi/2 pulse on qubit A and flipping the state of qubit B
mYI = expect(MYI, rho)  # measurement after doing Y -pi/2 pulse on qubit A and nothing on qubit B
mYX = expect(MYX, rho)  # measurement after doing Y -pi/2 pulse on qubit A and X pi/2 pulse on qubit B
mYY = expect(MYY, rho)  # measurement after doing Y -pi/2 pulse on qubit A and Y -pi/2 pulse on qubit B
mYZ = expect(MYZ, rho)  # measurement after doing Y -pi/2 pulse on qubit A and flipping the state of qubit B
mZI = expect(MZI, rho)  # measurement after flipping the state of qubit A and nothing on qubit B
mZX = expect(MZX, rho)  # measurement after flipping the state of qubit A and X pi/2 pulse on qubit B
mZY = expect(MZY, rho)  # measurement after flipping the state of qubit A and Y -pi/2 pulse on qubit B
mZZ = expect(MZZ, rho)  # measurement after flipping the state of qubit A and flipping the state of qubit B
# write the matrix in order II, IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ
# mij = matrix dot avg
measurement_matrix = 0.25 * np.array([[betaII, 0, 0, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, betaZI, 0, 0, betaZZ],  # mII
                                      [betaII, 0, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, betaZI, 0, betaZZ, 0],  # mIX
                                      [betaII, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, betaZI, betaZZ, 0, 0],  # mIY
                                      [betaII, 0, 0, -betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, betaZI, 0, 0, -betaZZ],  # mIZ
                                      [betaII, 0, 0, betaIZ, 0, 0, 0, 0, betaZI, 0, 0, betaZZ, 0, 0, 0, 0],  # mXI
                                      [betaII, 0, betaIZ, 0, 0, 0, 0, 0, betaZI, 0, betaZZ, 0, 0, 0, 0, 0],  # mXX
                                      [betaII, betaIZ, 0, 0, 0, 0, 0, 0, betaZI, betaZZ, 0, 0, 0, 0, 0, 0],  # mXY
                                      [betaII, 0, 0, -betaIZ, 0, 0, 0, 0, betaZI, 0, 0, -betaZZ, 0, 0, 0, 0],  # mXZ
                                      [betaII, 0, 0, betaIZ, betaZI, 0, 0, betaZZ, 0, 0, 0, 0, 0, 0, 0, 0],  # mYI
                                      [betaII, 0, betaIZ, 0, betaZI, 0, betaZZ, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # mYX
                                      [betaII, betaIZ, 0, 0, betaZI, betaZZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # mYY
                                      [betaII, 0, 0, -betaIZ, betaZI, 0, 0, -betaZZ, 0, 0, 0, 0, 0, 0, 0, 0],  # mYZ
                                      [betaII, 0, 0, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, -betaZI, 0, 0, -betaZZ],  # mZI
                                      [betaII, 0, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, -betaZI, 0, -betaZZ, 0],  # mZX
                                      [betaII, betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -betaZI, -betaZZ, 0, 0],  # mZY
                                      [betaII, 0, 0, -betaIZ, 0, 0, 0, 0, 0, 0, 0, 0, -betaZI, 0, 0, betaZZ]])  # mZZ
print(measurement_matrix)

avgII, avgIX, avgIY, avgIZ, \
avgXI, avgXX, avgXY, avgXZ, \
avgYI, avgYX, avgYY, avgYZ, \
avgZI, avgZX, avgZY, avgZZ \
    = np.linalg.inv(measurement_matrix).dot(np.array([
    mII, mIX, mIY, mIZ,
    mXI, mXX, mXY, mXZ,
    mYI, mYX, mYY, mYZ,
    mZI, mZX, mZY, mZZ]).transpose()).transpose()
rho_reconstructed = 0.25 * (
            avgII * tensor(qeye(2), qeye(2)) + avgIX * tensor(qeye(2), sigmax()) + avgIY * tensor(qeye(2),
                                                                                                  sigmay()) + avgIZ * tensor(
        qeye(2), sigmaz()) +
            avgXI * tensor(sigmax(), qeye(2)) + avgXX * tensor(sigmax(), sigmax()) + avgXY * tensor(sigmax(),
                                                                                                    sigmay()) + avgXZ * tensor(
        sigmax(), sigmaz()) +
            avgYI * tensor(sigmay(), qeye(2)) + avgYX * tensor(sigmay(), sigmax()) + avgYY * tensor(sigmay(),
                                                                                                    sigmay()) + avgYZ * tensor(
        sigmay(), sigmaz()) +
            avgZI * tensor(sigmaz(), qeye(2)) + avgZX * tensor(sigmaz(), sigmax()) + avgZY * tensor(sigmaz(),
                                                                                                    sigmay()) + avgZZ * tensor(
        sigmaz(), sigmaz()))
# plotting
matrix_histogram_complex(rho)
matrix_histogram_complex(rho_measure)
matrix_histogram_complex(rho_reconstructed)

plt.show()
