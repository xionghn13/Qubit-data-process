import numpy as np
from qutip import *
from matplotlib import pyplot as plt
import scipy.interpolate as itp



def vec2dm(vec):
    return 0.5 * (qeye(2) + vec[0] * sigmax() + vec[1] * sigmay() + vec[2] * sigmaz())


##############single qubit tomography##############
n_rho1 = np.array([0.8, 0.6, 0])
n_noise = np.random.randn(3) / 20
n_rho1_noise = n_rho1 + n_noise
# ideal measurement
rho1 = vec2dm(n_rho1)
rho1_measure = vec2dm(n_rho1_noise) + qeye(2) * np.random.randn() / 20
# print(n_rho1_noise)
print('Fidelity before construction', fidelity(rho1_measure, rho1))
# real measurement with voltages
psi_k_arr = np.concatenate((sigmaz().eigenstates()[1], sigmax().eigenstates()[1], sigmay().eigenstates()[1]))
# psi_k_arr = sigmaz().eigenstates()[1]
n_iter = 10
epsilon = 0.01
Delta = 0.001

dimension = rho1_measure.shape[0]
num_states = len(psi_k_arr)
Pi_k_arr = np.zeros_like(psi_k_arr)
f_k_arr = np.zeros(num_states)
p_k_arr = np.zeros_like(f_k_arr)
for i in range(len(psi_k_arr)):
    f_k_arr[i] = expect(rho1_measure, psi_k_arr[i])
    Pi_k_arr[i] = ket2dm(psi_k_arr[i])


rho_iter = qeye(dimension) / dimension
for i in range(n_iter):
    # print(rho_iter)
    for j in range(len(psi_k_arr)):
        p_k_arr[j] = expect(rho_iter, psi_k_arr[j])
        # print(psi_k_arr[j])
    # print(p_k_arr)
    R = np.sum(f_k_arr / p_k_arr * Pi_k_arr)

    if np.sqrt((((R - (R * rho_iter).tr()) * rho_iter) ** 2).tr()) < Delta:
        break
    else:
        epsilon_arr = np.array([[-epsilon], [0], [epsilon]])
        epsilon_mat = np.concatenate((epsilon_arr ** 2, epsilon_arr, np.ones_like(epsilon_arr)), axis=1)
        logL_arr = np.zeros_like(epsilon_arr)

        def new_rho(old_rho, epsilon_p):
            p_old_k_arr = np.zeros_like(f_k_arr)
            for j in range(len(psi_k_arr)):
                p_old_k_arr[j] = expect(old_rho, psi_k_arr[j])
            R_old = np.sum(f_k_arr / p_old_k_arr * Pi_k_arr)
            # print('test', (R_old - 1))
            # print(epsilon_p)
            temp = 1 + epsilon_p * (R_old - (R_old * old_rho).tr())
            # print('R_old', R_old)
            new_dm = temp * old_rho * temp
            new_dm /= new_dm.tr()
            return new_dm

        def logL(rho_test):
            p_old_k_arr = np.zeros_like(f_k_arr)
            # print(rho_test)
            for j in range(len(psi_k_arr)):
                # print('expect', expect(rho_test, psi_k_arr[j]))
                p_old_k_arr[j] = expect(rho_test, psi_k_arr[j])
            # print('p', np.log(p_old_k_arr))
            return np.sum(f_k_arr * np.log(p_old_k_arr))
        # print('rho_iter', rho_iter)

        for j in range(3):
            logL_arr[j] = logL(new_rho(rho_iter, epsilon_arr[j, 0]))

        a, b, c = np.matmul(np.linalg.inv(epsilon_mat), logL_arr)
        epsilon_opt = -b / 2 / a
        rho_iter = new_rho(rho_iter, epsilon_opt[0])
        fid = fidelity(rho_iter, rho1)
        print('Fidelity is', fid)
print('R\n', np.array(R))
print('ideal rho\n', np.array(rho1))
print('measured rho\n', np.array(rho1_measure))
print('reconstructed rho\n', np.array(rho_iter))
print('Fidelity is', fidelity(rho_iter, rho1))

##############two qubit tomography##############
# rho = ket2dm(bell_state())
# n_iter = 10
# epsilon = 0.01
# Delta = 0.001
#
# psi_k_arr = np.concatenate([sigmaz().eigenstates()[1], sigmax().eigenstates()[1], sigmay().eigenstates()[1]])
#
# rho_measure = 0.25 * (tensor(qeye(2), qeye(2))
#                       + expect(tensor(qeye(2), sigmax()), rho) * tensor(qeye(2), sigmax())
#                       + expect(tensor(qeye(2), sigmay()), rho) * tensor(qeye(2), sigmay())
#                       + expect(tensor(qeye(2), sigmaz()), rho) * tensor(qeye(2), sigmaz())
#                       + expect(tensor(sigmax(), qeye(2)), rho) * tensor(sigmax(), qeye(2))
#                       + expect(tensor(sigmax(), sigmax()), rho) * tensor(sigmax(), sigmax())
#                       + expect(tensor(sigmax(), sigmay()), rho) * tensor(sigmax(), sigmay())
#                       + expect(tensor(sigmax(), sigmaz()), rho) * tensor(sigmax(), sigmaz())
#                       + expect(tensor(sigmay(), qeye(2)), rho) * tensor(sigmay(), qeye(2))
#                       + expect(tensor(sigmay(), sigmax()), rho) * tensor(sigmay(), sigmax())
#                       + expect(tensor(sigmay(), sigmay()), rho) * tensor(sigmay(), sigmay())
#                       + expect(tensor(sigmay(), sigmaz()), rho) * tensor(sigmay(), sigmaz())
#                       + expect(tensor(sigmaz(), qeye(2)), rho) * tensor(sigmaz(), qeye(2))
#                       + expect(tensor(sigmaz(), sigmax()), rho) * tensor(sigmaz(), sigmax())
#                       + expect(tensor(sigmaz(), sigmay()), rho) * tensor(sigmaz(), sigmay())
#                       + expect(tensor(sigmaz(), sigmaz()), rho) * tensor(sigmaz(), sigmaz())
#                       )
#
#
#
# dimension = rho1_measure.shape[0]
# num_states = len(psi_k_arr)
# Pi_k_arr = np.zeros_like(psi_k_arr)
# f_k_arr = np.zeros(num_states)
# p_k_arr = np.zeros_like(f_k_arr)
# for i in range(len(psi_k_arr)):
#     f_k_arr[i] = expect(rho1_measure, psi_k_arr[i])
#     Pi_k_arr[i] = ket2dm(psi_k_arr[i])
#
#
# rho_iter = qeye(dimension) / dimension
# for i in range(n_iter):
#     # print(rho_iter)
#     for j in range(len(psi_k_arr)):
#         p_k_arr[j] = expect(rho_iter, psi_k_arr[j])
#         # print(psi_k_arr[j])
#     # print(p_k_arr)
#     R = np.sum(f_k_arr / p_k_arr * Pi_k_arr)
#
#     if np.sqrt((((R - (R * rho_iter).tr()) * rho_iter) ** 2).tr()) < Delta:
#         break
#     else:
#         epsilon_arr = np.array([[-epsilon], [0], [epsilon]])
#         epsilon_mat = np.concatenate((epsilon_arr ** 2, epsilon_arr, np.ones_like(epsilon_arr)), axis=1)
#         logL_arr = np.zeros_like(epsilon_arr)
#
#         def new_rho(old_rho, epsilon_p):
#             p_old_k_arr = np.zeros_like(f_k_arr)
#             for j in range(len(psi_k_arr)):
#                 p_old_k_arr[j] = expect(old_rho, psi_k_arr[j])
#             R_old = np.sum(f_k_arr / p_old_k_arr * Pi_k_arr)
#             # print('test', (R_old - 1))
#             # print(epsilon_p)
#             temp = 1 + epsilon_p * (R_old - (R_old * old_rho).tr())
#             # print('R_old', R_old)
#             new_dm = temp * old_rho * temp
#             new_dm /= new_dm.tr()
#             return new_dm
#
#         def logL(rho_test):
#             p_old_k_arr = np.zeros_like(f_k_arr)
#             # print(rho_test)
#             for j in range(len(psi_k_arr)):
#                 # print('expect', expect(rho_test, psi_k_arr[j]))
#                 p_old_k_arr[j] = expect(rho_test, psi_k_arr[j])
#             # print('p', np.log(p_old_k_arr))
#             return np.sum(f_k_arr * np.log(p_old_k_arr))
#         # print('rho_iter', rho_iter)
#
#         for j in range(3):
#             logL_arr[j] = logL(new_rho(rho_iter, epsilon_arr[j, 0]))
#
#         a, b, c = np.matmul(np.linalg.inv(epsilon_mat), logL_arr)
#         epsilon_opt = -b / 2 / a
#         rho_iter = new_rho(rho_iter, epsilon_opt[0])
#         fid = fidelity(rho_iter, rho1)
        # print('Fidelity is', fid)
# print('R\n', np.array(R))
# print('ideal rho\n', np.array(rho1))
# print('measured rho\n', np.array(rho1_measure))
# print('reconstructed rho\n', np.array(rho_iter))
# print('Fidelity is', fidelity(rho_iter, rho1))