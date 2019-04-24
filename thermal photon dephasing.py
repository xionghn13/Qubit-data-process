import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

n_eff = np.linspace(0, 0.2, 101)

chi = 1.3 * 2 * np.pi
T2 = 12.1
T1 = 77.8
gamma_phi = 1 / T2 - 1 / 2 / T1
kappa = 6.8 * 2 * np.pi
f = 7.79e9
kB = 1.38e-23
h = 6.63e-34


def gamma(n):
    return kappa / 2 * np.real(np.sqrt((1 + 1j * chi / kappa) ** 2 + 4j * chi * n / kappa) - 1)


def eqn(n):
    return gamma(n) - gamma_phi


sol = optimize.fsolve(eqn, x0=0)
n_sol = sol[0]
gamma_n = gamma(n_eff)
T = h * f / kB / np.log(1 + 1 / sol)

print(sol)
fig, ax = plt.subplots()
plt.plot(n_eff, gamma_n)
plt.plot(np.array([n_eff[0], n_eff[-1]]), np.array([1, 1]) * gamma_phi, '--')
plt.xlabel('n_eff', fontsize='x-large')
plt.ylabel('Gamma(MHz)', fontsize='x-large')
plt.title('crossing at n_eff = %.3G, temperature T = %.3GmK' % (sol, T * 1e3))
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()
