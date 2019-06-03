import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize


def gamma(n):
    return kappa / 2 * np.real(np.sqrt((1 + 1j * chi / kappa) ** 2 + 4j * chi * n / kappa) - 1)


def gamma_n_kappa_chi(n, kappa_0, chi_0):
    return kappa_0 / 2 * np.real(np.sqrt((1 + 1j * chi_0 / kappa_0) ** 2 + 4j * chi_0 * n / kappa_0) - 1)


def cavityThermalPhotonTemperature(chi, T1, T2, kappa, f, n_eff=np.linspace(0, 0.2, 101)):
    kB = 1.38e-23
    h = 6.63e-34
    gamma_phi = 1 / T2 - 1 / 2 / T1

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
    plt.ylabel('Gamma_phi(MHz)', fontsize='x-large')
    plt.title('crossing at n_eff = %.3G, temperature T = %.3GmK' % (sol, T * 1e3))
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

    kappa_array = np.linspace(0.1, 1000, 1001) * kappa

    fig, ax = plt.subplots()
    plt.plot(kappa_array / 2 / np.pi, gamma_n_kappa_chi(n_sol, kappa_array, chi))
    plt.plot(np.array([kappa_array[0], kappa_array[-1]]) / 2 / np.pi, np.array([1, 1]) * gamma_phi, '--')
    plt.xlabel('kappa/2pi (MHz)', fontsize='x-large')
    plt.ylabel('Gamma_phi(MHz)', fontsize='x-large')
    plt.title('crossing at n_eff = %.3G, temperature T = %.3GmK' % (sol, T * 1e3))
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()

    chi_array = np.linspace(0.1, 10, 1001) * chi

    fig, ax = plt.subplots()
    plt.plot(chi_array / 2 / np.pi, gamma_n_kappa_chi(n_sol, kappa, chi_array))
    plt.plot(np.array([chi_array[0], chi_array[-1]]) / 2 / np.pi, np.array([1, 1]) * gamma_phi, '--')
    plt.xlabel('chi/2pi (MHz)', fontsize='x-large')
    plt.ylabel('Gamma_phi(MHz)', fontsize='x-large')
    plt.title('crossing at n_eff = %.3G, temperature T = %.3GmK' % (sol, T * 1e3))
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    chi = 0.5 * 2 * np.pi
    T2 = 21
    T1 = 1000
    kappa = 4 * 2 * np.pi
    f = 8.894e9
    cavityThermalPhotonTemperature(chi, T1, T2, kappa, f)