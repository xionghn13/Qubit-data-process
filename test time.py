from qutip import *
import numpy as np
import Single_small_junction as ssj
import pprofile

prof = pprofile.Profile()
with prof():
    # Code to profile
    vec = np.linspace(1, 51, 50)
    matrix = np.tensordot(vec, vec, axes=0)

    E_l = 1.2384393267942007
    E_c = 2.4511689738320404
    # EC = 3.3
    E_j = 4.783739177406206
    I0 = 0.15267272001393203
    I_period = 10.132145138431971
    N = 12900
    level_num = 3
    Imin = 0
    Imax = 15
    I_array = np.linspace(Imin, Imax, 1)

    phi_ext = (I_array - I0) / I_period * 2 * np.pi
    if level_num < 3:
        level_num = 3
    phi_num = phi_ext.size

    n0x = np.zeros((phi_num, level_num - 1))
    n1x = np.zeros((phi_num, level_num - 2))
    n2x = np.zeros((phi_num, level_num - 3))
    phi0x = np.zeros((phi_num, level_num - 1))
    phi1x = np.zeros((phi_num, level_num - 2))
    phi2x = np.zeros((phi_num, level_num - 3))
    freq0x = np.zeros((phi_num, level_num - 1))
    freq1x = np.zeros((phi_num, level_num - 2))
    freq2x = np.zeros((phi_num, level_num - 3))

    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)

    for idx, phi0 in enumerate(phi_ext):
        matrix = ssj.bare_hamiltonian(N, E_l, E_c, E_j, phi0).full()
        w, v = np.linalg.eig(matrix)
        Hamiltonian = ssj.bare_hamiltonian(N, E_l, E_c, E_j, phi0)
        energies, states = Hamiltonian.eigenstates()
        # for idy in range(level_num - 1):
        #     freq0x[idx, idy] = energies[idy + 1] - energies[0]
        #     phi0x[idx, idy] = abs(phi.matrix_element(states[0], states[idy + 1]))
        #     n0x[idx, idy] = abs(na.matrix_element(states[0], states[idy + 1]))
        # for idy in range(level_num - 2):
        #     freq1x[idx, idy] = energies[idy + 2] - energies[1]
        #     phi1x[idx, idy] = abs(phi.matrix_element(states[1], states[idy + 2]))
        #     n1x[idx, idy] = abs(na.matrix_element(states[1], states[idy + 2]))
        # for idy in range(level_num - 3):
        #     freq2x[idx, idy] = energies[idy + 3] - energies[2]
        #     phi2x[idx, idy] = abs(phi.matrix_element(states[2], states[idy + 3]))
        #     n2x[idx, idy] = abs(na.matrix_element(states[2], states[idy + 3]))
prof.print_stats()
