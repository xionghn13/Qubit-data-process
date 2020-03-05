import numpy as np


def T1_curve(t, A, T1, B):
    r = A * np.exp(- t / T1) + B
    # print('T1 = %.3Gus, A = %.3G, B = %.3G' % (T1 / 1000, A, B))
    return r


def DoubleExp_curve(t, A, TR, B, Tqp, lamb):
    r = A * np.exp(- t / TR + lamb * (np.exp(- t / Tqp) - 1)) + B
    # print('T1 = %.3Gus, A = %.3G, B = %.3G' % (T1 / 1000, A, B))
    return r


def TwoExp_curve(t, A, T1, B, T2, C):
    r = A * np.exp(- t / T1) + B * np.exp(- t / T2) + C
    return r


def rabi_curve(t, A, T1, B, Tpi, phi0):
    r = A * np.exp(- t / T1) * np.cos(t / Tpi * np.pi + phi0) + B
    return r


def rabi_two_exp_curve(t, A, T1, B, Tpi, phi0, T_out, C):
    r = A * np.exp(- t / T1) * np.cos(t / Tpi * np.pi + phi0) + B * np.exp(-t / T_out) + C
    return r


def transient_time_curve(power_dBm_Gamma_r, Gamma_in, Gamma_out, A):
    power_dBm = power_dBm_Gamma_r[0]
    Gamma_r = power_dBm_Gamma_r[1]
    power = 1e-3 * 10 ** (power_dBm / 10)
    Gamma = Gamma_in + power * A * Gamma_out / (Gamma_r ** 2 + 2 * power * A)
    time = 1 / Gamma
    return time


def lorenztian(f, f0, kappa, A, B):
    t = A / ((f - f0) ** 2 + (kappa / 2) ** 2) + B
    return t

# move FitTransientTime to FittingFunc


# move AutoRotate to DataManipulationFunc