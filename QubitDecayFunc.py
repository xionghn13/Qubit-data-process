import numpy as np
from scipy.optimize import curve_fit


def T1_curve(t, A, T1, B):
    r = A * np.exp(- t / T1) + B
    # print('T1 = %.3Gus, A = %.3G, B = %.3G' % (T1 / 1000, A, B))
    return r

def DoubleExp_curve(t, A, TR, B, Tqp, lamb):
    r = A * np.exp(- t / TR + lamb * (np.exp( - t / Tqp) - 1)) + B
    # print('T1 = %.3Gus, A = %.3G, B = %.3G' % (T1 / 1000, A, B))
    return r

def rabi_curve(t, A, T1, B, Tpi, phi0):
    r = A * np.exp(- t / T1) * np.cos(t / Tpi * np.pi + phi0) + B
    return r


def transient_time_curve(power_dBm_Gamma_r, Gamma_in, Gamma_out, A):
    power_dBm = power_dBm_Gamma_r[0]
    Gamma_r = power_dBm_Gamma_r[1]
    power = 1e-3 * 10 ** (power_dBm / 10)
    Gamma = Gamma_in + power * A * Gamma_out / (Gamma_r ** 2 + 2 * power * A)
    time = 1 / Gamma
    return time


def FitTransientTime(DrivePowerArray, Gamma_r, OptMatrix, power_for_plot=[]):
    if len(power_for_plot) == 0:
        power_for_plot = DrivePowerArray
    x_data = [DrivePowerArray, Gamma_r]
    x_for_plot = [power_for_plot, Gamma_r]
    y_data = OptMatrix[1, :] / 1000
    pow_mean = 1e-3 * 10 ** (x_data[0].mean() / 10)
    pow_ratio_guess = Gamma_r ** 2 / 2 / pow_mean
    Gamma_in_guess = 1 / y_data.max()
    Gamma_out_guess = 1 / y_data.min()
    guess = ([Gamma_in_guess, Gamma_out_guess, pow_ratio_guess])
    bounds = (
        (0, 0, 0),
        (10 * Gamma_r, 10 * Gamma_r, np.inf)
    )
    opt, cov = curve_fit(transient_time_curve, x_data, y_data, p0=guess, bounds=bounds)
    Gamma_in_fit, Gamma_out_fit, pow_ratio_fit = opt
    FitTime = transient_time_curve(x_for_plot, Gamma_in_fit, Gamma_out_fit, pow_ratio_fit)
    return opt, cov, FitTime


def AutoRotate(complex):
    complex_out = complex
    disp_complex = complex_out - complex_out[0]
    for i in range(3):
        max_disp_ind = np.argmax(np.abs(disp_complex))
        disp_complex = complex_out - complex_out[max_disp_ind]
    max_disp_ind = np.argmax(np.abs(disp_complex))
    max_disp = disp_complex[max_disp_ind]
    normalized_max_disp = max_disp / np.abs(max_disp)
    # print(normalized_max_disp)
    if normalized_max_disp.real < 0:
        normalized_max_disp = -normalized_max_disp
    CM = np.mean(complex_out) * 0 # origin
    complex_out -= CM
    complex_out /= normalized_max_disp
    complex_out += CM
    return complex_out