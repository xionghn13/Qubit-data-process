from FunctionLib import *
from scipy.optimize import curve_fit



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



def fit_lorenztian(frequency, V_sq_abs):
    MaxAbs_0 = np.max(V_sq_abs)
    V_sq_abs_fit = V_sq_abs / MaxAbs_0
    # V_sq_abs_fit = V_sq_abs
    MaxAbs = np.max(V_sq_abs_fit)
    MinAbs = np.min(V_sq_abs_fit)
    MaxInd = V_sq_abs_fit.argmax()
    # print(MaxInd)
    f0_guess = frequency[MaxInd]
    kappa_guess = (frequency[-1] - frequency[0]) / 4
    B_guess = MinAbs
    A_guess = (MaxAbs - MinAbs) * (kappa_guess / 2) ** 2

    guess = ([f0_guess, kappa_guess, A_guess, B_guess])
    # print(guess)
    bounds = (
        (frequency[0], 0, 0, 0),
        (frequency[-1], kappa_guess * 4, MaxAbs * 10, MaxAbs)
    )

    opt, cov = curve_fit(lorenztian, frequency, V_sq_abs_fit, guess, bounds=bounds)
    # f0_fit, kappa_fit, A_fit, B_fit = qopt
    # print(qopt)
    err = np.sqrt(np.diag(cov))
    opt[2:] *= MaxAbs_0
    err[2:] *= MaxAbs_0
    # print(qopt)
    fit_abs = lorenztian(frequency, *opt)

    return opt, err, fit_abs