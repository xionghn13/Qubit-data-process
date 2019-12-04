from QubitDataProcessPackages import *


def thermal_photon(frequency, temperature):
    photon = 1 / (np.exp(h * frequency / k_B / temperature) - 1)
    return photon


def effective_temperature(frequency, photon_number):
    temperature_effective = h * frequency / k_B / np.log(1 + 1 / photon_number)
    return temperature_effective


def attenuated_photon(att_vec, n_ph_mat):
    num_f = n_ph_mat.shape[1]
    # print(np.shape(att_vec))
    num_T = np.shape(att_vec)[1] + 1
    n_ph_attenuated = np.zeros((num_T, num_f))
    # print(num_T)
    for i in range(num_T):
        if i == 0:
            n_ph_attenuated[i, :] = n_ph_mat[i, :] / np.prod(att_vec, axis=1)
        else:
            n_ph_attenuated[i, :] = n_ph_mat[i, :] / np.prod(att_vec[:, i - 1:], axis=1) * (att_vec[:, i - 1] - 1)
    return n_ph_attenuated


def pow_law(x, a, b):
    t = - a * x ** b
    return t

def component_att_from_measurement(data_file, freq_vec, lowT_factor, interp_for_att=True):
    if data_file.endswith('hdf5'):
        if data_file.endswith('2.hdf5'):
            [freq, S21] = edf.readVNA4portS21(data_file)
        else:
            [freq, S21] = edf.readVNAS21(data_file)
    else:
        data = np.genfromtxt(data_file, skip_header=3, delimiter=';')
        freq = data[:, 0] / 1e9
        S21 = data[:, 1] + 1j * data[:, 2]
    freq *= 1e9
    S21_dB = 20 * np.log10(abs(S21))
    if interp_for_att:
        att_itp = itp.interp1d(freq, S21_dB, fill_value='extrapolate')
        S21_dB_lowT = att_itp(freq_vec) * lowT_factor
    else:
        b_guess = 0.5
        a_guess = - S21_dB[-1] / freq[-1] ** b_guess

        guess = ([a_guess, b_guess])
        bounds = (
            (0, 0),
            (10, np.inf)
        )

        qopt, qcov = curve_fit(pow_law, freq, S21_dB, guess, bounds=bounds)
        a_fit, b_fit = qopt
        S21_dB_fit = pow_law(freq, a_fit, b_fit)
        S21_dB_lowT = pow_law(f, lowT_factor * a_fit, b_fit)

    # plot attenuators
    fig, ax = plt.subplots()
    ax.grid(linestyle='--')
    leg = []
    plt.plot(freq / 1e9, S21_dB, 'r.')
    plt.plot(f / 1e9, S21_dB_lowT, 'g:')
    leg += ['Measurement', 'Low T']
    if not interp_for_att:
        plt.plot(freq / 1e9, S21_dB_fit, 'b-', Linewidth=2)
        leg += ['High T Fit']
        plt.title('attenuation=-a*f^b\na=%.3G, b=%.3G' % (a_fit, b_fit))
    plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim((S21_dB.min(), S21_dB.max()))
    plt.xlim((freq.min() / 1e9, freq.max() / 1e9))
    plt.xlabel('Freq(GHz)', fontsize='x-large')
    plt.ylabel('Attenuation(dB)', fontsize='x-large')
    plt.tick_params(axis='both', which='major', labelsize='x-large')
    plt.tight_layout()
    return S21_dB_lowT


fridge_att_config_list = [
    {
        'temperature': 4,
        # 'attenuation': ['22dB', 'ecco_filter2']
        # 'attenuation': ['2dB', 'ecco_filter2']
        'attenuation': ['20dB']
    },
    {
        'temperature': 0.06,
        'attenuation': ['20dB']
    },
    # {
    #     'temperature': 0.056,
    #     'attenuation': ['20dB'],
    # },
    {
        'temperature': 0.008,
        # 'attenuation': ['att1'],
        'attenuation': ['20dB'],
        # 'attenuation': ['att1', 'att2'],
    }
]

f_start = 100e6
f_stop = 20e9
fc = 6.48e9
LogX = False
if LogX:
    f = np.logspace(np.log10(f_start), np.log10(f_stop), num=101)
else:
    f = np.linspace(f_start, f_stop, 101)

f = f.reshape((len(f), 1))


att1_file = 'C:/SC Lab\\Labber\Data\\2019\\07\\Data_0709\\home made attenuator.hdf5'
S21_att1dB_lowT = component_att_from_measurement(att1_file, f, 0.74, True)

att2_file = 'C:/SC Lab\\Labber\Data\\2019\\07\Data_0731\homemade_attenuator2.hdf5'
S21_att2dB_lowT = component_att_from_measurement(att2_file, f, 0.48, True)

ecco_filter2_file = 'C:\\SC Lab\Projects\Fluxonium\data_process\Fluxonium032619/filter2.csv'
S21_ecco_filter2dB_lowT = component_att_from_measurement(ecco_filter2_file, f, 1, True)

components_dict = {
    'att1': - S21_att1dB_lowT,
    'att2': - S21_att2dB_lowT,
    'ecco_filter2': - S21_ecco_filter2dB_lowT,
}

temp = np.array([300])  # 300K, 4K, 100mK, 10mK
for i, stage in enumerate(fridge_att_config_list):
    temp = np.append(temp, stage['temperature'])
    components_list = stage['attenuation']
    tot_att = 0
    for comp in components_list:
        if comp.endswith('dB'):
            comp_att = float(comp[:-2])
            tot_att += comp_att * np.ones((len(f), 1))
        else:
            tot_att += components_dict[comp]
    # print(tot_att.shape)
    if i == 0:
        att_dB = tot_att
    else:
        att_dB = np.append(att_dB, tot_att, axis=1)
# att_dB = np.array([20, 20, 20])  # 4K, 100mK, 10mK
# temp = np.array([300, 4, 0.1, 0.02])  # 300K, 4K, 100mK, 10mK
att = 10 ** (att_dB / 10)

n_ph = np.zeros((len(temp), len(f)))
for i in range(len(temp)):
    n_ph[i, :] = thermal_photon(f[:, 0], temp[i])

# n_ph_att = np.array([n_ph[0] / att[0] / att[1] / att[2], n_ph[1] * (1 - 1 / att[0]) / att[1] / att[2],
#                      n_ph[2] * (1 - 1 / att[1]) / att[2], n_ph[3] * (1 - 1 / att[2])])
n_ph_att = attenuated_photon(att, n_ph)

# n_ph_c = (1 - 1 / att[-1]) * n_ph[-1, :] + 1 / att[-1] * (
#         (1 - 1 / att[-2]) * n_ph[-2, :] + 1 / att[-2] * ((1 - 1 / att[-3]) * n_ph[-3, :] + 1 / att[-3] * n_ph[0, :]))
n_ph_c = n_ph_att.sum(axis=0)
n_ph_c_itp = itp.interp1d(f[:, 0], n_ph_c)

Teff = effective_temperature(f[:, 0], n_ph_c)
Teff_c = effective_temperature(fc, n_ph_c_itp(fc))


if LogX:
    f_plot = f
    fc_plot = fc
else:
    f_plot = f / 1e9
    fc_plot = fc / 1e9
fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
for i in range(len(temp)):
    plt.plot(f_plot, n_ph_att[i, :], '--')
    leg += [str(temp[i]) + 'K']
plt.plot(f_plot, n_ph_c)
# plt.plot(f_plot, n_ph_eff)
leg += ['total photon number']
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
plt.axvline(x=fc_plot, color='g', linestyle=':')
plt.ylim((1e-2 * n_ph_c.min(), 10 * n_ph_c.max()))
plt.ylabel('n_eff', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.title('Thermal photon number that the cavity sees\n from all attenuators after attenuation')
ax.set_yscale('log')
if LogX:
    ax.set_xscale('log')
    plt.xlabel('Freq(Hz)', fontsize='x-large')
else:
    plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.tight_layout()


fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(f_plot, Teff)
# leg += ['outside cavity', 'effective T = %.3GK' % Teff]
# plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Effect T at %.3GGHz = %.3GK' % (fc / 1e9, Teff_c))
plt.axvline(x=fc_plot, color='g', linestyle=':')
# plt.ylim((1e-9, 1e4))
plt.ylabel('Effective temperature(K)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
ax.set_yscale('log')
if LogX:
    ax.set_xscale('log')
    plt.xlabel('Freq(Hz)', fontsize='x-large')
else:
    plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.tight_layout()

plt.show()
