from QubitDataProcessPackages import *


def thermal_photon(frequency, temperature):
    photon = 1 / (np.exp(h * frequency / k_B / temperature) - 1)
    return photon

def effective_temperature(frequency, photon_number):
    temperature_effective = h * frequency / k_B / np.log(1 + 1 / photon_number)
    return temperature_effective


def attenuated_photon(att_vec, n_ph_mat):
    num_f = n_ph_mat.shape[1]
    num_T = len(att_vec) + 1
    n_ph_attenuated = np.zeros((num_T, num_f))
    for i in range(num_T):
        if i == 0:
            n_ph_attenuated[i, :] = n_ph_mat[i, :] / np.prod(att_vec)
        else:
            n_ph_attenuated[i, :] = n_ph_mat[i, :] / np.prod(att_vec[i - 1:]) * (att_vec[i - 1] - 1)
    return n_ph_attenuated


def pow_law(f, a, b):
    t = - a * f ** b
    return t

# att_dB_list =
att_dB = np.array([20, 20, 20])  # 4K, 100mK, 10mK
temp = np.array([300, 4, 0.1, 0.02])  # 300K, 4K, 100mK, 10mK
att = 10 ** (att_dB / 10)

f = np.logspace(6, 12, num=101)
# f = f.reshape((1, len(f)))
fc = 12e9
att1 = 'C:/Users\\admin\Labber\Data\\2019\\07\\Data_0709\\home made attenuator.hdf5'
[freq_att1, S21_att1] = edf.readVNAS21(att1)
freq_att1 *= 1e9
S21_att1dB = 20 * np.log10(abs(S21_att1))
b_guess = 0.5
a_guess = - S21_att1dB[-1] / freq_att1[-1] ** b_guess

guess = ([a_guess, b_guess])
bounds = (
    (0, 0),
    (10, np.inf)
)

qopt, qcov = curve_fit(pow_law, freq_att1, S21_att1dB, guess, bounds=bounds)
a_att1_fit, b_att1_fit = qopt
S21_att1dB_fit = pow_law(freq_att1, a_att1_fit, b_att1_fit)


att2 = 'C:/Users\\admin\Labber\Data\\2019\\07\Data_0731\homemade_attenuator2.hdf5'
[freq_att2, S21_att2] = edf.readVNA4portS21(att2)
freq_att2 *= 1e9
# att2_itp = itp.interp1d(freq_att1, S21_att1)
S21_att2dB = 20 * np.log10(abs(S21_att2))
b_guess = 0.5
a_guess = - S21_att2dB[-1] / freq_att2[-1] ** b_guess

guess = ([a_guess, b_guess])
bounds = (
    (0, 0),
    (10, np.inf)
)

print(freq_att1.shape)
qopt, qcov = curve_fit(pow_law, freq_att2, S21_att2dB, guess, bounds=bounds)
a_att2_fit, b_att2_fit = qopt
S21_att2dB_fit = pow_law(freq_att2, a_att2_fit, b_att2_fit)

n_ph = np.zeros((len(temp), len(f)))
for i in range(len(temp)):
    n_ph[i, :] = thermal_photon(f, temp[i])

# n_ph_att = np.array([n_ph[0] / att[0] / att[1] / att[2], n_ph[1] * (1 - 1 / att[0]) / att[1] / att[2],
#                      n_ph[2] * (1 - 1 / att[1]) / att[2], n_ph[3] * (1 - 1 / att[2])])
n_ph_att = attenuated_photon(att, n_ph)

# n_ph_c = (1 - 1 / att[-1]) * n_ph[-1, :] + 1 / att[-1] * (
#         (1 - 1 / att[-2]) * n_ph[-2, :] + 1 / att[-2] * ((1 - 1 / att[-3]) * n_ph[-3, :] + 1 / att[-3] * n_ph[0, :]))
n_ph_c = n_ph_att.sum(axis=0)
n_ph_c_itp = itp.interp1d(f, n_ph_c)
# print(f)

Teff = effective_temperature(f, n_ph_c)
Teff_c = effective_temperature(fc, n_ph_c_itp(fc))

# plot attenuators
fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(freq_att1 / 1e9, S21_att1dB, '.')
plt.plot(freq_att1 / 1e9, S21_att1dB_fit, '-')
# plt.plot(f, n_ph_eff)
leg += ['Measurement', 'Fit']
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
# plt.axvline(x=fc, color='g', linestyle=':')
plt.title('attenuation=-a*f^b\na=%.3G, b=%.3G' % (a_att1_fit, b_att1_fit))
# plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.ylabel('Attenuation(dB)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
# ax.set_yscale('log')
# ax.set_xscale('log')

fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(freq_att2 / 1e9, S21_att2dB, '.')
plt.plot(freq_att2 / 1e9, S21_att2dB_fit, '-')
# plt.plot(f, n_ph_eff)
leg += ['Measurement', 'Fit']
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
# plt.axvline(x=fc, color='g', linestyle=':')
plt.title('attenuation=-a*f^b\na=%.3G, b=%.3G' % (a_att2_fit, b_att2_fit))
# plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.ylabel('Attenuation(dB)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
# ax.set_yscale('log')
# ax.set_xscale('log')


fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
for i in range(len(temp)):
    plt.plot(f, n_ph_att[i, :], '--')
    leg += [str(temp[i]) + 'K']
plt.plot(f, n_ph_c)
# plt.plot(f, n_ph_eff)
leg += ['outside cavity']
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
plt.axvline(x=fc, color='g', linestyle=':')
plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(Hz)', fontsize='x-large')
plt.ylabel('Thermal photon', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_yscale('log')
ax.set_xscale('log')

fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(f, Teff)
# leg += ['outside cavity', 'effective T = %.3GK' % Teff]
# plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Effect T at %.3GGHz = %.3GK' % (fc / 1e9, Teff_c))
plt.axvline(x=fc, color='g', linestyle=':')
# plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(Hz)', fontsize='x-large')
plt.ylabel('Effective temperature(K)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()
