from QubitDataProcessPackages import *


def thermal_photon(frequency, temperature):
    photon = 1 / (np.exp(h * frequency / k_B / temperature) - 1)
    return photon


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


att_dB = np.array([20, 40, 20])  # 4K, 100mK, 10mK
temp = np.array([300, 4, 0.1, 0.06])  # 300K, 4K, 100mK, 10mK
att = 10 ** (att_dB / 10)

f = np.logspace(6, 12, num=101)
# f = f.reshape((1, len(f)))
fc = 11e9

n_ph = np.zeros((len(temp), len(f)))
for i in range(len(temp)):
    n_ph[i, :] = thermal_photon(f, temp[i])

# n_ph_att = np.array([n_ph[0] / att[0] / att[1] / att[2], n_ph[1] * (1 - 1 / att[0]) / att[1] / att[2],
#                      n_ph[2] * (1 - 1 / att[1]) / att[2], n_ph[3] * (1 - 1 / att[2])])
n_ph_att = attenuated_photon(att, n_ph)

n_ph_c = (1 - 1 / att[-1]) * n_ph[-1, :] + 1 / att[-1] * (
        (1 - 1 / att[-2]) * n_ph[-2, :] + 1 / att[-2] * ((1 - 1 / att[-3]) * n_ph[-3, :] + 1 / att[-3] * n_ph[0, :]))
# n_ph_c = n_ph_att.sum(axis=0)
n_ph_c_itp = itp.interp1d(f, n_ph_c)
# print(f)
Teff = h * fc / k_B / np.log(1 + 1 / n_ph_c_itp(fc))
n_ph_eff = thermal_photon(f, Teff)

fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
for i in range(len(temp)):
    plt.plot(f, n_ph_att[i, :], '--')
    leg += [str(temp[i]) + 'K']
plt.plot(f, n_ph_c)
plt.plot(f, n_ph_eff)
leg += ['outside cavity', 'effective T = %.3GK' % Teff]
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
plt.axvline(x=fc, color='g', linestyle=':')
plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(Hz)', fontsize='x-large')
plt.ylabel('thermal photon', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()
