import numpy as np
import matplotlib.pyplot as plt
folder = 'C:\\Projects\\Fluxonium\\data_process\\Fluxonium032619/'
file = 'driving2_2.csv'
data = np.genfromtxt(folder + file, skip_header=3, delimiter=';')
freq = data[:, 0]
S21 = data[:, 1] + 1j * data[:, 2]
S21dB = 20 * np.log10(np.abs(S21))

fig, ax = plt.subplots()
ax.grid(linestyle='--')
leg = []
plt.plot(freq / 1e9, S21dB, '-')
# leg += ['Measurement', 'Fit']
# plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5))
# plt.axvline(x=fc, color='g', linestyle=':')
# plt.title('attenuation=-a*f^b\na=%.3G, b=%.3G' % (a_att1_fit, b_att1_fit))
# plt.ylim((1e-9, 1e4))
plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.ylabel('S21(dB)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()
print(data)
print(freq)

# ------------------------------------------------------------------
