# This program computes quality factor from reflection measurement
# from sympy import*
from scipy import *
import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
# import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
import ExtractDataFunc as edf

DataFolderName = '10092019_wg5 in 8.5GHz cavity (add coax atts, eccosorb ...)'
DataPath = 'C:/SC Lab\\Labber\\' + DataFolderName + '/2019/10\Data_1010/'
S11File = 'cavity Q_14.hdf5'


[Freq, RComplex] = edf.readVNAS11(DataPath + S11File)
print(Freq)
print(RComplex)
TruncateFreq = False

StartFreq = 7.7
EndFreq = 7.8
if TruncateFreq:
    FreqInd = (EndFreq >= Freq) == (Freq >= StartFreq)
    Freq = Freq[FreqInd]
    RComplex = RComplex[FreqInd]

# plt.close('all')

# directory = 'G:\Projects\Cavity and Wave Guide\Copper cavity\\2019_01\\7p5GHz_4mm2_in_'

# fname = 'Mag.csv'
# path = directory + fname
# data = np.genfromtxt(path, delimiter=',', skip_header=3)
# freq = data[:, 0] / 1e9
# mag_data = data[:, 1]
#
# fname = 'Phase.csv'
# path = directory + fname
# data = np.genfromtxt(path, delimiter=',', skip_header=3)
# phase_data = data[:, 1]

mag_data = np.abs(RComplex)
phase_data = np.unwrap(np.angle(RComplex))
freq = Freq
#
# directory = 'G:\Projects\Cavity and Wave Guide\Copper cavity\\2018_11_27_4pieces\test'
# fname = 'test'
# path = directory + '\\' + fname
#
# phase = np.genfromtxt(path+"_PHASEMAG.csv")[1:,0]
# mag_data = np.genfromtxt(path+"_PHASEMAG.csv")[1:,1]/0.003
# freq = np.genfromtxt(path+"_FREQ.csv")[1:]
# phase = np.unwrap(phase)*180/np.pi
# offset = (phase[-1] - phase[0])/(freq[-1] - freq[0])*freq
# phase = phase - offset
# phase_data = phase - np.mean(phase)
# phase_data = phase_data - phase_data[0]

# Find resonant frequency and offsets
index = 0
fo = 0
amp_cal = 0
indexo = 0
# phase_cal = (phase_data[0] + phase_data[len(phase_data) - 1]) / 2.0
phase_cal = 0
print(phase_cal)
# for amp in mag_data:
#     if amp_cal > amp:
#         amp_cal = amp
#         indexo = index
#     index = index + 1
# for amp in phase_data:
#    if abs(amp - phase_cal) < 2:
#        indexo=index
#    index = index+1
indexo = np.argmin(mag_data)
fo = freq[indexo]

print(fo)
# offset=phase_data[indexo]
# phase_data = phase_data - offset
# mag_data = mag_data - mag_data[0]

# Convert to appropriate units
# mag_dataN = 10.0 ** (mag_data / 20.0)
# phase_dataN = phase_data * pi / 180.0
mag_dataN = mag_data
phase_dataN = phase_data
# If real and imaginary data are not available, compute them here. Skip if data are available
real_data = mag_dataN * cos(phase_dataN)
imag_data = mag_dataN * sin(phase_dataN)
all_data = np.array([real_data, imag_data]).flatten()

f, axarr = plt.subplots(2, sharex=True)

# Theoretical model

def S_11(f, q_ext1, q_ext2, q_int, tau1, tau2, tau3):
    S = (2j * (f - fo) / fo - (q_ext1 + 1j * q_ext2) ** (-1) + q_int ** (-1)) / (
            2j * (f - fo) / fo + (q_ext1 + 1j * q_ext2) ** (-1) + q_int ** (-1))
    S = S * tau1 * exp(1j * (tau2 + f * tau3))
    return np.array([S.real, S.imag]).flatten()


# Make a guess here for the values of Q_ext, Q_int,
# tau1 (magnitude offset), tau2 (linear phase offset), tau3 (freq dependent phase offset)

guess = ([5e3, -1e3, 1e4, 0, 0, 0])

qopt, qcov = curve_fit(S_11, freq, all_data, guess)
q_external1 = qopt[0]
q_external2 = qopt[1]
q_internal = qopt[2]
tau_final1 = qopt[3]
tau_final2 = qopt[4]
tau_final3 = qopt[5]

# Check

S_final = (2j * (freq - fo) / fo - (q_external1 + 1j * q_external2) ** (-1) + q_internal ** (-1)) / (
        2j * (freq - fo) / fo + (q_external1 + 1j * q_external2) ** (-1) + q_internal ** (-1))
S_final = S_final * tau_final1 * exp(1j * (tau_final2 + freq * tau_final3))
mag_final = abs(S_final)
phase_final = angle(S_final)
phase_final = unwrap(phase_final)
real_final = S_final.real
imag_final = S_final.imag

# Plotting
axarr[0].plot(freq, mag_dataN, '-')
axarr[1].plot(freq, phase_dataN, '-')

axarr[0].plot(freq, mag_final, '--', alpha=0.5)
axarr[1].plot(freq, phase_final, '--', alpha=0.5)

# Specify title and axes
# axarr[0].set_title('06/04/2015: 52mK, -80dBm \n' + str(qopt4))
# axarr[0].set_ylabel('dB')
# axarr[1].set_ylabel('degrees')
axarr[1].set_xlabel('Frequency(GHz)')

# axarr[0].set_title(
#     'Q_ext= ' + str(qopt[0]) + ", Q_int= " + str(qopt[2]) + ", f_o= " + str(fo))  # + " , error= " + str(q_ext_error))
axarr[0].set_title(
    'Q_ext= %.3G, Q_int= %.3G, f_o= %.3G' % (qopt[0], qopt[2], fo))

plt.tight_layout()

plt.show()