import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=False)
plt.close("all")

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline_qubit_n30dBm_cav_10dBm&n30dB_YOKO_0to1p5V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)


directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline_qubit_n27dBm_cav_10dBm&n30dB_YOKO_1p5to3V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline_qubit_n20dBm_cav_10dBm&n30dB_YOKO_3to4V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
# Voltage = np.linspace(3,4,1000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline_qubit_n10dBm_cav_10dBm&n30dB_YOKO_4to5'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
# Voltage = np.linspace(3,4,1000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline_qubit_n10dBm_cav_10dBm&n30dB_YOKO_4to5 cont'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
# Voltage = np.linspace(3,4,1000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21_2tone_plasmonline2_qubit_n5dBm_cav_1dBm&n30dB_YOKO_0to3V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'phase slip test4 cav 1dBm qubit n10dBm'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'phase slip test4 cav 1dBm qubit 0dBm p3'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'phase slip test4 cav 1dBm qubit n10dBm p2'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'phase slip test4 cav 1dBm qubit 0dBm'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit 0dBm YOKO 6to7V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit 3dBm YOKO 5p6to6V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit n10dBm YOKO 6p5to7V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit 0dBm YOKO 6p2to6p4V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit n10dBm YOKO 4p8to4p9V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse 3rd CD'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit n30dBm YOKO 0to3V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0, alpha =0.4)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse 3rd CD'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit 0dBm YOKO 4p5to5p5V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0, alpha =0.5)

directory = 'G:\Projects\Fluxonium & qubits\Data\Fluxonium #10\Pulse 3rd CD'
measurement = 'S21 2tones PlasmonLine cav 1dBm qubit 0dBm new YOKO 4p8to5p6V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Voltage = np.genfromtxt(path_vol, delimiter =',')
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Voltage)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        V = Voltage[idx:idx+11]
        X, Y = np.meshgrid(V, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0, alpha =0.4)


#####################################################################################
# qu_V = np.linspace(0,10,1000)
# f = 9.38*np.cos(0.24*qu_V) + 0.3
# g = 5-f
# plt.plot(qu_V,g, qu_V,f)
# plt.plot(qu_V, g-0.6, '--', qu_V, g+0.6,'--', qu_V, f-0.6,'--', qu_V, f+0.6,'--')
#####################################################################################
#####################################################################################                                                                                                
plt.grid()
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
#plt.colorbar()
#plt.xlabel("Voltage")
#plt.ylabel("Frequency (GHz)")                
plt.show()

# Z1 = RawData[0].transpose()
# Z2 = RawData[1].transpose()
# print np.vstack((Z1,Z2)).transpose()