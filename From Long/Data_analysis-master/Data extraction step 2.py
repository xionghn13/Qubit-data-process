# Plot intensity data

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from scipy import interpolate

rc('text', usetex=False)

#Enter directory and name of measurement
directory = 'G:\Projects\Fluxonium & qubits\Data\\2016_01\\13'
measurement = 'S21_Phase2tones_ZNB_0&n40_BW100Hz_SMB_n20dBm_6to8GHz_YOKO_n7ton10V'
path_data = directory + '\\' + measurement + '_Phase_diff.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
V = np.genfromtxt(path_vol, delimiter =',')

Z = RawData.transpose()
          
#Optional: calculate differential
Z_diff = np.diff(Z.transpose()) 
Z_diff = Z_diff.transpose()

#Plot the intensity map
X, Y = np.meshgrid(V,Freq)
fig=plt.figure()
#plt.pcolormesh(X, Y, Z, cmap=cm.BuGn, vmin =-4 , vmax =1)
plt.title(measurement)
plt.xlabel("Voltage (V)")
plt.ylabel("Frequency (GHz)")
#plt.colorbar()

##################################Approximation line#################################
A = np.array([
[-9.977319, 7.857741],
[-9.841230, 7.852510],
[-9.561492, 7.821130],
[-9.221270, 7.721757],
[-8.744960, 7.559623],
[-8.283770, 7.298117],
[-8.011593, 7.115063],
[-7.792339, 6.984310],
[-7.610887, 6.832636],
[-7.474798, 6.748954],
[-7.300907, 6.607741],
[-7.202621, 6.513598],
[-7.127016, 6.414226],
[-7.058972, 6.278243],
[-7.028730, 6.231172]

])
x=A[:,0]
y=A[:,1]
#plt.plot(x,y,'bo')

################################Interpolation##################################

tck = interpolate.splrep(x, y, s=1)
f = interpolate.splev(V, tck, der=0)
#plt.plot(V,f,'r-.')

################################Find minimum###################################
seek_range = 0.2
Freq_min = np.zeros(len(V))
for idy in range(len(Z[0,:])):
    min_val = 0.0
    for idx in range(len(Z[:,0])):
        if (Freq[idx] < (f[idy] + seek_range)) and (Freq[idx] > (f[idy] - seek_range)):
            if (Z[idx,idy] < min_val):
                min_val = Z[idx,idy]
                idxo=idx
    Freq_min[idy] = Freq[idxo]

plt.plot(V,Freq_min,'k--',linewidth = '2')
#############################################
directory = 'G:\Projects\Fluxonium & qubits\Data\\2016_01\\12'
measurement = 'S21_Phase_ZNB_0dBm&n40dB_YOKO_n11p4ton11p9'
path_data = directory + '\\' + measurement + '_Phase_diff.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
V = np.genfromtxt(path_vol, delimiter =',')
Z = RawData.transpose()          
X, Y = np.meshgrid(V,Freq)
plt.pcolormesh(X, Y, Z)

plt.show()
