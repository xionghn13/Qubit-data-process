import numpy as np
from matplotlib import pyplot as plt


#File path
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_38to40mA_Qubit_3.5to5GHz_10dBm_Cav_10.3045GHz_5dBm_IF_0.05GHz_measTime_500ns_avg_50000'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1::]
freq = np.genfromtxt(path + '_FREQ.dat')
freq = freq[1::]
data = np.genfromtxt(path + '_PHASEMAG.dat')
phase = data[1::,0] #phase is recorded in rad
phase = phase#
mag = data[1::,0]

Z = np.zeros((len(current),len(freq)))
for idx in range(len(current)):
    temp = np.unwrap(phase[idx*len(freq):(idx+1)*len(freq)])
    Z[idx,:] = temp - np.average(temp)
Z = Z*180/(np.pi)
X,Y = np.meshgrid(current,freq)
plt.figure(1)
plt.pcolormesh(X,Y,Z.transpose(), cmap= 'GnBu_r', vmin = -5, vmax=0)
plt.colorbar()
plt.xlabel('Current (mA)')
plt.ylabel('Qubit freq (GHz)')

plt.figure(1)
for idx in range(len(current)):
    Z[idx,:] = mag[idx*len(freq):(idx+1)*len(freq)]
    # Z[idx,:] = temp - np.average(temp)

plt.figure(2)
plt.pcolormesh(X,Y,Z.transpose(), cmap= 'GnBu_r')#, vmin = -5, vmax=0)
plt.colorbar()
plt.xlabel('Current (mA)')
plt.ylabel('Qubit freq (GHz)')
plt.show()
plt.ion()
# print phase
