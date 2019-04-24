import numpy as np
from matplotlib import pyplot as plt


#File path
directory = 'D:\Data\Fluxonium #10_New software'
measurement = 'onetone_phase_zoomed.dat'
path = directory + '\\' + measurement

data = np.genfromtxt(path)
data = data[1::,:]
freq = data[:,0]
phase = np.unwrap(data[:,1])*180/np.pi
phase = phase - (phase[-1] - phase[0])/(freq[-1] - freq[0])*freq
mag = data[:,2]
mag=20*np.log10(np.sqrt(mag))

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(freq, phase)
axarr[0].set_title('Phase')
axarr[1].plot(freq,mag)
axarr[1].set_title('Mag')

plt.show()