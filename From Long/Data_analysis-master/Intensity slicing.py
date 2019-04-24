import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider

#Enter directory and name of measurement
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_41to43mA_Qubit_3to4GHz_5dBm_Cav_10.3039GHz_8dBm_IF_0.05GHz_measTime_500ns_avg_50000'
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

plt.figure(1)
Z = np.zeros((len(current),len(freq)))
for idx in range(len(current)):
    temp = np.unwrap(phase[idx*len(freq):(idx+1)*len(freq)])
    Z[idx,:] = temp - np.mean(temp)
Z = Z*180/(np.pi)
X,Y = np.meshgrid(current-0.041,freq)
# plt.figure(1)
# plt.pcolormesh(X,Y,Z.transpose(), cmap= 'GnBu_r')#, vmin = -5, vmax=0)
# plt.show()
# RawData = Z
Freq = freq
V = current

Z = Z.transpose()
          
#Optional: calculate differential
# Z_diff = np.diff(Z.transpose())
#Z_diff = Z_diff.transpose()
# Freq_diff = Freq[0:len(Freq)-1]

#Plot the data for V =0 here
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
l, =plt.plot(Freq,Z[:,0])
#Slider defined here
axcolor = 'lightgoldenrodyellow'
axVol = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
sVol = Slider(axVol, 'Voltage index', 0 , len(V)-1, valinit =0, valfmt='%0.0f')

def update(val):
    vol = sVol.val
    l.set_ydata(Z[:,vol])

sVol.on_changed(update)

plt.show()