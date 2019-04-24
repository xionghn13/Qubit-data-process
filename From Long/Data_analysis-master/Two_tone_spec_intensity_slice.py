import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

#File path
directory = 'D:\Data\Fluxonium #10_New software'
measurement = 'Two tone spec_YOKO41mA_46mA_Qubit2 to 4GHz 5dBm_Cav 10p304GHz 5dBm'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1::]
freq = np.genfromtxt(path + '_FREQ.dat')
freq = freq[1::]
data = np.genfromtxt(path + '_PHASEMAG.dat')
phase = data[1::,0] #phase is recorded in rad
phase = phase
mag = data[1::,0]
# plt.figure(1)
Z = np.zeros((len(current),len(freq)))
for idx in range(len(current)):
    temp = np.unwrap(phase[idx*len(freq):(idx+1)*len(freq)])
    Z[idx,:] = temp - np.average(temp)

Z = Z*180/np.pi #Convert to degrees
#Default plot at index 0
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
l, =plt.plot(freq,Z[0])
ax.set_xlabel('Freq (GHz)')
ax.set_ylabel('Phase (deg)')

#Slider defined here
axcolor = 'lightgoldenrodyellow'
axFlux = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
sFlux = Slider(axFlux, 'Current', current[0] , current[-1], valinit = current[0], valfmt='%0.00f')

def update(flux):
    flux = sFlux.val
    idx = current.tolist().index(flux)
    l.set_ydata(Z[idx])

sFlux.on_changed(update)

plt.show()
