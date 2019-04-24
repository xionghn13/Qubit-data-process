from matplotlib import pyplot as plt
import numpy as np
import h5py
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)
warnings.simplefilter("error", RuntimeWarning)

def func(x, a, b, c, d):
    return a*np.exp(-(x-c)/b) + d

#Path for saving file
directory = 'G:\Projects\Fluxonium\Data\Summary of T1_T2_vs flux_Fluxonium#10\Automation code'
fname = 'T1_rabi_44to45mA.csv'
path_s = directory + '\\' + fname

#Path to read data
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'T1_auto3.h5'
path = directory + '\\' + measurement
# time_fname = 'T1_delays_35to36p1mA.txt'
# path_t = directory + '\\' + time_fname
t_inc = 10e3 #unit in ns
t1_guess = 30e-6 #Unit in s

#Read data and fit
time = np.linspace(0, 20*t_inc, 20)
with h5py.File(path,'r') as hf:
    print('List of arrays in this file: \n', hf.keys())
    count = np.array(hf.get('count'))
    # print count
    phase_raw = hf.get('demod_Phase0')
    freq_raw = np.array(hf.get('freq'))
    flux_raw = np.array(hf.get('flux'))
    total_count = len(flux_raw)
    T1_array = []
    T1_err_array = []
    freq = []
    flux = []
    RabiA = []
    for idx in range(total_count):
        phase = -phase_raw[idx, 0]
        phase = np.unwrap(phase)
        phase = phase - np.min(phase)
        phase = phase*180/np.pi
        guessA= np.max(phase) - np.min(phase)
        # plt.plot(phase)
        # print freq_raw[idx]
        if guessA < 0.5 or guessA > 10:
            continue
        guess = [guessA, t1_guess, 0, 0]
        try:
            popt, pcov = curve_fit(func, time*1e-9, phase, guess)
        except RuntimeError:
            # print "Cannot fit entry " + str(idx)
            continue
        except OptimizeWarning:
            # print "Doesn't fit well entry " + str(idx)
            continue
        except RuntimeWarning:
            # print "Doesn't fit well entry " + str(idx)
            continue
        a,b,c,d = popt #b is T1
        phase_fit = func(time*1e-9, a, b, c, d)
        perr = np.sqrt(abs(np.diag(pcov)))
        T1 = b*1e6 #unit is us
        T1_err = perr[1]*1e6
        if b*1e6 < t_inc/1e3 or b*1e6 > t_inc/1e3*19 or a < 0.5 or a > 10 or T1_err > T1*0.15:
            continue
        T1_err_array = np.append(T1_err_array, T1_err)
        T1_array = np.append(T1_array, T1) #T1 in us
        flux = np.append(flux, flux_raw[idx])
        freq = np.append(freq, freq_raw[idx])
        RabiA = np.append(RabiA, a)
        print T1, T1_err
        plt.plot(time,phase)
        # if b*1e6 > 20:
        #     plt.plot(time/1e3, phase, time/1e3, phase_fit)
        #     plt.title(str(b*1e6))
        #     plt.show()

###############################################################################
#Save file here
###############################################################################
data = np.zeros((len(flux),5))
data[:,0] = flux
data[:,1] = freq
data[:,2] = T1_array
data[:,3] = T1_err_array
data[:,4] = RabiA
print 'Total points taken ' + str(total_count)
print 'pts fitted well ' + str(len(T1_array))
#Check
np.savetxt(path_s, data, delimiter =',', header='YOKO, Freq, T1, T1_err, Amp')
# data = np.genfromtxt(path_s)
# plt.plot(data[:,0], data[:,2], 'ro')
plt.show()