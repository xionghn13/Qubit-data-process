import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings
warnings.simplefilter('error', OptimizeWarning)
plt.figure(figsize=(5, 5))
plt.rc('font', family='serif')

#Fitting function
def T1_func(x,a,b,c,d):
    return a*np.exp(-(x-c)/b) + d

#Directory for reading data
directory = "G:\Projects\Fluxonium\Data\Fluxonium #10_python code by Jon\T1_T2_(0to1)_YOKO_43p75to45p4mA\T1_T2_YOKO_45p4mA\T1s"
measurement = 'generic_sequence_84'

#Information to save data
current_uncorr = 45.4
current = 45.361
freq = 1.558
fname = 'T1_T1err_'+str(current_uncorr)+'mA.csv'
#Optional: print flux here to double check
print str(current_uncorr) +'mA'

#Read data and fit here
path = directory + '\\' + measurement
time = np.genfromtxt(path+'_time0.csv')
phase_avg = np.zeros(len(time))
T1_guess = 100e-6 #Provide guess here
count = 0
T1_array = []
T1_err_array = []
for idx in range(0,10): #can change the range and indices here
    phase = np.genfromtxt(path + '_phase' + str(idx)+'.csv')
    phase = -phase*180/np.pi
    phase = phase - np.min(phase)
    guessA = np.max(phase)-np.min(phase)
    guess =[guessA, T1_guess, 0, 0]
    if guessA < 1: #Change as necessary. Observation of the raw data would suffice in determining the threshold.
        continue
    try:
        popt,pcov = curve_fit(T1_func, time*1e-9, phase, guess)
    except RuntimeError:
        continue
    except OptimizeWarning:
        continue
    a,b,c,d = popt
    perr = np.sqrt(abs(np.diag(pcov)))
    T1 = b*1e6
    T1_err = perr[1]*1e6
    T1_array = np.append(T1_array, T1)
    T1_err_array = np.append(T1_err_array, T1_err)
    phase_avg = phase_avg + phase
    count = count + 1
    #Plot raw data
    print T1, ',', T1_err
    plt.plot(time,phase)

###################################################################
#Save file
directory = 'G:\Projects\Fluxonium\Data\Summary of T1_T2_vs flux_Fluxonium#10\Corrected flux\Individual fits'
data = np.zeros((count, 5))
data[0:,0] = np.ones(count)*current
data[0:,1] = np.ones(count)*freq
data[0:,4] = np.ones(count)*current_uncorr
data[0:,2] = T1_array
data[0:,3] = T1_err_array
path = directory + '\\' + fname
np.savetxt(path, data, delimiter =',', header='YOKO, YOKO_uncorr, T1, T1_err')
#####################################################################
############################AVG PHASE################################
#####################################################################
#For comparision and conclusion purposes, calculate the average phase and fit it here
phase_avg = phase_avg / count
guessA = np.max(phase_avg)-np.min(phase_avg)
guess =[guessA, T1_guess, 0, 0]
popt,pcov = curve_fit(T1_func, time*1e-9, phase_avg, guess)
a,b,c,d = popt
perr = np.sqrt(abs(np.diag(pcov)))
T1 = b*1e6
T1_err = perr[1]*1e6
print "Average: ", T1,',', T1_err
plt.show()