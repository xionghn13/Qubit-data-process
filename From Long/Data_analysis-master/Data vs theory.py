# Plot intensity data

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from qutip import*
fig=plt.figure(figsize=(20, 10))
plt.rc('font', family='serif')
rc('text', usetex=False)
# plt.close("all")

ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)


#####################################################################################################################################################################################
#####################################################################################################################################################################################
# '''
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#Plasmon line scan
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_0to20mA_currentMode_qubit_n30dBm_cav_1dBm_avg50K_pulse25'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)
########################################################################
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_20to30mA_currentMode_qubit_n5dBm_cav_5dBm_avg20K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)
########################################################################
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_30to32mA_currentMode_qubit_n5dBm_cav_5dBm_avg20K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)
########################################################################
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_32to39mA_currentMode_qubit_n5dBm_cav_5dBm_avg20K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-5, vmax = 0)
########################################################################
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_39to50mA_currentMode_qubit_0dBm_cav_1dBm_avg20K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-6, vmax = -2)

########################################################################
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_46to48mA_currentMode_qubit_0dBm_cav_5dBm_avg20K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_cur = directory + '\\' + measurement + '_Current.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
Current = np.genfromtxt(path_cur, delimiter =',')*1e3
#Voltage = np.linspace(0,3,3000)
for idx in range(len(Current)-1):
    if (idx%10) == 0:
        f = Freq[idx]
        Z = RawData[idx:idx+11].transpose()
        I = Current[idx:idx+11]
        X, Y = np.meshgrid(I, f)
        plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin =-4, vmax = -1)
#Small scan
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_43to44mA_currentMode_qubit_2p5to3p2GHz_0dBm_cav5dBm_avg50K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_current = directory + '\\' + measurement + '_I.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
I = np.genfromtxt(path_current, delimiter =',')*1e3
Z = RawData.transpose()
X, Y = np.meshgrid(I,Freq)
plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin = -4 , vmax = -1)
########################################################################
#Small scan
directory = 'C:\Data\Fluxonium #10'
measurement = 'S21_43p15to43p85mA_currentMode_qubit_1p5to2p5GHz_0dBm_cav5dBm_avg50K_pulse(test)'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_current = directory + '\\' + measurement + '_I.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
I = np.genfromtxt(path_current, delimiter =',')*1e3
Z = RawData.transpose()
X, Y = np.meshgrid(I,Freq)
plt.pcolormesh(X, Y, Z, cmap=cm.GnBu_r, vmin = -4 , vmax = -1)
########################################################################
#Blue side band
directory = 'C:\Data\Fluxonium #10'
measurement = 'Two tune spectroscopy_YOKO 43p4to43p6mA_ qubit tone 10p5to11p2GHz_5dBm_Cav_10p304GHz_8dBm_pulse 34us duty2_avg5K'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_current = directory + '\\' + measurement + '_I.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
I = (np.genfromtxt(path_current, delimiter =',')-0.00003)*1e3
Z = RawData.transpose()
X, Y = np.meshgrid(I,Freq)
plt.pcolormesh(X, Y, Z, cmap=cm.RdBu_r, vmin = -5 , vmax = 5)
########################################################################
#Red side band
directory = 'C:\Data\Fluxonium #10'
measurement = 'Two tune spectroscopy_YOKO 43p4to43p6mA_ qubit tone 8p5to10p2GHz_5dBm_Cav_10p304GHz_8dBm_pulse 34us duty2_avg5K'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_current = directory + '\\' + measurement + '_I.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
I = (np.genfromtxt(path_current, delimiter =',')-0.00003)*1e3
Z = RawData.transpose()
X, Y = np.meshgrid(I,Freq)
plt.pcolormesh(X, Y, Z, cmap=cm.RdBu_r, vmin = -5 , vmax = 5)

# '''
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#Plotting data taken with new software
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_41to43mA_Qubit_3to4GHz_5dBm_Cav_10.3039GHz_8dBm_IF_0.05GHz_measTime_500ns_avg_50000'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1::]-0.04
freq = np.genfromtxt(path + '_FREQ.dat')
freq = freq[1::]
data = np.genfromtxt(path + '_PHASEMAG.dat')
phase = data[1::,0] #phase is recorded in rad
phase = phase#
mag = data[1::,0]

# plt.figure(1)
Z = np.zeros((len(current),len(freq)))
for idx in range(len(current)):
    temp = np.unwrap(phase[idx*len(freq):(idx+1)*len(freq)])
    Z[idx,:] = temp - np.average(temp)
Z = Z*180/(np.pi)
X,Y = np.meshgrid(current,freq[0:len(freq)/2+2])
Z1= Z.transpose()[0:len(freq)/2+2]
plt.pcolormesh(X,Y,Z1, cmap= 'GnBu_r', vmin = -4, vmax=-1, alpha = 1)

X,Y = np.meshgrid(current,freq[len(freq)/2+2:len(freq)-1])
Z2= Z.transpose()[len(freq)/2+2:len(freq)-1]
plt.pcolormesh(X,Y,Z2, cmap= 'GnBu_r', vmin = -4, vmax=-1, alpha = 1)
# '''
########################################################################
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_38.1to40mA_Qubit_3.5to5GHz_5dBm_Cav_10.3039GHz_8dBm_IF_0.05GHz_measTime_500ns_avg_25000'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1::]-0.04
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
Z = Z.transpose()[1:len(freq)-1]
X,Y = np.meshgrid(current,freq[1:len(freq)-1])
# plt.figure(1)
plt.pcolormesh(X,Y,Z, cmap= 'GnBu_r', vmin = -4, vmax=-1, alpha = 1)

#####################################################################
# high power scan
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_38to40mA_Qubit_3.5to5GHz_10dBm_Cav_10.3045GHz_5dBm_IF_0.05GHz_measTime_500ns_avg_50000'
path = directory + '\\' + measurement

# Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1::]-0.04
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
Z = Z.transpose()[1:len(freq)-1]
X,Y = np.meshgrid(current,freq[1:len(freq)-1])
# plt.figure(1)
# plt.pcolormesh(X,Y,Z, cmap= 'GnBu_r', vmin = -4, vmax=-1, alpha = 1)

#####################################################################
directory = 'G:\Projects\Fluxonium\Data\Fluxonium #10_New software'
measurement = 'Two_tone_spec_YOKO_38.56to38.66mA_Qubit_4.2to5.1GHz_-6dBm_Cav_10.3045GHz_5dBm_IF_0.05GHz_measTime_500ns_avg_20000'
path = directory + '\\' + measurement

#Read data
current = np.genfromtxt(path + '_CURR.dat')
current = current[1:-1] - 0.04
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
Z = Z.transpose()[1:len(freq)-1]
X,Y = np.meshgrid(current,freq[1:len(freq)-1])
plt.figure(1)
# plt.pcolormesh(X,Y,Z, cmap= 'Reds_r', vmin = -4, vmax=-0.5, alpha = 0.2)

#####################################################################################################################################################################################
#####################################################################################################################################################################################

#Define constants
e = 1.602e-19    #Fundamental charge
h = 6.62e-34    #Placnk's constant
phi_o = h/(2*e) #Flux quantum

"""
First section of the script attempts to plot the energies vs external flux
"""
#Hamiltonian definition
def Ho(N,E_l, E_c, E_j_sum, d, phi_squid, phi_ext):
    E_j1 = 0.5*E_j_sum*(1+d)
    E_j2 = 0.5*E_j_sum*(1-d)
    a = tensor(destroy(N))
    mass = 1.0/(8.0*E_c)
    w = sqrt(8.0*E_c*E_l)
    phi = (a+a.dag())*(8*E_c/E_l)**(0.25)/np.sqrt(2)
    na = 1j*(a.dag()-a)*(E_l/(8*E_c))**(0.25)/np.sqrt(2)
    ope1 = 1j*(-phi + phi_ext)
    ope2 = 1j*(phi - phi_ext + phi_squid) # phi_squid and phi_ext here are the external phases, or normalized flux, = flux*2pi/phi_o
    H = 4.0*E_c*na**2 + 0.5*E_l*(phi)**2 - 0.5*E_j1*(ope1.expm() + (-ope1).expm()) - 0.5*E_j2*(ope2.expm() + (-ope2).expm())
    return H.eigenenergies()

def coupled_H(Na, E_l, E_c, E_j_sum, d, phi_squid, phi_ext, Nr, wr, g):
    E_j1 = 0.5*E_j_sum*(1 + d)
    E_j2 = 0.5*E_j_sum*(1 - d)
    a = tensor(destroy(Na), qeye(Nr))
    b = tensor(qeye(Na), destroy(Nr))
    phi = (a + a.dag())*(8.0*E_c/E_l)**(0.25)/np.sqrt(2.0)
    na = 1.0j*(a.dag() - a)*(E_l/(8*E_c))**(0.25)/np.sqrt(2.0)
    ope1 = 1.0j*(phi_ext - phi)
    ope2 = 1.0j*(phi + phi_squid - phi_ext)
    H_f = 4.0*E_c*na**2 + 0.5 * E_l*(phi)** 2 - 0.5*E_j1*(ope1.expm()+(-ope1).expm()) - 0.5*E_j2*(ope2.expm()+(-ope2).expm())
    H_r = wr*(b.dag()*b + 1.0/2)
    H_c = -2*g * na * (b.dag + b)
    H = H_f + H_r + H_c
    return H.eigenenergies()


def trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState):  
    B_field = current*B_coeff*1e-4  # in T, this depends on a seperate measurement 
    phi_squid = B_field*A_j # these are flux, not normalized
    phi_ext = B_field*A_c
    trans_energy = np.zeros((level_num-iState,len(phi_ext))) 
    for idx in range(len(phi_ext)):
        energies = Ho(N,E_l, E_c, E_j_sum, d, 2*np.pi*(phi_squid[idx]/phi_o - beta_squid), 2*np.pi*(phi_ext[idx]/phi_o - beta_ext)) #normalize the flux -> phase here    
        for level in range(iState+1,level_num):
            trans_energy[level-iState,idx]=energies[level]-energies[iState]    
    return trans_energy

def coupled_trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState, Nr, wr, g):
    B_field = current*B_coeff*1e-4  # in T, this depends on a seperate measurement
    phi_squid = B_field*A_j # these are flux, not normalized
    phi_ext = B_field*A_c
    trans_energy = np.zeros((level_num-iState,len(phi_ext)))
    for idx in range(len(phi_ext)):
        energies = coupled_H(N, E_l, E_c, E_j_sum, d, 2*np.pi*(phi_squid[idx]/phi_o - beta_squid), 2*np.pi*(phi_ext[idx]/phi_o - beta_ext), Nr, wr, g) #normalize the flux -> phase here
        for level in range(iState+1,level_num):
            trans_energy[level-iState,idx]=energies[level]-energies[iState]
    return trans_energy
    

########################################################################
#Fitting for bottom spectrum
N = 50
E_l=0.722729827116
E_c=0.552669197076
E_j_sum=17.61374383
A_j=4.76321410213e-12
A_c=1.50075181762e-10
d=0.125005274368
beta_squid=0.129912406349
beta_ext=0.356925557542

current = np.linspace(0.038,0.047,9000)
B_coeff = 60
level_num = 5

iState = 0
spectrum = trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState)
for idx in range(iState,3):
    line = plt.plot(current*1e3, spectrum[idx,:])  # transition from state (iState)
    plt.setp(line,linewidth=1.0, linestyle ='-', color = "black", alpha=0.2)
    # line = plt.plot(current, spectrum[idx,:]+10.304)  # transition from state (iState)
    # plt.setp(line,linewidth=2.0, linestyle ='--', color = "black", alpha=0.5)
    # line = plt.plot(current, -spectrum[idx,:]+10.304)  # transition from state (iState)
    # plt.setp(line,linewidth=2.0, linestyle ='--', color = "black", alpha=0.5)


# iState = 1
# spectrum = trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState)
# for idx in range(iState,level_num):
#     line = plt.plot(current*1e3, spectrum[idx-iState,:])  # transition from state (iState)
#     plt.setp(line,linewidth=1.0, linestyle ='--', color = "red", alpha=0.5)
#     line = plt.plot(current, spectrum[idx-iState,:]+10.304)  # transition from state (iState)
#     plt.setp(line,linewidth=2.0, linestyle ='--', color = "red", alpha=0.5)
#     line = plt.plot(current, -spectrum[idx-iState,:]+10.304)  # transition from state (iState)
#     plt.setp(line,linewidth=2.0, linestyle ='-.', color = "red", alpha=0.5)

'''
#Coupled Transition energy calculation and fitting for top spectrum here
N = 40
Nr = 10
E_l=0.735773762652
E_c=0.537375025825
E_j_sum=22.3
A_j=3.83424869313e-12
A_c=1.46689233147e-10
d=0.185865262485
beta_squid=-2.58488114861e-05
beta_ext=-0.0251115059548
B_coeff = 60
g=0.0845608058905
wr = 10.304
current = np.linspace(0.02,0.03,100)
level_num = 7

iState = 0
spectrum = coupled_trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState, Nr, wr, g)
for idx in range(iState,level_num):
    line = plt.plot(current*1e3, spectrum[idx,:])  # transition from state (iState)
    plt.setp(line,linewidth=1.0, linestyle ='-', color = "black", alpha=0.5)
    line = plt.plot(current*1e3, spectrum[idx,:]/2)  # transition from state (iState)
    plt.setp(line,linewidth=1.0, linestyle ='-.', color = "black", alpha=0.9)


iState = 1
spectrum = coupled_trans_energies(N, E_l, E_c, E_j_sum, d, A_j, A_c, B_coeff, beta_squid, beta_ext, level_num, current, iState, Nr, wr, g)
for idx in range(iState,level_num):
    line = plt.plot(current*1e3, spectrum[idx-iState,:])  # transition from state (iState)
    plt.setp(line,linewidth=1.0, linestyle ='--', color = "red", alpha=0.5)
'''
# plt.grid("on")
# plt.xlabel("YOKO I (mA)")
# plt.ylabel("Freq (GHz)")
#plt.title(measurement)
# plt.xlim([38.3,46.7])
# plt.ylim([2,5])
# plt.tick_params(labelsize=18)
# plt.xlim([38.546, 38.63])
# plt.ylim([4.2,5.1])
plt.xlim([38.3,46.7])
# plt.xlim([38.3,39.6])
plt.ylim([2,5])
plt.tick_params(labelsize=16)
# plt.colorbar()
def onclick(event):
    print '[%f, %f],'%(event.xdata, event.ydata)
cid = fig.canvas.mpl_connect('button_press_event', onclick)

click_data = np.array([[38.345726, 4.775000],
[38.323589, 4.782484],
[38.342809, 4.774688],
[38.364651, 4.766892],
[38.394355, 4.743503],
[38.456384, 4.708420],
[38.484341, 4.688929],
[38.514046, 4.665541],
[38.534140, 4.642152],
[38.555108, 4.618763],
[38.563844, 4.591476],
[38.583065, 4.427755],
[38.604032, 4.225052],
[38.625000, 4.010655],
[38.633737, 3.909304],
[38.653831, 3.694906],
[39.364987, 3.776767],
[39.372849, 3.870322],
[39.385081, 3.967775],
[39.392944, 4.065229],
[39.405175, 4.162682],
[39.413038, 4.232848],
[39.422648, 4.295218],
[39.434005, 4.326403],
[39.457594, 4.349792],
[39.478562, 4.357588],
[39.506519, 4.369283],
[39.541465, 4.377079],
[39.582527, 4.380977],
[39.617339, 4.381250],
[39.667419, 4.395313],
[39.743629, 4.390625],
[39.806774, 4.385938],
[39.887339, 4.385938],
[39.976613, 4.371875],
[40.024516, 4.367188],
[40.135565, 4.343750],
[40.196532, 4.315625],
[40.235726, 4.301563],
[40.359839, 4.250000],
[40.499194, 4.179688],
[40.597177, 4.114063],
[40.675565, 4.057813],
[40.727823, 4.010938],
[40.790968, 3.935938],
[40.836694, 3.856250],
[40.845403, 3.837500],
[41.583871, 3.050000],
[41.591935, 3.125000],
[41.604839, 3.200000],
[41.612903, 3.265625],
[41.622581, 3.326563],
[41.632258, 3.378125],
[41.640323, 3.415625],
[41.661290, 3.485938],
[41.682258, 3.537500],
[41.712903, 3.579688],
[41.743548, 3.607813],
[41.796774, 3.650000],
[41.835484, 3.664063],
[41.875806, 3.682813],
[41.924194, 3.696875],
[41.974194, 3.706250],
[42.022581, 3.715625],
[42.074194, 3.720313],
[42.132258, 3.725000],
[42.214516, 3.725000],
[42.308065, 3.701563],
[42.388710, 3.682813],
[42.459677, 3.664063],
[42.504839, 3.650000],
[42.572581, 3.607813],
[42.635484, 3.579688],
[42.737097, 3.509375],
[42.858065, 3.387500],
[42.943548, 3.284375],
[42.983871, 3.214063],
[42.973871, 3.237500],
[43.004032, 3.176563],
[43.030645, 3.110938],
[43.082097, 2.951563],
[43.124677, 2.773438],
[43.156613, 2.604688],
[43.192097, 2.328125],
[43.222258, 2.098438],
[43.234677, 2.023438],
[43.713710, 2.056250],
[43.735000, 2.182812],
[43.754516, 2.318750],
[43.784677, 2.478125],
[43.804194, 2.567188],
[43.834355, 2.679688],
[43.884032, 2.825000],
[43.912419, 2.885938],
[43.955000, 2.960938],
[44.004677, 3.026563],
[44.073871, 3.092188],
[44.114677, 3.125000],
[44.196290, 3.162500],
[44.272581, 3.195312],
[44.359516, 3.214063],
[44.465968, 3.218750],
[44.531613, 3.209375],
[44.595484, 3.200000],
[44.648710, 3.190625],
[44.744516, 3.143750],
[44.817258, 3.092188],
[44.882903, 3.045313],
[44.962742, 2.960938],
[45.024677, 2.890625],
[45.054839, 2.829688],
[45.105565, 2.735938],
[45.148065, 2.637500],
[45.952823, 2.576563],
[45.972016, 2.642188],
[46.011774, 2.721875],
[46.066613, 2.834375],
[46.294194, 3.101562],
[46.375081, 3.153125],
[46.435403, 3.171875],
[46.525887, 3.204688],
[46.628710, 3.218750]
])
# plt.plot(click_data[:,0],click_data[:,1], 'ro')
plt.show()


