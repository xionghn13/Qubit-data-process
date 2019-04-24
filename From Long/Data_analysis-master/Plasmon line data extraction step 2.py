from scipy.optimize import curve_fit
import numpy as np
from qutip import*
from matplotlib import pyplot as plt
import matplotlib.cm as cm


#####################################################################################
#####################################################################################
fig =plt.figure()
click_data = np.array([
[38.345726, 4.775000],
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
flux_points = click_data[:,0]*1e-3
freq_points = click_data[:,1]
plt.plot(flux_points, freq_points, 'r.')
#Define constants
e = 1.602e-19    #Fundamental charge
h = 6.62e-34    #Placnk's constant
phi_o = h/(2*e) #Flux quantum

"""
First section of the script attempts to plot the energies vs external flux
"""
#Hamiltonian definition
def Ho(N, E_l, E_c, E_j_sum, d, phi_squid, phi_ext):
    E_j1 = 0.5*(E_j_sum + d*E_j_sum)
    E_j2 = 0.5*(E_j_sum - d*E_j_sum)
    a = tensor(destroy(N))
    mass = 1.0/(8.0*E_c)
    w = sqrt(8.0*E_c*E_l)
    phi = (a+a.dag())*(8.0*E_c/E_l)**(0.25)/np.sqrt(2.0)
    na = 1j*(a.dag()-a)*(E_l/(8.0*E_c))**(0.25)/np.sqrt(2.0)
    ope1 = 1j*(-phi + phi_ext)
    ope2 = 1j*(phi - phi_ext + phi_squid)
    H = 4.0*E_c*na**2 + 0.5*E_l*(phi)**2.0 - 0.5*E_j1*(ope1.expm() + (-ope1).expm()) - 0.5*E_j2*(ope2.expm() + (-ope2).expm())
    return H.eigenenergies() 
    
def Ho_alt(N,E_l, E_c, E_j_sum, d, phi_squid, phi_ext):
    E_j = E_j_sum*np.cos(phi_squid/2.0)*np.sqrt(1+(d*np.tan(phi_squid/2.0))**2.0)
    theta = np.arctan(d*np.tan(phi_squid/2.0))
    a = tensor(destroy(N))
    mass = 1.0/(8.0*E_c)
    w = sqrt(8.0*E_c*E_l)
    phi = (a+a.dag())*(8*E_c/E_l)**(0.25)/np.sqrt(2)
    na = 1j*(a.dag()-a)*(E_l/(8*E_c))**(0.25)/np.sqrt(2)
    ope = 1j*(phi - phi_ext + theta + phi_squid/2.0) # phi_squid and phi_ext here are the external phases, or normalized flux, = flux*2pi/phi_o
    H = 4.0*E_c*na**2 + 0.5*E_l*(phi)**2 - 0.5*E_j*(ope.expm() + (-ope).expm())
    return H.eigenenergies()
        
def trans_energies(current, E_l, E_c, E_j_sum, A_j, A_c, d, offset_squid, offset_ext): 
    N=50    
    B_coeff = 60
    B_field = current*B_coeff*1e-4
    phi_squid = B_field*A_j
    phi_ext = B_field*A_c
    trans_energy = np.zeros(len(current))     
    for idx in range(len(current)):
        energies = Ho(N, E_l, E_c, E_j_sum, d, 2*np.pi*(phi_squid[idx]/phi_o - offset_squid), 2*np.pi*(phi_ext[idx]/phi_o - offset_ext))     
        trans_energy[idx]=energies[1]-energies[0]
    return trans_energy    
    
def plot_trans_energies(N, E_l, E_c, E_j_sum, d, phi_squid, phi_ext, offset_squid, offset_ext, level_num, current):  
    energy0 = np.empty((level_num,len(phi_ext)),dtype=float) 
    energy1 = np.empty((level_num,len(phi_ext)),dtype=float)
    energy2 = np.empty((level_num,len(phi_ext)),dtype=float) 
    for idx in range(len(phi_ext)):
        energies = Ho(N,E_l, E_c, E_j_sum, d, 2*np.pi*(phi_squid[idx]/phi_o- offset_squid), 2*np.pi*(phi_ext[idx]/phi_o- offset_ext))
        #1 photon
        for level in range(0,level_num):
            energy0[level,idx]=energies[level]-energies[0]
            #energy1[level,idx]=energies[level]-energies[1]
            #energy2[level,idx]=energies[level]-energies[2]
    for idx in range(0,level_num):    
        line = plt.plot(current, energy0[idx,:])
        plt.setp(line,linewidth=1.0, linestyle ='-', color = "black", alpha=0.7)
        #line = plt.plot(current, energy1[idx,:])
        #plt.setp(line,linewidth=1.0, linestyle ='-', color = "red", alpha=0.7)  
        
    return   

##########################################################################################################################################################################
##########################################################################################################################################################################
#Energy scale in GHz                                                            
E_l_guess = 0.7
E_c_guess = 0.6
E_j_sum_guess = 22

#Define external parameters
A_j_guess = 3.8e-12  #in m
A_c_guess = 150e-12
d_guess = 0
offset_squid_guess = 0
offset_ext_guess=0                                                  
guess = ([E_l_guess, E_c_guess, E_j_sum_guess, A_j_guess, A_c_guess, d_guess, offset_squid_guess, offset_ext_guess])                                                           
opt, cov = curve_fit(trans_energies, flux_points, freq_points, guess) 
                                                                     
################################################################################################################################################################################

E_l_fit, E_c_fit , E_j_fit , A_j_fit , A_c_fit, d_fit, offset_squid_fit, offset_ext_fit = opt

print 'E_l=' + str(E_l_fit) + ', E_c=' + str(E_c_fit) + ', E_j_sum=' + str(E_j_fit) + '\n' + 'A_j=' + str(A_j_fit) + ', A_c=' + str(A_c_fit) + ', d=' + str(d_fit) + \
', beta_squid=' + str(offset_squid_fit) + ', beta_ext=' + str(offset_ext_fit)
#Only plot the first excited state transition

current = np.linspace(0.03,0.05,1000)
# plt.plot (current, trans_energies(current, E_l_fit, E_c_fit, E_j_fit, A_j_fit, A_c_fit, d_fit, offset_squid_fit, offset_ext_fit), color = "black", alpha=0.7)
#plt.plot (current, trans_energies(current, E_l_guess, E_c_guess, E_j_sum_guess, A_j_guess, A_c_guess, d_guess, beta_squid_guess, beta_ext_guess)) 

#Plot other transition    
N=50
level_num = 5
B_coeff = 60
B_field = current*B_coeff*1e-4
phi_squid = B_field*A_j_fit
phi_ext = B_field*A_c_fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
fig = plt.figure(1)
plot_trans_energies(N, E_l_fit, E_c_fit, E_j_fit, d_fit, phi_squid, phi_ext, offset_squid_fit, offset_ext_fit, level_num, current)
# plt.ylim([0,12])
plt.grid()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
plt.show()

