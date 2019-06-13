from QubitDataProcessPackages import *
Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium waveguide 1\T1/'
File = '121218_T1_YOKO_1.534mA_Cav7.0803GHz_-20dBm_Qubit1.059GHz_10dBm_PiPulse488ns_Count20_TimeStep5000.h5'
ParamList = File.split('_')
Count = int(ParamList[-2][5:])
TimeStep = float(ParamList[-1][8:-3])
f = h5py.File(Folder + File, 'r')
key_list = list(f.keys())

# Get the data
t1_amp = list(f[key_list[0]])[0][0]
t1_phase = list(f[key_list[1]])[0][0]
t1_complex = t1_amp * np.exp(1j * t1_phase)

t1_time = np.linspace(TimeStep, Count * TimeStep, Count)

plt.plot(t1_time, t1_phase)
plt.show()