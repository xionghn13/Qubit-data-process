from QubitDataProcessPackages import *

I0 = 2.351
I_period = 3.744 * 2

FreqSingle = np.array([109.6, 109.8, 110.2, 109.6, 109.6])
CurrentSingle = np.array([6.095, 21.042, 28.5006, 43.383, 50.787])
T1Single = np.array([27.3, 30.1, 28.8, 28.4, 31.9])
T1ErrSingle = np.array([3.6, 2.5, 0.41, 0.44, 0.95])
T2Single = np.array([6.97, 6.23, np.nan, 4.74, 5.3])
T2ErrSingle = np.array([0.82, 0.56, np.NaN, 0.28, 0.32])
FluxSingle = np.array([0.5, 2.5, 3.5, 5.5, 6.5])

FreqRepeated = np.array([109.5])
CurrentRepeated = np.array([13.5705])
T1Repeated = np.array([36.5])
T1ErrRepeated = np.array([1.9])
T2Repeated = np.array([6.63])
T2ErrRepeated = np.array([0.72])
FluxRepeated = np.array([1.5])

FreqMerge = np.concatenate((FreqRepeated, FreqSingle))
FluxMerge = np.concatenate((FluxRepeated, FluxSingle))
T1Merge = np.concatenate((T1Repeated, T1Single))
T2Merge = np.concatenate((T2Repeated, T2Single))
SortInd = FluxMerge.argsort()
FreqMerge = FreqMerge[SortInd]
FluxMerge = FluxMerge[SortInd]
T1Merge = T1Merge[SortInd]
T2Merge = T2Merge[SortInd]

fig, ax = plt.subplots()
ax.errorbar(FreqSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
ax.errorbar(FreqSingle, T2Single, yerr=T2ErrSingle, fmt='b^')

ax.errorbar(FreqRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
ax.errorbar(FreqRepeated, T2Repeated, yerr=T2ErrRepeated, fmt='r^')

plt.legend(['T1 - Single', 'T2echo - Single', 'T1 - Averge', 'T2echo - Averge'])

plt.plot(FreqMerge, T1Merge, ':')
plt.plot(FreqMerge, T2Merge, ':')

plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.ylabel('Decay time(us)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
ax.errorbar(FluxSingle, T1Single, yerr=T1ErrSingle, fmt='bo')
ax.errorbar(FluxSingle, T2Single, yerr=T2ErrSingle, fmt='b^')

ax.errorbar(FluxRepeated, T1Repeated, yerr=T1ErrRepeated, fmt='ro')
ax.errorbar(FluxRepeated, T2Repeated, yerr=T2ErrRepeated, fmt='r^')

plt.legend(['T1 - Single', 'T2echo - Single', 'T1 - Averge', 'T2echo - Averge'])

plt.plot(FluxMerge, T1Merge, ':')
plt.plot(FluxMerge, T2Merge, ':')
plt.xlabel('Flux/Phi_0', fontsize='x-large')
plt.ylabel('Decay time(us)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
