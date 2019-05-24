from QubitDataProcessPackages import *
I0 = 2.199
I_period = 2.469 * 2

FreqSingle = np.array(
    [6.6199, 4.5596, 4.2657, 3.9635, 3.6598, 3.3561, 1.8472, 1.5516, 1.2614, 0.9822, 0.7298, 0.54, 0.4864])
CurrentSingle = np.array([2.5, 3.2, 3.3, 3.4, 3.5, 3.6, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.674])
T1Single = np.array([8.81, 19.3, 11.9, 29.6, 30.2, 27.3, 25.2, 23, 28, 41.6, 42.3, 48.1, 11.5])
T1ErrSingle = np.array([0.17, 0.49, 0.57, 0.66, 0.37, 0.5, 1.1, 1.9, 3.9, 4.9, 1.2, 1.2, 2.4])
T2Single = np.array([3.18, 4.85, 3.42, 3.8, 3.79, 3.75, 3.27, 4.11, 3.34, 4.29, 5.66, 5.33, 7.73])
T2ErrSingle = np.array([0.11, 0.42, 0.25, 0.26, 0.26, 0.35, 0.26, 0.48, 0.3, 0.54, 0.64, 0.27, 0.95])
FluxSingle = (CurrentSingle - I0) / I_period

FreqRepeated = np.array(
    [7.18915, 7.1005, 6.8807, 6.338, 6.0534, 5.7595, 5.462, 5.169, 4.87, 3.0522, 2.7488, 2.4463, 2.1462])
CurrentRepeated = np.array([2.22, 2.3, 2.4, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.7, 3.8, 3.9, 4])
T1Repeated = np.array([2.87, 3.74, 1.61, 7.37, 17.3, 21.2, 21.4, 17.1, 11.5, 8.98, 34.9, 17.1, 45.4])
T1ErrRepeated = np.array([0.22, 0.4, 0.068, 1.1, 1.3, 4.6, 0.78, 0.82, 0.68, 1, 2.5, 3.1, 5.3])
T2Repeated = np.array([4.34, 3.62, 2.05, 3.01, 3.77, 3.56, 3.85, 3.17, 2.97, 3.24, 3.6, 2.72, 3.46])
T2ErrRepeated = np.array([0.32, 0.41, 0.075, 0.76, 0.072, 0.17, 0.088, 0.11, 0.075, 0.4, 0.07, 0.32, 0.077])
FluxRepeated = (CurrentRepeated - I0) / I_period

FreqMerge = np.concatenate((FreqRepeated, FreqSingle))
FluxMerge = np.concatenate((FluxRepeated, FluxSingle))
T1Merge = np.concatenate((T1Repeated, T1Single))
T2Merge = np.concatenate((T2Repeated, T2Single))
SortInd = FreqMerge.argsort()
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