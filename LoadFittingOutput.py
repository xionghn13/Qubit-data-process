import pickle
import matplotlib.pyplot as plt
import numpy as np

with open ('FittingOutput', 'rb') as fp:
    itemlist = pickle.load(fp)


[CounterArray, OptMatrixList, ErrMatrixList] = itemlist

fig, ax = plt.subplots()
ax.grid(linestyle='--')
plotInd = 3
ind = 1
Tpi_us = OptMatrixList[ind][plotInd, :] / 1000
Tpi_us_err = ErrMatrixList[ind][plotInd, :] / 1000
delta_f_kHz = 1000 / 2 / Tpi_us
delta_f_kHz_err = Tpi_us_err / 2 / Tpi_us ** 2 * 1000
avgList = []
stdList = []
avg_std_display = []
avgList += [np.nanmean(delta_f_kHz)]
stdList += [np.nanstd(delta_f_kHz)]
ax.errorbar(CounterArray[:], delta_f_kHz, yerr=delta_f_kHz_err, fmt='o')
avg_std_display += [avgList[0], stdList[0]]
tit = '$\delta f$=%.3G$\pm$%.2GkHz'
plt.title(tit % tuple(avg_std_display))
plt.plot(CounterArray[:],
         delta_f_kHz, '--')
plt.xlabel('Trial #', fontsize='x-large')
plt.ylabel('$\delta f$(kHz)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

t2leg = ['T2Ramsey', 'T2echo']
t2title = 'T2Ramsey=%.3G$\pm$%.2Gus, T2echo=%.3G$\pm$%.2Gus'
fig, ax = plt.subplots()
ax.grid(linestyle='--')
plotInd = 1
avgList = []
stdList = []
avg_std_display = []
num_meas = len(OptMatrixList)
for ind in range(num_meas):
    avgList += [np.nanmean(OptMatrixList[ind][plotInd, :] / 1000)]
    stdList += [np.nanstd(OptMatrixList[ind][plotInd, :] / 1000)]
    ax.errorbar(CounterArray[:],
                OptMatrixList[ind][plotInd, :] / 1000,
                yerr=ErrMatrixList[ind][plotInd, :] / 1000, fmt='o')
    avg_std_display += [avgList[ind], stdList[ind]]
T1Name = 'T1'
plt.legend([T1Name] + t2leg)
tit = T1Name + '=%.3G$\pm$%.2Gus, ' + t2title
plt.title(tit % tuple(avg_std_display))
num_meas = len(OptMatrixList)
for ind in range(num_meas):
    plt.plot(CounterArray,
             OptMatrixList[ind][plotInd, :] / 1000, '--')

plt.show()