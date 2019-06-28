from QubitDataProcessPackages import *
import Single_small_junction as ssj
from xlrd import open_workbook


def getRowData(sheet, rowInd, length=None):
    Row = sheet.row(rowInd)[1:]
    RowData = np.array([])
    if length == None:
        for ele in Row:
            if ele.ctype == 0:
                break
            RowData = np.insert(RowData, len(RowData), ele.value)
    else:
        for i in range(length):
            ele = Row[i]
            if ele.ctype == 0:
                data = np.nan
            else:
                data = ele.value
            RowData = np.insert(RowData, len(RowData), data)
    return RowData




I0 = 5.285
hI = 2.57  # mA
I_period = hI * 2  # mA
SortWithFlux = True
PlotT1diel = False

N = 50
EL = 0.437
EC = 2.265
EJ = 6.487
loss_tan = 3.33e-6
Q_cap = 1 / loss_tan
T = 1e-3
f01 = 0.4864
Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
File = 'wg5 in 8.5GHz cavity 0612 cd.xlsx'
ReadDataFromExcel = True
# CurrentSingle = np.array([2.5, 3.2, 3.3, 3.4, 3.5, 3.6, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.674])
# T1Single = np.array([8.81, 19.3, 11.9, 29.6, 30.2, 27.3, 25.2, 23, 28, 41.6, 42.3, 48.1, 11.5])
# T1ErrSingle = np.array([0.17, 0.49, 0.57, 0.66, 0.37, 0.5, 1.1, 1.9, 3.9, 4.9, 1.2, 1.2, 2.4])
# T2Single = np.array([3.18, 4.85, 3.42, 3.8, 3.79, 3.75, 3.27, 4.11, 3.34, 4.29, 5.66, 5.33, 7.73])
# T2ErrSingle = np.array([0.11, 0.42, 0.25, 0.26, 0.26, 0.35, 0.26, 0.48, 0.3, 0.54, 0.64, 0.27, 0.95])
#
# FreqRepeated = np.array(
#     [7.18915, 7.1005, 6.8807, 6.338, 6.0534, 5.7595, 5.462, 5.169, 4.87, 3.0522, 2.7488, 2.4463, 2.1462])
# CurrentRepeated = np.array([2.22, 2.3, 2.4, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.7, 3.8, 3.9, 4])
# T1Repeated = np.array([2.87, 3.74, 1.61, 7.37, 17.3, 21.2, 21.4, 17.1, 11.5, 8.98, 34.9, 17.1, 45.4])
# T1ErrRepeated = np.array([0.22, 0.4, 0.068, 1.1, 1.3, 4.6, 0.78, 0.82, 0.68, 1, 2.5, 3.1, 5.3])
# T2Repeated = np.array([4.34, 3.62, 2.05, 3.01, 3.77, 3.56, 3.85, 3.17, 2.97, 3.24, 3.6, 2.72, 3.46])
# T2ErrRepeated = np.array([0.32, 0.41, 0.075, 0.76, 0.072, 0.17, 0.088, 0.11, 0.075, 0.4, 0.07, 0.32, 0.077])
SpuriousMode = np.array([])
if ReadDataFromExcel:
    book = open_workbook(Folder + File, on_demand=True)
    sheet = book.sheet_by_index(0)
    FSrow = sheet.row(0)[1:]
    FreqSingle = getRowData(sheet, 0)
    Len = len(FreqSingle)
    CurrentSingle = getRowData(sheet, 1, Len)
    T1Single = getRowData(sheet, 2, Len)
    T1ErrSingle = getRowData(sheet, 3, Len)
    T2Single = getRowData(sheet, 4, Len)
    T2ErrSingle = getRowData(sheet, 5, Len)

    FreqRepeated = getRowData(sheet, 6)
    Len = len(FreqRepeated)
    CurrentRepeated = getRowData(sheet, 7, Len)
    T1Repeated = getRowData(sheet, 8, Len)
    T1ErrRepeated = getRowData(sheet, 9, Len)
    T2Repeated = getRowData(sheet, 10, Len)
    T2ErrRepeated = getRowData(sheet, 11, Len)

FluxSingle = (CurrentSingle - I0) / I_period
FluxRepeated = (CurrentRepeated - I0) / I_period

FreqMerge = np.concatenate((FreqRepeated, FreqSingle))
FluxMerge = np.concatenate((FluxRepeated, FluxSingle))
T1Merge = np.concatenate((T1Repeated, T1Single))
T2Merge = np.concatenate((T2Repeated, T2Single))
if SortWithFlux:
    SortInd = FluxMerge.argsort()
else:
    SortInd = FreqMerge.argsort()
FreqMerge = FreqMerge[SortInd]
FluxMerge = FluxMerge[SortInd]
T1Merge = T1Merge[SortInd]
T2Merge = T2Merge[SortInd]

if PlotT1diel:
    MinFlux = np.min(FluxMerge)
    MaxFlux = np.max(FluxMerge)
    FluxPlot = np.linspace(MinFlux, MaxFlux, 21)
    T1diel = FluxPlot * 0
    for i, flux in enumerate(FluxPlot):
        [pem, freq] = np.abs(ssj.phase_matrix_element_freq(N, EL, EC, EJ, flux * 2 * np.pi, 0, 1))
        # print('flux=%.3G, pem=%.3G, freq=%.3G' % (flux, pem, freq))
        T1diel[i] = 10e6 / ssj.relaxation_rate_cap(EL, EC, EJ, Q_cap, freq, pem, T)
        qsf.printPercent(i, 21)

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

if PlotT1diel:
    plt.plot(FluxPlot, T1diel)

if PlotT1diel:
    plt.legend(['Dielectric loss tangent = %.3G' % loss_tan, 'T1 - Single', 'T2echo - Single', 'T1 - Averge',
                'T2echo - Averge'])
else:
    plt.legend(['T1 - Single', 'T2echo - Single', 'T1 - Averge', 'T2echo - Averge'])
plt.plot(FluxMerge, T1Merge, ':')
plt.plot(FluxMerge, T2Merge, ':')
plt.xlabel('Flux/Phi_0', fontsize='x-large')
plt.ylabel('Decay time(us)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
if PlotT1diel:
    ax.set_yscale('log')

fig, ax = plt.subplots()
plt.plot(FluxMerge, FreqMerge, 'o:')
plt.xlabel('Flux/Phi_0', fontsize='x-large')
plt.ylabel('Freq(GHz)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
