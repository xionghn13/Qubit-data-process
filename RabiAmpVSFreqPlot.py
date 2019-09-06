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


I0 = -82.5e-3  # mA
hI = 2.571  # mA
I_period = hI * 2  # mA
SortWithFlux = False
PlotT1diel = False


Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
File = 'wg5 in 8.5GHz cavity rabi 0710 cd.xlsx'
ReadDataFromExcel = True

SpuriousMode = [0.51373, 0.51491, 0.51508, 0.51596, 0.51883, 0.52205, 0.52262]
PlotSpuriousMode = False

if ReadDataFromExcel:
    book = open_workbook(Folder + File, on_demand=True)
    sheet = book.sheet_by_index(0)
    FSrow = sheet.row(0)[1:]
    Freq = getRowData(sheet, 0)
    Len = len(Freq)

    Current = getRowData(sheet, 1, Len)
    A = getRowData(sheet, 3, Len)
    AErr = getRowData(sheet, 4, Len)


Flux = (Current - I0) / I_period

FluxSort_Ind = Flux.argsort()
FreqSort_Ind = Freq.argsort()
FreqSort_Freq = Freq[FreqSort_Ind]
FluxSort_Freq = Freq[FluxSort_Ind]
FluxSort_Flux = Flux[FluxSort_Ind]
FreqSort_A = A[FreqSort_Ind]
FluxSort_A = A[FluxSort_Ind]


fig, ax = plt.subplots()
ax.errorbar(Freq, A, yerr=AErr, fmt='bo')
plt.plot(FreqSort_Freq, FreqSort_A, ':')
if PlotSpuriousMode:
    for mode in SpuriousMode:
        plt.axvline(x=mode, color='g', linestyle='--')
plt.xlabel('Freq(GHz)', fontsize='x-large')
plt.ylabel('Rabi Amplitude', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
leg = []
ax.errorbar(Flux, A, yerr=AErr, fmt='bo')
plt.plot(FluxSort_Flux, FluxSort_A, ':')
plt.xlabel('Flux/Phi_0', fontsize='x-large')
plt.ylabel('Rabi Amplitude', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
plt.plot(FluxSort_Flux, FluxSort_Freq, 'o:')
plt.xlabel('Flux/Phi_0', fontsize='x-large')
plt.ylabel('Freq(GHz)', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

plt.show()
