from QubitDataProcessPackages import *
import Single_small_junction as ssj
from xlrd import open_workbook
from CoherenceVSFreqPlot import sortPlot

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


# I0 = -82.5e-3  # mA
# hI = 2.571  # mA
# I_period = hI * 2  # mA
# SortWithFlux = False
# PlotT1diel = False


Folder = 'E:\Projects\Fluxonium\data_process\Fluxonium032619/'
File = 'wg5 in 8.5GHz cavity rabi CH1 pumped 0830 cd.xlsx'
ReadDataFromExcel = True


if ReadDataFromExcel:
    book = open_workbook(Folder + File, on_demand=True)
    sheet = book.sheet_by_index(0)
    FSrow = sheet.row(0)[1:]
    Freq = getRowData(sheet, 0)
    Len = len(Freq)

    Current = getRowData(sheet, 1, Len)
    Delay = getRowData(sheet, 2, Len) * 1e6
    A = getRowData(sheet, 3, Len)
    AErr = getRowData(sheet, 4, Len)



fig, ax = plt.subplots()
ax.grid(linestyle='--')
ax.errorbar(Delay, A, yerr=AErr, fmt='bo')
sortPlot(Delay, A, ':')
plt.xlabel('Delay(us)', fontsize='x-large')
plt.ylabel('Rabi Amplitude', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()
