from QubitDataProcessPackages import*

FolderName = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619\\Shared log\\'
FileName = 'magnet_20190813.txt'
DataTable = np.genfromtxt(FolderName + FileName, skip_header=23, delimiter=',')
lines = []
File = open(FolderName + FileName)
for i in range(22):
    lines += [File.readline()]
print(DataTable[0,:])