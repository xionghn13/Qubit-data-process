import numpy as np

Folder = 'E:\\Projects\\Fluxonium\\data_process\\Fluxonium032619/'
file = open(Folder + 'filenames0626')
names = file.readlines()
for n in names:
    if n.startswith('Current'):
        print('[')
    elif n.startswith('t1'):
        print('\t\'' + n[:-1] + '\',')
    elif n.startswith('t2'):
        print('\t\'' + n[:-1] + '\'')
    elif n.startswith('two'):
        print(']')
print(']')