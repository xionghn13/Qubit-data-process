import ExtractDataFunc as edf
import matplotlib.pyplot as plt
import numpy as np

test_path = 'E:/Projects\Fluxonium\data_process/test_data/'

ATS_var = 'Alazar - Channel A - Average demodulated value'
# test_out = os.path.join(test_path, 'test digitizer_37.hdf5')
test_out = test_path + 'test digitizer_37.hdf5'
[BackFreq, BackComplex] = edf.readFSweepLabber(test_out)

fig, ax = plt.subplots()
plt.plot(BackFreq, np.abs(BackComplex))
plt.xlabel('freq/GHz', fontsize='x-large')
plt.ylabel('Abs', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.show()