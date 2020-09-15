import os
import sys
import numpy as np
from matplotlib import pyplot as plt

import Labber

if __file__ in [f for f in os.listdir('.') if os.path.isfile(f)]:
    SCRIPT_PATH = os.path.dirname(os.getcwd())
else:
    SCRIPT_PATH = os.path.dirname(__file__)
LOCAL_PATH = os.path.join(SCRIPT_PATH.rsplit('fluxonium-waveguide', 1)[0],
                          'fluxonium-waveguide')
if LOCAL_PATH not in sys.path:
    sys.path.append(LOCAL_PATH)

from plotting import correct_grid
from local_settings import samples_path

sample = 'Augustus 17'


def plot_spectrum():
    # I0 = -0.36  # mA
    # I0 = 175e-3  # mA
    # hI = 502e-3 - I0  # mA
    hI = (1.501 - 0.6675)  # mA
    I0 = 0.6675  # mA
    path_filenames = [
        # '023bis.hdf5',
        '024.hdf5',
        # '024c.hdf5',
        # '024d.hdf5',
        # '062e.hdf5',
        'Two_tone_check_0.hdf5',
        'Two_tone_check.hdf5',
        # 'Two_tone_check_2.hdf5',
        # 'Two_tone_check_3.hdf5',
        'Two_tone_check_4.hdf5',
    ]

    path_filenames = [os.path.join(samples_path, sample, pf)
                      for pf in path_filenames]
    for path_filename in path_filenames:
        # data_name = 'Alazar - Channel A - Average demodulated value'
        data_name = 'Signal Demodulation - Value'
        f = Labber.LogFile(path_filename)
        # print(path_filename.split('\\')[-1])
        # freq_var = 'Drive2 - Frequency'
        freq_var = 'Qubit1 - Frequency'
        bias_var = 'Yoko2FluxoniumExperiment - Current'
        freq_2d = f.getData(freq_var) / 1.e9
        data_2d = data = f.getData(f.getLogChannels()[0]['name'])
        bias_2d = f.getData(bias_var) / 1.e-3

        n_entries = f.getNumberOfEntries()

        bias_size = bias_2d.shape[0]
        freq_size = bias_2d.shape[1]

        if n_entries != bias_size:
            bias_size = n_entries
            for k in range(bias_size):
                d = f.getEntry(k)
                # print(d)
                freq_1d = d[freq_var] / 1.e9
                bias_1d = d[bias_var] / 1.e-3

                freq_size = freq_1d.size

                if k == 0:
                    data_2d = np.empty((bias_size, freq_size),
                                       dtype=complex)
                    freq_2d = np.empty_like(data_2d, dtype=float)
                    bias_2d = np.empty_like(data_2d, dtype=float)
                    plot_2d = np.empty_like(data_2d, dtype=float)

                data_2d[k] = d[data_name]
                freq_2d[k] = freq_1d
                bias_2d[k] = bias_1d
        else:
            bias_size = bias_2d.shape[0]
            freq_size = bias_2d.shape[1]

            freq_2d.shape = (bias_size, freq_size)
            bias_2d.shape = (bias_size, freq_size)
            data_2d.shape = (bias_size, freq_size)
            plot_2d = np.zeros_like(data_2d, dtype=np.float)

        bias_2d -= I0  # mA
        bias_2d /= (2. * hI)  # mA

        if path_filename.endswith('024.hdf5'):
            bias_2d += 0.5 - 0.4964
        # print(abs(np.mean(data_2d, axis=1)))
        # print(abs(np.var(data_2d, axis=1)))
        plot_2d = np.transpose((data_2d.transpose() / np.nanmean(data_2d, axis=1)))
        # plot_2d = data_2d
        # plot_2d = plot_2d / np.mean(plot_2d)
        plot_2d_max_contrast = np.zeros(plot_2d.shape)
        if path_filename.endswith('024.hdf5'):
            plot_2d_max_contrast = np.imag(plot_2d)
        else:
            for k in range(plot_2d.shape[0]):
                if np.var(np.imag(plot_2d[k, :])) < np.var(np.real(plot_2d[k, :])):
                    plot_2d_max_contrast[k, :] = np.real(plot_2d[k, :])
                else:
                    plot_2d_max_contrast[k, :] = np.imag(plot_2d[k, :])
        plot_2d = plot_2d_max_contrast
        plot_2d = np.transpose((plot_2d.transpose() - np.mean(plot_2d, axis=1)))
        # plot_2d = np.transpose((plot_2d.transpose() / np.var(plot_2d, axis=1)))

        # plot_2d = np.abs(data_2d - np.mean(data_2d))

        # for k in range(data_2d.shape[1]):
        #     plot_2d[:, k] /= np.max(plot_2d[:, k])
        # plot_2d = 1. - np.abs(plot_2d)
        plot_2d -= np.nanmean(plot_2d)
        plot_2d = np.abs(plot_2d)
        # print(np.max(plot_2d))


        # plot_2d_average_flux = np.average(plot_2d, axis=0)
        # bad_ind = plot_2d_average_flux > 0.05
        # plot_2d[:, bad_ind] *= 1e-2
        # # print(plot_2d_average_flux[bad_ind])
        # plot_2d /= (np.nanmax(plot_2d))


        # print(np.max(plot_2d))
        # print(np.min(plot_2d))
        # print(np.mean(plot_2d))
        # print(np.max(plot_2d))
        bias_2d, freq_2d = correct_grid(bias_2d, freq_2d)
        # colormap = 'PuBu'
        colormap = 'jet'
        # if path_filename.endswith('024.hdf5'):
        if 0:
            vmax = 0.2
            vmin = 0
            print(np.nanmax(plot_2d))
            print(np.nanmin(plot_2d))
        else:
            plot_2d /= (np.nanmax(plot_2d))
            plot_2d = np.log(plot_2d)
            vmax = np.nanmax(plot_2d)
            vmin = np.nanmin(plot_2d)
        plt.pcolormesh(bias_2d, freq_2d, plot_2d, cmap=colormap, vmin=vmin, vmax=vmax)
        # plt.pcolormesh(bias_2d, freq_2d, plot_2d)

        # print(np.min(plot_2d))
        # print(np.max(plot_2d))


def label_axes(title='', xlim=[0., 1.], ylim=[0, 7],
               title_color='k'):
    labelsize = 18
    plt.xlabel('$\Phi_\mathrm{ext}/\Phi_0$', fontsize=labelsize)
    plt.ylabel('Frequency (GHz)', fontsize=labelsize)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title, fontsize=labelsize - 3, color=title_color)
    plt.gcf().set_size_inches(7.5, 10, forward=True)
    plt.gca().tick_params(which='major', direction='in', right=True,
                          top=True, labelsize=labelsize, length=6, width=1)
    plt.gca().tick_params(which='minor', direction='in', right=True,
                          top=True, length=3, width=1)
    for axis in ['left', 'right', 'top', 'bottom']:
        plt.gca().spines[axis].set_linewidth(1)
    plt.tight_layout()


if __name__ == '__main__':
    FigPath = os.path.join(samples_path, sample, 'Plots')
    if not os.path.exists(FigPath):
        os.makedirs(FigPath)
    filename = os.path.join(FigPath, 'spectrum.png')
    plot_spectrum()
    label_axes()


    def onclick(event):
        print('{\'transition\': (0, 2),  # 00-10\n\'external flux quanta\': %.5G,\n\'frequency\': %.5G,\n},' % (
            event.xdata, event.ydata))


    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.savefig(filename, dpi=600)
    plt.show()
