import os
import sys
import numpy as np
from matplotlib import pyplot as plt

if __file__ in [f for f in os.listdir('.') if os.path.isfile(f)]:
    SCRIPT_PATH = os.path.dirname(os.getcwd())
else:
    SCRIPT_PATH = os.path.dirname(__file__)
LOCAL_PATH = os.path.join(SCRIPT_PATH.rsplit('fluxonium-waveguide', 1)[0],
                                             'fluxonium-waveguide')
if LOCAL_PATH not in sys.path:
    sys.path.append(LOCAL_PATH)

from local_settings import samples_path
import utilities
from plotting import colors
import anchor_points
import plot_spectrum
import plot_anchor_points
from constants import GHz2nH, GHz2fF


def main():
    sample = 'ziggy4'
    fitpath = 'Processed Data/Fluxonium/'
    pltpath = os.path.join(samples_path, sample, 'Plots/Fits/Fluxonium')
    if not os.path.exists(pltpath):
        os.makedirs(pltpath)
    filename_in = 'ziggy4_in_6.5GHz_waveguide_1.hdf5'
    
    filename_out = os.path.splitext(filename_in)[0]

    filename = os.path.join(samples_path, sample, fitpath, filename_in)
    params = utilities.load_fit(filename)

    utilities.print_params(params)

    plot_spectrum.plot_spectrum()
    
    title = ('\n$E_L/h=$%.3f GHz, $E_C/h=$%.3f GHz, $E_J/h=$%.3f GHz\n'
           % (params['E_L'], params['E_C'], params['E_J']))
    plot_spectrum.label_axes(title, title_color='w')

    fig_path = os.path.join(pltpath, '%s_spectrum.png' % filename_out)
    plt.savefig(fig_path, dpi=600)

    # plot_anchor_points.plot_anchor_points(params['data_set'])
    # fig_path = os.path.join(samples_path, sample, 'Fits/Fluxonium',
            # '%s_spectrum_anchor_points.png' % filename_out)
    # plt.savefig(fig_path, dpi=600)

    phi_ext = params['phi_ext']
    levels = params['levels']

    for idx in range(1, params['num_qbt']-1):
        freq = levels[:,idx] - levels[:,0]
        plt.plot(phi_ext, freq,
                color=colors[idx % len(colors)], lw=1.5, ls='-')

    # for idx in range(2, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,1]),
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')
    # for idx in range(2, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,1]) / 2,
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')

    # for idx in range(3, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:, idx] - levels[:, 2]),
    #              color=colors[idx % len(colors)], lw=1.0, ls='-.')
    # for idx in range(3, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,2]) / 2,
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')
    #
    # # print(params['num_qbt']-4)
    # for idx in range(4, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,3]),
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')
    # for idx in range(4, params['num_qbt']):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,3]) / 2,
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')

    # for idx in range(1, params['num_qbt']-1):
    #     plt.plot(phi_ext, (levels[:,idx] - levels[:,0]) / 2.,
    #         color=colors[idx % len(colors)], lw=1.0, ls='-.')

    plot_spectrum.label_axes(title, title_color='k')
    plt.tight_layout()
    fig_path = os.path.join(pltpath, '%s_spectrum_fit.png' % filename_out)
    plt.savefig(fig_path, dpi=600)

    plt.show()


if __name__ == '__main__':
    main()
