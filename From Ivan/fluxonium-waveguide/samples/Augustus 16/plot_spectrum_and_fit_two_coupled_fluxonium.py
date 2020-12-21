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
from plotting import colored_lines
import plot_spectrum
import plot_anchor_points


def main():
    sample = 'Augustus 16'
    fitpath = 'Processed Data/Two Coupled Fluxoniums/'
    pltpath = os.path.join(samples_path, sample,
                           'Plots/Fits/Two Coupled Fluxoniums')
    if not os.path.exists(pltpath):
        os.makedirs(pltpath)

    filename_in = 'Augustus_16.hdf5'
    filename_out = os.path.splitext(filename_in)[0]

    filename = os.path.join(samples_path, sample, fitpath, filename_in)
    params = utilities.load_fit(filename)

    utilities.print_params(params)

    plot_spectrum.plot_spectrum()
    plot_anchor_points.plot_anchor_points(params['data_set'])
    # plot_anchor_points.plot_anchor_points('data2')

    title = ('$E_L1/h=$%.3f GHz, $E_C1/h=$%.3f GHz, $E_J1/h=$%.3f GHz\n'
             '$E_L2/h=$%.3f GHz, $E_C2/h=$%.3f GHz, $E_J2/h=$%.3f GHz\n'
             '$J_C/h=$%.3f GHz' % (
                 params['E_L1'], params['E_C1'], params['E_J1'], params['E_L2'], params['E_C2'], params['E_J2'],
                 params['E_int_chg']))


    # plt.xlim(-.05, 0.525)
    plot_spectrum.label_axes(title, title_color='w')
    fig_path = os.path.join(pltpath, '%s_spectrum.png' % filename_out)
    plt.savefig(fig_path, dpi=600)

    phi_ext = params['phi_ext']
    levels = params['levels']
    weights = params['weights']
    num_tot = params['num_tot']


    number_photons = 1
    lines = colored_lines(phi_ext, levels, weights, 0, num_tot, color2QB=True, num_photons=number_photons)
    for idx in range(params['num_tot'] - 1):
        plt.gca().add_collection(lines[idx])

    lines = colored_lines(phi_ext, levels, weights, 1, num_tot, color2QB=True, num_photons=number_photons)
    for idx in range(params['num_tot'] - 2):
        plt.gca().add_collection(lines[idx])

    lines = colored_lines(phi_ext, levels, weights, 2, num_tot, color2QB=True, num_photons=number_photons)
    for idx in range(params['num_tot'] - 3):
        plt.gca().add_collection(lines[idx])

    lines = colored_lines(phi_ext, levels, weights, 3, num_tot, color2QB=True, num_photons=number_photons)
    for idx in range(params['num_tot'] - 4):
        plt.gca().add_collection(lines[idx])

    # lines = colored_lines(phi_ext, levels, weights, 7, num_tot, color2QB=True)
    # for idx in range(params['num_tot'] - 8):
    #     plt.gca().add_collection(lines[idx])
    #
    # lines = colored_lines(phi_ext, levels, weights, 3, num_tot)
    # for idx in range(params['num_tot']-4):
    #     plt.gca().add_collection(lines[idx], '-')
    # lines = colored_lines(phi_ext, levels, weights, 0, num_tot, 2)
    # for idx in range(params['num_tot']-1):
    # plt.gca().add_collection(lines[idx])

    plot_spectrum.label_axes(title, title_color='k')
    plt.tight_layout()
    fig_path = os.path.join(pltpath, '%s_spectrum_fit.png' % filename_out)
    plt.savefig(fig_path, dpi=600)

    plt.show()

if __name__ == '__main__':
    main()
