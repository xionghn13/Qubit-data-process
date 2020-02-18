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


    sample = 'fluxonium032619'
    fitpath = 'Processed Data/Fluxonium Coupled to Harmonic Modes/'
    pltpath = os.path.join(samples_path, sample,
            'Plots/Fits/Fluxonium Coupled to Harmonic Modes')

    filename_in = 'two_mode_in_8.5GHz_cavity_0830cd.hdf5'
    if not os.path.exists(pltpath):
        os.makedirs(pltpath)


    filename_out = os.path.splitext(filename_in)[0]

    filename = os.path.join(samples_path, sample, fitpath, filename_in)
    params = utilities.load_fit(filename)

    utilities.print_params(params)

    # plot_spectrum.plot_spectrum()
    # plot_anchor_points.plot_anchor_points(params['data_set'])
    # plot_anchor_points.plot_anchor_points('data2')

    str_freq = ', '.join(['%.3f' % freq for freq in params['frequencies']])
    str_coup = ', '.join(['%.3f|%.3f' % (c, f) for c, f in
            zip(params['n_couplings'], params['phi_couplings'])])

    title = ('$E_L/h=$%.3f GHz, $E_C/h=$%.3f GHz, $E_J/h=$%.3f GHz\n'
            '$f_i\in${%s} GHz\n$g_{i,n|\phi}\in${%s} GHz' % (params['E_L'],
            params['E_C'], params['E_J'], str_freq, str_coup))

    # plt.xlim(-.05, 0.525)
    plot_spectrum.label_axes(title, title_color='w', xlim=[-0.5, 0.6])
    fig_path = os.path.join(pltpath, '%s_spectrum.png' % filename_out)
    plt.savefig(fig_path, dpi=600)

    phi_ext = params['phi_ext']
    levels = params['levels']
    weights = params['weights']
    num_tot = params['num_tot']

    lines = colored_lines(phi_ext, levels, weights, 0, num_tot)
    for idx in range(params['num_tot']-1):
        plt.gca().add_collection(lines[idx])

    # lines = colored_lines(phi_ext, levels, weights, 1, num_tot)
    # for idx in range(params['num_tot']-2):
    #     plt.gca().add_collection(lines[idx], '-')

    # lines = colored_lines(phi_ext, levels, weights, 2, num_tot)
    # for idx in range(params['num_tot']-3):
    #     plt.gca().add_collection(lines[idx], '-')
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
