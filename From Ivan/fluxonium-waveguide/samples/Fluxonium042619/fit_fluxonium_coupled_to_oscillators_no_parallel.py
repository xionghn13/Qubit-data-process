import os
import sys
import time
import numpy as np
import scipy.optimize
from qutip.parallel import parallel_map
# from qutip.parallel import serial_map as parallel_map
from matplotlib import pyplot as plt

if __file__ in [f for f in os.listdir('.') if os.path.isfile(f)]:
    SCRIPT_PATH = os.path.dirname(os.getcwd())
else:
    SCRIPT_PATH = os.path.dirname(__file__)
LOCAL_PATH = os.path.join(SCRIPT_PATH.rsplit('fluxonium-waveguide', 1)[0],
                          'fluxonium-waveguide')
if LOCAL_PATH not in sys.path:
    sys.path.append(LOCAL_PATH)

import utilities
from local_settings import samples_path
from simple_hamiltonians import Fluxonium
from model_hamiltonians import (UncoupledOscillators,
                                FluxoniumCoupledToOscillators)
from plotting import colored_lines

import anchor_points
import plot_spectrum
import plot_anchor_points


def diagonalization(phi_ext, qubit, oscillators, params):
    # start_time = time.time()
    qubit.phi_ext = phi_ext
    system = FluxoniumCoupledToOscillators(qubit, oscillators, params)
    system.levels()
    # print('FluxoniumCoupledToOscillators: %f s' % (time.time() - start_time))
    return system


def error_function(x, qubit, params):
    # total_time = time.time()
    # qubit.E_L = np.max([x[0], .05])
    # qubit.E_C = np.max([x[1], .05])
    # qubit.E_J = np.max([x[2], .05])

    n_modes = len(params['frequencies'])
    params['frequencies'] = x[3:3 + n_modes]
    # params['phi_couplings'][0] = x[3+n_modes]
    params['n_couplings'] = np.abs(x[3 + n_modes:])
    # start_time = time.time()
    oscillators = UncoupledOscillators(params)
    # print('UncoupledOscillators: %f s' % (time.time() - start_time))

    data = anchor_points.data[params['data_set']]
    error_type = params['error_type']
    error = 0.
    N = len(data)
    # print('Oscillator cut-off energy: %f GHz.' % np.real(oscillators.H()[-1,-1]))
    # print('Qubit cut-off energy: %f GHz.' % np.real(qubit.H()[-1,-1]))
    prev_phi_ext = np.nan
    for k in range(N):
        # start_time = time.time()
        phi_ext = data[k]['external flux quanta']
        if prev_phi_ext != phi_ext:
            qubit.phi_ext = phi_ext
            system = FluxoniumCoupledToOscillators(qubit,
                                                   oscillators, params)
            prev_phi_ext = phi_ext
        # run_time = time.time() - start_time
        # print('%s %f' % (str(data[k]['transition']), run_time))
        freq_sim = system.frequency(*data[k]['transition'])
        freq_exp = data[k]['frequency']
        # print('%s error: %f GHz' % (data[k]['transition'],
        # freq_sim - freq_exp))
        if error_type == 'absolute_error':
            error += (freq_sim - freq_exp) ** 2.
        else:
            error += ((freq_exp - freq_sim) / freq_exp) ** 2.
    error = np.sqrt(error / float(N))
    if error_type == 'absolute_error':
        sys.stdout.write('Absolute error: %.3f MHz.    \r' % (1.e3 * error))
    else:
        sys.stdout.write('Relative error: %.5f%%.      \r' % (1.e2 * error))
    return error


def fit(params):
    start_time = time.time()
    qubit = Fluxonium(params)
    res = scipy.optimize.minimize(fun=error_function,
                                  method='Nelder-Mead',
                                  x0=[params['E_L'], params['E_C'], params['E_J'], *params['frequencies'],
                                      *params['n_couplings']],
                                  args=(qubit, params), )
    run_time = time.time() - start_time
    print(res)

    params['E_L'] = res.x[0]
    params['E_C'] = res.x[1]
    params['E_J'] = res.x[2]
    n_modes = len(params['frequencies'])
    params['frequencies'] = res.x[3:3 + n_modes]
    # params['phi_couplings'][0] = res.x[3+n_modes]
    params['n_couplings'] = np.abs(res.x[3 + n_modes:])

    if params['error_type'] == 'absolute_error':
        params['absolute_error'] = res.fun
        if 'relative_error' in params:
            del params['relative_error']
    else:
        params['relative_error'] = res.fun
        if 'absolute_error' in params:
            del params['absolute_error']

    print('Optimization run time: %f sec.' % run_time)
    params['run_time'] = run_time

    return params


def main():
    sample = 'fluxonium042619'
    subpath = 'Processed Data/Fluxonium Coupled to Harmonic Modes/'
    path = os.path.join(samples_path, sample, subpath)

    filename_in = 'one_mode_anticrossing.hdf5'
    filename_fit = 'one_mode_anticrossing_2.hdf5'

    params = utilities.load_fit(os.path.join(path, filename_in))

    # Initial guess for the qubit parameters.
    # params = {'E_L': 0.66,  # The inductive energy.
    #           'E_C': 3,  # The charging energy.
    #           'E_J': 9.4,  # The Josephson energy.
    #           'phi_ext': [0.],
    #           'frequencies': [7.45],  # Single mode frequencies.
    #           'n_couplings': [0.],  # Charge couplings.
    #           'phi_couplings': [0],  # Flux couplings.
    #           'num_mod': [10],
    #           'num_osc': 100,
    #           'num_qbt': 10,
    #           'num_tot': 10,
    #           'num_cpl': 10,
    #           'error_type': 'relative_error',
    #           'data_set': 'data0'
    #           }

    # params['E_L'] = 1.016
    # params['E_C'] = 1.421
    # params['E_J'] = 14.858
    # params['num_qbt'] = 10
    # params['num_tot'] = 10
    # params['num_cpl'] = 10
    params['frequencies'] = np.array([7.473])
    # params['num_mod'] = [10]
    params['n_couplings'] = np.array([0.4])
    # params['phi_couplings'] = np.array([0.])
    params['data_set'] = 'data5'
    # params['error_type'] = 'absolute_error'

    phi_ext1 = np.linspace(0.297, 0.2985, 51)
    phi_ext2 = np.linspace(0., 0.7, 51)
    phi_ext = np.sort(np.concatenate((phi_ext1, phi_ext2)))
    # phi_ext = np.linspace(0.297, 0.2985, 51)

    utilities.print_params(params)
    params = fit(params)
    utilities.print_params(params)
    utilities.save_fit(os.path.join(path, filename_fit), params)

    # Compute energy levels.
    params['num_tot'] = 50
    levels = np.zeros((phi_ext.size, params['num_tot']))
    weights = np.zeros((phi_ext.size, params['num_tot'], params['num_qbt']))

    qubit = Fluxonium(params)
    oscillators = UncoupledOscillators(params)
    for idx_phi_ext, value in enumerate(phi_ext):
        qubit.phi_ext = value
        system = FluxoniumCoupledToOscillators(qubit, oscillators, params)
        levels[idx_phi_ext] = system.levels()
        weights[idx_phi_ext] = system.weights()
        sys.stdout.write('Progress: %5.1f%%\r'
                         % (100. * (idx_phi_ext + 1) / len(phi_ext)))

    params['levels'] = levels
    params['weights'] = weights
    params['phi_ext'] = phi_ext

    utilities.save_fit(os.path.join(path, filename_fit), params)

    plot_spectrum.plot_spectrum()
    plot_anchor_points.plot_anchor_points(params['data_set'])

    lines = colored_lines(phi_ext, levels, weights, 0, params['num_tot'])
    for idx in range(params['num_tot'] - 1):
        plt.gca().add_collection(lines[idx])

    str_freq = ', '.join(['%.3f' % freq for freq in params['frequencies']])
    str_coup = ', '.join(['%.3f|%.3f' % (c, f) for c, f in
                          zip(params['n_couplings'], params['phi_couplings'])])

    title = ('$E_L/h=$%.3f GHz, $E_C/h=$%.3f GHz, $E_J/h=$%.3f GHz\n'
             '$f_i\in${%s} GHz\n$g_{i,n|\phi}\in${%s} GHz' % (params['E_L'],
                                                              params['E_C'], params['E_J'], str_freq, str_coup))
    plot_spectrum.label_axes(title)

    plt.show()


if __name__ == '__main__':
    main()
