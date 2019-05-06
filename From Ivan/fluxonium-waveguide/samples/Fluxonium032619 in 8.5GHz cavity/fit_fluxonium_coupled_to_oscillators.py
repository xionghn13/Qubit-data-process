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
    qubit.E_L = np.max([x[0], .05]) 
    qubit.E_C = np.max([x[1], .05])
    qubit.E_J = np.max([x[2], .05])

    n_modes = len(params['frequencies'])
    # params['frequencies'] = x[3:3+n_modes]
    # params['phi_couplings'][0] = x[3+n_modes]
    params['n_couplings'] = np.abs(x[3:])
    # start_time = time.time()
    oscillators = UncoupledOscillators(params)
    # print('UncoupledOscillators: %f s' % (time.time() - start_time))

    data = anchor_points.data[params['data_set']]
    flux = np.unique([point['external flux quanta'] for point in data])
    
    systems = parallel_map(diagonalization, flux,
            task_args=(qubit, oscillators, params))
    
    error_type = params['error_type']
    error = 0.
    N = len(data)
    # print('Oscillator cut-off energy: %f GHz.' % oscillators.levels()[-1])
    # print('Qubit cut-off energy: %f GHz.' % qubit.levels()[-1])
    for k in range(N):
        system = systems[np.where(flux == data[k]['external flux quanta'])[0][0]]
        freq_sim = system.frequency(*data[k]['transition'])
        freq_exp = data[k]['frequency']
        # print('%s error: %f GHz' % (data[k]['transition'],
                                    # freq_sim - freq_exp))
        if error_type == 'absolute_error':
            error += (freq_sim - freq_exp)**2.
        else:
            error += ((freq_exp - freq_sim) / freq_exp)**2.
    error = np.sqrt(error / float(N))
    if error_type == 'absolute_error':
        sys.stdout.write('Absolute error: %.3f MHz.    \r' % (1.e3 * error))
    else:
        sys.stdout.write('Relative error: %.5f%%.      \r' % (1.e2 * error))
    # print('\nTotal time: %f s' % (time.time() - total_time))

    return error
    
def fit(params):
    start_time = time.time()
    qubit = Fluxonium(params)
    res = scipy.optimize.minimize(fun=error_function,
            method='Nelder-Mead',
            x0=[params['E_L'], params['E_C'], params['E_J'],
                *params['n_couplings']],
            args=(qubit, params),)
    run_time = time.time() - start_time
    print(res)
    
    params['E_L'] = res.x[0]
    params['E_C'] = res.x[1]
    params['E_J'] = res.x[2]
    # n_modes = len(params['frequencies'])
    # params['frequencies'] = res.x[3:3+n_modes]
    # params['phi_couplings'][0] = res.x[3+n_modes]
    params['n_couplings'] = np.abs(res.x[3:])

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
    sample = 'fluxonium022319'
    subpath = 'Processed Data/Fluxonium Coupled to Harmonic Modes/'
    path = os.path.join(samples_path, sample, subpath)

    filename_in = 'one_mode_2.hdf5'
    filename_fit = 'one_mode_3.hdf5'

    params = utilities.load_fit(os.path.join(path, filename_in))

    # Initial guess for the qubit parameters.
    # params = {'E_L': 0.200, # The inductive energy.
              # 'E_C': 6.842, # The charging energy.
              # 'E_J': 4.575, # The Josephson energy.
              # 'phi_ext': 0.,
              # 'frequencies': [5.979, 8.548, 9.823], # Single mode frequencies.
              # 'n_couplings': [0., 5.397, 2.667], # Charge couplings.
              # 'phi_couplings': [0.153, 0., 0.], # Flux couplings.
              # 'num_mod': [10, 7, 7],
              # 'num_osc': 100,
              # 'num_qbt': 15,
              # 'num_tot': 21,
              # 'num_cpl': 50,
              # 'error_type': 'relative_error',
              # 'data_set': 'data4'
             # }

    params['num_qbt'] = 10
    params['num_tot'] = 10
    params['num_cpl'] = 10
    params['num_mod'] = [10]
    params['frequencies'] = np.array([7.])
    params['n_couplings'] = np.array([0.])
    params['phi_couplings'] = np.array([0.])
    params['data_set'] = 'data2'
    
    phi_ext = np.linspace(-1., 1., 11)

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
    for idx in range(params['num_tot']-1):
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
