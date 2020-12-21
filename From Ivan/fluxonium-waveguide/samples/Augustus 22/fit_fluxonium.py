import os
import sys
import time
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize

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
from model_hamiltonians import TwoCoupledFluxoniums
from plotting import colors

import anchor_points
import plot_spectrum
import plot_anchor_points


def error_function(x, qubit, params):
    qubit.E_L = np.max([0.001, x[0]])
    qubit.E_C = np.max([0.001, x[1]])
    qubit.E_J = np.max([0.001, x[2]])

    data = anchor_points.data[params['data_set']]
    error_type = params['error_type']
    error = 0.
    N = len(data)
    phi_ext_prev = np.nan
    for k in range(N):
        phi_ext = data[k]['external flux quanta']
        if phi_ext != phi_ext_prev:
            qubit.phi_ext = phi_ext
            phi_ext_prev = phi_ext
        freq_sim = qubit.frequency(*data[k]['transition'])
        freq_exp = data[k]['frequency']
        if error_type == 'absolute_error':
            error += (freq_sim - freq_exp)**2.
        else:
            error += ((freq_exp - freq_sim) / freq_exp)**2.
    error = np.sqrt(error / float(N))
    if error_type == 'absolute_error':
        sys.stdout.write('Absolute error: %.3f MHz.    \r' % (1.e3 * error))
    else:
        sys.stdout.write('Relative error: %.5f%%.      \r' % (1.e2 * error))
    return error

def fit(qubit, params):
    start_time = time.time()
    res = scipy.optimize.minimize(fun=error_function,
            method='Nelder-Mead',
            x0=[params['E_L'], params['E_C'], params['E_J']],
            args=(qubit, params),)
    print(res)
    run_time = time.time() - start_time
    print('Optimization run time: %f sec.' % run_time)
    params['run_time'] = run_time

    params['E_L'] = res.x[0]
    params['E_C'] = res.x[1]
    params['E_J'] = res.x[2]

    if params['error_type'] == 'absolute_error':
        params['absolute_error'] = res.fun
        if 'relative_error' in params:
            del params['relative_error']
    else:
        params['relative_error'] = res.fun
        if 'absolute_error' in params:
            del params['absolute_error']

    return params

def main():
    sample = 'ziggy4'
    subpath = 'Processed Data/Fluxonium/'
    path = os.path.join(samples_path, sample, subpath)

    filename_in = 'ziggy4_in_6.5GHz_waveguide_1.hdf5'
    filename_fit = 'ziggy4_in_6.5GHz_waveguide_0628.hdf5'
    
    params = utilities.load_fit(os.path.join(path, filename_in))

    # Initial guess for the qubit parameters.
    # params = {'E_L': 0.617, # The inductive energy.
    #           'E_C': 1.172, # The charging energy.
    #           'E_J': 2.024, # The Josephson energy.
    #           'num_osc': 100,
    #           'num_qbt': 10,
    #           'phi_ext': [0.]
    #          }
             
    # params['num_osc'] = 50
    # params['error_type'] = 'absolute_error'
    params['data_set'] = 'data2'
    
    # params['E_L'] = .1
    phi_ext1 = np.linspace(0.495, 0.535, 101)
    phi_ext = np.linspace(0.4, 1.35, 96)
    phi_ext = np.sort(np.concatenate((phi_ext1, phi_ext)))

    utilities.print_params(params)
    qubit = Fluxonium(params)

    # params = fit(qubit, params)
    utilities.print_params(params)
    utilities.save_fit(os.path.join(path, filename_fit), params)

    plot_spectrum.plot_spectrum()
    plot_anchor_points.plot_anchor_points(params['data_set'])

    # Compute energy levels.
    levels = np.zeros((phi_ext.size, params['num_qbt']))
    freqs = np.zeros((phi_ext.size, params['num_qbt']))

    system = qubit
    for idx_phi_ext, value in enumerate(phi_ext):
        system.phi_ext = value
        levels[idx_phi_ext] = system.levels()
        freqs[idx_phi_ext] = system.levels() - system.levels()[0]
        sys.stdout.write('Progress: %5.1f%%\r'
                % (100. * (idx_phi_ext + 1) / len(phi_ext)))
    
    params['levels'] = levels
    params['phi_ext'] = phi_ext

    utilities.save_fit(os.path.join(path, filename_fit), params)

    for idx in range(1, params['num_qbt']-1):
        plt.plot(phi_ext, levels[:,idx] - levels[:,0],
                color=colors[idx % len(colors)])

    title = ('$E_L/h=$%.3f GHz, $E_C/h=$%.3f GHz, $E_J/h=$%.3f GHz'
           % (params['E_L'], params['E_C'], params['E_J']))
    plot_spectrum.label_axes(title, ylim=[0, 10])

    plt.show()


if __name__ == '__main__':
    main()
