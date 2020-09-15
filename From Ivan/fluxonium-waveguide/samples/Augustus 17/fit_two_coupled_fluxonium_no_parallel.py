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
                                FluxoniumCoupledToOscillators,
                                TwoCoupledFluxoniums, )
from plotting import colored_lines

import anchor_points
import plot_spectrum
import plot_anchor_points


# def diagonalization(phi_ext, qubit1, qubit2, params):
#     # start_time = time.time()
#     qubit1.phi_ext = phi_ext
#     dphi_ext = 0
#     ratio_phi = 1
#     qubit2.phi_ext = ratio_phi * phi_ext + dphi_ext
#     system = TwoCoupledFluxoniums(qubit1, qubit2, params)
#     system.levels()
#     # print('FluxoniumCoupledToOscillators: %f s' % (time.time() - start_time))
#     return system

def error_function(x, qubit1, qubit2, params):
    # total_time = time.time()
    qubit1.E_L = np.max([x[0], .05])
    qubit1.E_C = np.max([x[1], .05])
    qubit1.E_J = np.max([x[2], .05])

    qubit2.E_L = np.max([x[3], .05])
    qubit2.E_C = np.max([x[4], .05])
    qubit2.E_J = np.max([x[5], .05])

    params['E_int_chg'] = np.abs(x[6])

    data = anchor_points.data[params['data_set']]
    error_type = params['error_type']
    dphi_ext = params['dphi_ext']
    ratio_phi = params['ratio_phi']
    error = 0.
    N = len(data)
    # print('Oscillator cut-off energy: %f GHz.' % np.real(oscillators.H()[-1,-1]))
    # print('Qubit cut-off energy: %f GHz.' % np.real(qubit.H()[-1,-1]))
    prev_phi_ext = np.nan
    for k in range(N):
        # start_time = time.time()
        phi_ext = data[k]['external flux quanta']
        if prev_phi_ext != phi_ext:
            qubit1.phi_ext = ratio_phi * phi_ext + dphi_ext
            qubit2.phi_ext = phi_ext
            system = TwoCoupledFluxoniums(qubit1, qubit2, params)
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
    params1 = {}
    params1['E_L'] = params['E_L1']
    params1['E_C'] = params['E_C1']
    params1['E_J'] = params['E_J1']
    params1['phi_ext'] = params['phi_ext1']
    params1['num_osc'] = params['num_osc1']
    params1['num_qbt'] = params['num_qbt1']

    params2 = {}
    params2['E_L'] = params['E_L2']
    params2['E_C'] = params['E_C2']
    params2['E_J'] = params['E_J2']
    params2['phi_ext'] = params['phi_ext2']
    params2['num_osc'] = params['num_osc2']
    params2['num_qbt'] = params['num_qbt2']

    qubit1 = Fluxonium(params1)
    qubit2 = Fluxonium(params2)

    res = scipy.optimize.minimize(fun=error_function,
                                  method='Nelder-Mead',
                                  x0=[params['E_L1'], params['E_C1'], params['E_J1'], params['E_L2'], params['E_C2'],
                                      params['E_J2'],
                                      params['E_int_chg']],
                                  args=(qubit1, qubit2, params), )
    run_time = time.time() - start_time
    print(res)

    params['E_L1'] = res.x[0]
    params['E_C1'] = res.x[1]
    params['E_J1'] = res.x[2]
    params['E_L2'] = res.x[3]
    params['E_C2'] = res.x[4]
    params['E_J2'] = res.x[5]
    params['E_int_chg'] = res.x[6]

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
    sample = 'Augustus 17'
    subpath = 'Processed Data/Two Coupled Fluxoniums/'
    path = os.path.join(samples_path, sample, subpath)

    filename_in = 'Augustus_17_2.hdf5'
    filename_fit = 'Augustus_17_2.hdf5'

    params = utilities.load_fit(os.path.join(path, filename_in))

    # Initial guess for the qubit parameters.
    # params = {'E_L1': 0.83,  # The inductive energy.
    #           'E_C1': 1.03,  # The charging energy.
    #           'E_J1': 5.35,  # The Josephson energy.
    #           'phi_ext1': [0.],
    #           'num_osc1': 50,
    #           'num_qbt1': 10,
    #           'E_L2': 0.81,  # The inductive energy.
    #           'E_C2': 1.06,  # The charging energy.
    #           'E_J2': 3.93,  # The Josephson energy.
    #           'phi_ext2': [0.],
    #           'num_osc2': 50,
    #           'num_qbt2': 10,
    #           'E_int_chg': 0.27,
    #           'num_tot': 10,
    #           'dphi_ext': 0.,
    #           'ratio_phi': 1,
    #           'error_type': 'absolute_error',
    #           'data_set': 'data1'
    #           }

    # params['E_L'] = 0.8
    # params['E_C'] = 2.26
    # params['E_J'] = 8
    params['E_int_chg'] = 0.27
    params['num_osc1'] = 50
    params['num_osc2'] = 50
    params['num_tot'] = 20
    # params['dphi_ext'] = 0.0012
    params['dphi_ext'] = 0
    # params['num_mod'] = [10]
    # params['frequencies'] = np.array([10.898])
    # params['n_couplings'] = np.array([0.])
    # params['phi_couplings'] = np.array([0.])
    # params['error_type'] = 'absolute_error'
    params['data_set'] = 'data2'

    # phi_ext = np.linspace(0., 1., 51)
    phi_ext = np.linspace(0.43, 0.57, 101)
    # phi_ext2 = np.linspace(0., 1, 51)
    # phi_ext = np.sort(np.concatenate((phi_ext, phi_ext2)))

    utilities.print_params(params)
    params = fit(params)
    utilities.print_params(params)
    utilities.save_fit(os.path.join(path, filename_fit), params)

    # Compute energy levels.
    # params['num_tot'] = 50
    levels = np.zeros((phi_ext.size, params['num_tot']))
    weights = np.zeros((phi_ext.size, params['num_tot'], params['num_qbt1'], params['num_qbt2']))

    params1 = {}
    params1['E_L'] = params['E_L1']
    params1['E_C'] = params['E_C1']
    params1['E_J'] = params['E_J1']
    params1['phi_ext'] = params['phi_ext1']
    params1['num_osc'] = params['num_osc1']
    params1['num_qbt'] = params['num_qbt1']

    params2 = {}
    params2['E_L'] = params['E_L2']
    params2['E_C'] = params['E_C2']
    params2['E_J'] = params['E_J2']
    params2['phi_ext'] = params['phi_ext2']
    params2['num_osc'] = params['num_osc2']
    params2['num_qbt'] = params['num_qbt2']

    qubit1 = Fluxonium(params1)
    qubit2 = Fluxonium(params2)
    # qubit1.phi_ext = 0.5
    # qubit2.phi_ext = 0.5
    # print(qubit1.levels())
    # print(qubit2.levels())
    for idx_phi_ext, value in enumerate(phi_ext):
        qubit1.phi_ext = value
        qubit2.phi_ext = params['ratio_phi'] * value + params['dphi_ext']
        # print(qubit2.phi_ext)
        system = TwoCoupledFluxoniums(qubit1, qubit2, params)
        levels[idx_phi_ext] = system.levels()
        # print(system.weights().shape)
        # print(weights[idx_phi_ext].shape)
        weights[idx_phi_ext] = system.weights()
        sys.stdout.write('Progress: %5.1f%%\r'
                         % (100. * (idx_phi_ext + 1) / len(phi_ext)))

    params['levels'] = levels
    params['weights'] = weights
    params['phi_ext'] = phi_ext

    utilities.save_fit(os.path.join(path, filename_fit), params)

    plot_spectrum.plot_spectrum()
    plot_anchor_points.plot_anchor_points(params['data_set'])

    lines = colored_lines(phi_ext, levels, weights, 0, params['num_tot'], color2QB=True)
    for idx in range(params['num_tot'] - 1):
        plt.gca().add_collection(lines[idx])

    title = ('$E_L1/h=$%.3f GHz, $E_C1/h=$%.3f GHz, $E_J1/h=$%.3f GHz\n'
             '$E_L2/h=$%.3f GHz, $E_C2/h=$%.3f GHz, $E_J2/h=$%.3f GHz\n'
             '$J_C/h=$%.3f GHz' % (
             params['E_L1'], params['E_C1'], params['E_J1'], params['E_L2'], params['E_C2'], params['E_J2'],
             params['E_int_chg']))
    plot_spectrum.label_axes(title)

    plt.show()


if __name__ == '__main__':
    main()
