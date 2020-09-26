import h5py
import numpy as np
from devices import devices

units = {
    'E_L': 'GHz',
    'E_C': 'GHz',
    'E_J': 'GHz',
    'E_B': 'GHz',
    'E_Bs': 'GHz',
    'offset': 'GHz',
    'phi_ext': 'Phi_0',
    'n_g': '2e',
    'f_m': 'GHz',
    'f_m1': 'GHz',
    'f_m2': 'GHz',
    'f_r': 'GHz',
    'g_m_J': 'GHz',
    'g_r_J': 'GHz',
    'g_m_r': 'GHz',
    'g_m1_J': 'GHz',
    'g_m1_r': 'GHz',
    'g_m2_J': 'GHz',
    'g_m2_r': 'GHz',
    'g_m1_m2': 'GHz',
    'frequencies': 'GHz',
    'couplings': 'GHz',
    'n_couplings': 'GHz',
    'phi_couplings': 'GHz',
    'n_cross_couplings': 'GHz',
    'phi_cross_couplings': 'GHz',
    'cutoff_cpl': 'GHz',
    'levels': 'GHz',
    'weights': '',
    'E_C_m': 'GHz',
    'E_L_m': 'GHz',
    'E_C_m1': 'GHz',
    'E_L_m1': 'GHz',
    'E_C_m2': 'GHz',
    'E_L_m2': 'GHz',
    'E_C_r': 'GHz',
    'E_L_r': 'GHz',
    'Ea_C': 'GHz',
    'Ea_g': 'GHz',
    'Eb_g': 'GHz',
    'LJ': 'nH',
    'L': 'nH',
    'L1': 'nH',
    'L2': 'nH',
    'L3': 'nH',
    'M': 'nH',
    'Lr': 'nH',
    'CJ': 'fF',
    'Cm': 'fF',
    'Cm1': 'fF',
    'Cm2': 'fF',
    'Cr': 'fF',
    'mu': '',
    'rho': '',
    'absolute_error': 'GHz',
    'relative_error': '',
    'run_time': 's'
}


def get_device(sample):
    not_found = True
    for device in devices:
        if device['device'] == sample:
            not_found = False
            break
    if not_found:
        raise ValueError("Device '%s' could not be found." % sample)
    return device


def save_fit(filename, data):
    with h5py.File(filename, 'w') as f:
        if 'comment' in data:
            f.attrs['comment'] = data['comment']

        grp_fit = f.create_group('fit')
        grp_units = grp_fit.create_group('units')

        for key in data:
            if key in ['E_J', 'E_C', 'E_L', 'E_J1', 'E_C1', 'E_L1', 'E_J2', 'E_C2', 'E_L2', 'E_B', 'offset',
                       'f_m', 'E_C_m', 'E_L_m',
                       'f_m1', 'E_C_m1', 'E_L_m1',
                       'f_m2', 'E_C_m2', 'E_L_m2',
                       'f_r', 'E_C_r', 'E_L_r',
                       'Ea_C', 'Ea_g', 'Eb_g',
                       'g_m_J', 'g_r_J', 'g_m_r',
                       'g_m1_J', 'g_m1_r', 'g_m2_J', 'g_m2_r', 'g_m1_m2',
                       'LJ', 'L', 'L1', 'L2', 'L3', 'M', 'Lr',
                       'mu', 'rho',
                       'CJ', 'Cm', 'Cm1', 'Cm2', 'Cr',
                       'num_osc', 'num_osc1', 'num_osc2', 'num_chg', 'num_flx', 'num_hrm',
                       'num_qbt', 'num_qbt1', 'num_qbt2', 'num_chn', 'num_res', 'num_tot',
                       'num_chn1', 'num_chn2', 'num_cpl', 'cutoff_cpl',
                       'dphi_ext', 'ratio_phi', 'E_int_chg',
                       'N_chain',
                       'absolute_error', 'relative_error',
                       'error_type', 'coupling_type',
                       'data_set', 'run_time']:
                grp_fit.attrs[key] = data[key]
                if key in units:
                    grp_units.attrs[key] = units[key]

            elif key in ['E_Bs', 'frequencies', 'couplings',
                         'phi_ext', 'phi_ext1', 'phi_ext2', 'n_g', 'levels', 'weights',
                         'num_mod',
                         'phi_couplings', 'n_couplings',
                         'phi_cross_couplings', 'n_cross_couplings']:
                grp_fit.create_dataset(key, data=data[key])
                if key in units:
                    grp_units.attrs[key] = units[key]
            else:
                print("Variable '%s' has not been saved." % key)


def load_fit(filename):
    data = {}
    with h5py.File(filename, 'r') as f:
        grp_fit = f['fit']
        for key in grp_fit.attrs.keys():
            data[key] = grp_fit.attrs[key]

        for key in grp_fit.keys():
            if key != 'units':
                data[key] = np.array(grp_fit[key])

    return data


def print_param(name, value, array=False):
    if name.startswith('num_') or name.startswith('N_'):
        if isinstance(value, (list, np.ndarray)):
            print('%s: %s' % (name, value))
        else:
            print('%s: %d' % (name, value))
        return
    space = ' '
    if name in units:
        unit = units[name]
    else:
        unit = ''
    if name == 'absolute_error':
        value *= 1.e3
        unit = 'MHz'
    if name == 'relative_error':
        value *= 1.e2
        space = ''
        unit = '%'
    if not array:
        if unit != '':
            print('%s: %.4f%s%s' % (name, value, space, unit))
        else:
            print('%s: %.6f' % (name, value))
    elif isinstance(value, np.ndarray) and len(np.shape(value)) == 2:
        print('%s:\n%s%s%s' % (name, value, space, unit))
    elif len(value):
        if unit != '':
            print('%s: %s%s%s' % (name, value, space, unit))
        else:
            print('%s: %s' % (name, array))
    else:
        print('%s: None' % name)


def print_params(params):
    print('===Parameters===')
    for key, value in params.items():
        if isinstance(value, np.ndarray) or isinstance(value, list):
            if len(value) == 1:
                print_param(key, np.array(value).flatten()[0])
            elif key in ('frequencies', 'couplings',
                         'n_couplings', 'phi_couplings',
                         'n_cross_couplings', 'phi_cross_couplings',
                         'num_mod', 'E_Bs'):
                print_param(key, value, array=True)
        elif isinstance(value, list):
            print_param(key, value, array=True)
        elif isinstance(value, str):
            print('%s: %s' % (key, value))
        elif not isinstance(value, list):
            print_param(key, value)
    print('================')
