import numpy as np
import qutip as qt
from qutip.parallel import parallel_map

try:
    from simple_hamiltonians import Fluxonium, HarmonicOscillator
except:
    from .simple_hamiltonians import Fluxonium, HarmonicOscillator


class Model(object):
    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        self._eigvals = None
        self._eigvecs = None
        self._eye_op = None
        self._H_op = None
        self._weights = None
        self._purity = None
        
    @property
    def num_tot(self):
        return self._num_tot

    def _spectrum(self):
        """Eigen-energies and eigenstates in the oscillator basis."""
        if self._eigvals is None or self._eigvecs is None:
            self._eigvals, self._eigvecs = \
                    self._hamiltonian.eigenstates(sparse=True,
                    eigvals=self._num_tot)
        return self._eigvals, self._eigvecs

    def levels(self):
        """Eigen-energies of the coupled system.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenvalues.
        """
        return self._spectrum()[0][:self._num_tot]

    def _check_level(self, level):
        if level < 0 or level >= self._num_tot:
            raise ValueError('The level is out of bounds: 0 and %d.'
                    % self._num_tot)

    def level(self, level):
        """Energy of a single level of the qubit.

        Parameters
        ----------
        level: int
            Qubit level.

        Returns
        -------
        float
            Energy of the level.
        """
        self._check_level(level)
        return self._spectrum()[0][level]

    def frequency(self, level1, level2):
        """Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2` defined
            as the difference of energies. Positive if `level1` < `level2`.
        """
        self._check_level(level1)
        self._check_level(level2)
        return self.level(level2) - self.level(level1)
    
    def states(self):
        """Eigenstates of the system.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenstates.
        """
        return self._spectrum()[1]
        
    def eye(self):
        """Identity operator.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Identity operator.
        """
        if self._eye_op is None:
            self._eye_op = qt.qeye(self._num_tot)
        return self._eye_op
        
    def H(self):
        """Hamiltonian in its eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Hamiltonian operator.
        """
        if self._H_op is None:
            self._H_op = qt.Qobj(np.diag(self.levels()[:self._num_tot]))
        return self._H_op

    def weights(self):
        if self._weights is None:
            evecs = self._spectrum()[1]
            num_qbt = self._fluxonium.num_qbt
            self._weights = np.zeros((self._num_tot, num_qbt))
            for idx in range(self._num_tot):
                w = np.abs(np.array(evecs[idx].data.todense()))**2.
                w.shape = (num_qbt, -1)
                self._weights[idx] = np.sum(w, axis=1)
        return self._weights
        
    def purity(self):
        if self._purity is None:
            evecs = self._spectrum()[1]
            num_qbt = self._fluxonium.num_qbt
            self._purity = np.zeros((self._num_tot, num_qbt))
            for idx in range(self._num_tot):
                w = np.abs(np.array(evecs[idx].data.todense()))**2.
                w.shape = (num_qbt, -1)
                self._purity[idx] = w[:,0]
        return self._purity


class TwoCoupledFluxoniums(Model):
    """Fluxonium Hamiltonian coupled to n harmonic modes."""
    def __init__(self, fluxonium1, fluxonium2, params):

        self._num_tot = fluxonium1.num_qbt * fluxonium2.num_qbt
        
        self._num_wq1 = fluxonium1._num_qbt
        if 'num_wq1' in params:
            self._num_wq1 = int(params['num_wq1'])
        self._num_wq2 = fluxonium2._num_qbt
        if 'num_wq2' in params:
            self._num_wq2 = int(params['num_wq2'])

        num_tot = int(params['num_tot'])
        if num_tot > fluxonium1.num_qbt * fluxonium2.num_qbt:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor(fluxonium1.H(), fluxonium2.eye())
        H += qt.tensor(fluxonium1.eye(), fluxonium2.H())

        if 'E_int_chg' in params:
            H += params['E_int_chg'] * qt.tensor(fluxonium1.n(),
                                                 fluxonium2.n())
        if 'E_int_flx' in params:
            H += params['E_int_flx'] * qt.tensor(fluxonium1.phi(),
                                                 fluxonium2.phi())                                    
        self._hamiltonian = H
        # eigvals, eigvecs = \
        #     self._hamiltonian.eigenstates(sparse=True,
        #                                   eigvals=4)
        # print(eigvals)
        # print(fluxonium1)
        # print(fluxonium2)
        self._fluxonium1 = fluxonium1
        self._fluxonium2 = fluxonium2
        self._reset_cache()
        
    def frequency(self, level1, level2):
        """Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int or tuple
            The qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2`
            defined as the difference of energies.
            Positive if `level1` < `level2`.
        """
        if isinstance(level1, tuple):
            level1 = self.weights(labels=True)[level1[0],level1[1]]
        if isinstance(level2, tuple):
            level2 = self.weights(labels=True)[level2[0],level2[1]]
        self._check_level(level1)
        self._check_level(level2)
        return self.level(level2) - self.level(level1)
        
    def weights(self, labels=False):
        if self._weights is None:
            evecs = self.states()

            qbt1 = self._fluxonium1._num_qbt
            qbt2 = self._fluxonium2._num_qbt
            wq1 = self._num_wq1
            wq2 = self._num_wq2

            weights = np.zeros((self._num_tot, wq1, wq2), dtype=np.complex)
            for idx in range(self._num_tot):
                # print(evecs[idx].data.todense().shape)
                weights[idx] = \
                        evecs[idx].data.todense().reshape(qbt1, qbt2)[:wq1,:wq2]
            self._weights = np.abs(weights)**2.

            # Each label must be used only once.
            self._labels = np.zeros((wq1, wq2), dtype=np.int)
            w = self._weights.copy()
            for idx2 in range(wq2):
                for idx1 in range(wq1):
                    level_idx = np.argmax(w[:,idx1,idx2])
                    w[level_idx] = -np.ones((wq1, wq2))
                    self._labels[idx1][idx2] = level_idx
        if not labels:
            return self._weights
        else:
            return self._labels


class _ReducedHilbertSpace(Model):
    def __init__(self, system1, system2, num_tot):
        H = qt.tensor([system1.H(), system2.eye()])
        H += qt.tensor([system1.eye(), system2.H()])

        if num_tot > system1.H().shape[0] * system2.H().shape[0]:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        self._hamiltonian = H
        self._system1 = system1
        self._system2 = system2
        self._reset_cache()
        
    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        Model._reset_cache(self)
        self._b_ops = None

    def b(self):
        """Annihilation operators in the combined eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            List of annihilation operators.
        """
        if self._b_ops is None:
            self._b_ops = []
            evecs = self.states()
            num_tot = self.num_tot
            num_combined = (self._system1.H().shape[0]
                          * self._system2.H().shape[0])
            evecs_padded = np.hstack([evecs,
                    [0. * evecs[0]] * (num_combined - num_tot)])
            for k, system in enumerate([self._system1, self._system2]):
                bs = system.b()
                if type(bs) != list:
                    bs = [bs]
                for b_idx in bs:
                    if k:
                        b = qt.tensor([self._system1.eye(), b_idx])
                    else:
                        b = qt.tensor([b_idx, self._system2.eye()])

                    b_op = b.transform(evecs_padded)[:num_tot,:num_tot]
                    self._b_ops.append(qt.Qobj(b_op))
        return self._b_ops


class ParallelUncoupledOscillators(Model):
    """Coupled oscillators."""
    def __init__(self, params, pairwise=True):
        if len(params['num_mod']) != len(params['frequencies']):
            raise ValueError('Oscillators are not properly defined.')

        oscillators = []
        for idx in range(len(params['frequencies'])):
            oscillators.append(HarmonicOscillator(
                    frequency=params['frequencies'][idx],
                    num_osc=params['num_mod'][idx]))

        num_tot = int(params['num_cpl'])
        if num_tot > np.prod(params['num_mod']):
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        self._reset_cache()

        if pairwise:
            while len(oscillators) > 1:
                if len(oscillators) == 2:
                    cutoff = num_tot
                else:
                    cutoff = int(2.5 * np.sqrt(num_tot))
                if len(oscillators) % 2 == 0:
                    reduced = []
                    idx = 0
                else:
                    reduced = [oscillators[0]]
                    idx = 1
                reduced = reduced + parallel_map(self._parallel,
                            range(idx, len(oscillators), 2),
                            task_args=(oscillators, cutoff),
                            num_cpus=1)
                oscillators = reduced

            self._hamiltonian = oscillators[0].H()
            self._b_ops = oscillators[0].b()
            if type(self._b_ops) != list:
                self._b_ops = [self._b_ops]
        else:
            if len(oscillators) == 1:
                self._hamiltonian = oscillators[0].H()
                self._b_ops = [oscillators[0].b()]
            else:
                system = oscillators[0]
                for k, osc in enumerate(oscillators[1:]):
                    system = _ReducedHilbertSpace(system, osc, num_tot)
                self._hamiltonian = system.H()
                self._b_ops = system.b()
                
    @classmethod
    def _parallel(cls, k, oscillators, cutoff):
        return _ReducedHilbertSpace(oscillators[k],
                                    oscillators[k+1],
                                    cutoff)

    def b(self):
        """Annihilation operators in the multi-oscillator eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            List of annihilation operators.
        """
        return self._b_ops


class UncoupledOscillators(Model):
    """Coupled oscillators."""
    def __init__(self, params, pairwise=True):
        if len(params['num_mod']) != len(params['frequencies']):
            raise ValueError('Oscillators are not properly defined.')

        oscillators = []
        for idx in range(len(params['frequencies'])):
            oscillators.append(HarmonicOscillator(
                    frequency=params['frequencies'][idx],
                    num_osc=params['num_mod'][idx]))

        num_tot = int(params['num_cpl'])
        if num_tot > np.prod(params['num_mod']):
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        self._reset_cache()

        if pairwise:
            while len(oscillators) > 1:
                if len(oscillators) == 2:
                    cutoff = num_tot
                else:
                    cutoff = int(2.5 * np.sqrt(num_tot))
                if len(oscillators) % 2 == 0:
                    reduced = []
                    idx = 0
                else:
                    reduced = [oscillators[0]]
                    idx = 1
                for k in range(idx, len(oscillators), 2):
                    reduced.append(_ReducedHilbertSpace(oscillators[k],
                                                        oscillators[k+1],
                                                        cutoff))
                oscillators = reduced

            self._hamiltonian = oscillators[0].H()
            self._b_ops = oscillators[0].b()
            if type(self._b_ops) != list:
                self._b_ops = [self._b_ops]
        else:
            if len(oscillators) == 1:
                self._hamiltonian = oscillators[0].H()
                self._b_ops = [oscillators[0].b()]
            else:
                system = oscillators[0]
                for k, osc in enumerate(oscillators[1:]):
                    system = _ReducedHilbertSpace(system, osc, num_tot)
                self._hamiltonian = system.H()
                self._b_ops = system.b()

    def b(self):
        """Annihilation operators in the multi-oscillator eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            List of annihilation operators.
        """
        return self._b_ops


class CoupledOscillators(Model):
    """Coupled oscillators."""
    def __init__(self, params):
        if len(params['num_mod']) != len(params['frequencies']):
            raise ValueError('Oscillators are not properly defined.')

        oscillators = []
        for idx in range(len(params['frequencies'])):
            oscillators.append(HarmonicOscillator(
                    frequency=params['frequencies'][idx],
                    num_osc=params['num_mod'][idx]))
    
        if 'n_cross_couplings' in params:
            dim = len(oscillators)
            if params['n_cross_couplings'].shape != (dim, dim):
                raise ValueError('The charge cross coupling matrix'
                        ' is not properly defined.')
        
        if 'phi_cross_couplings' in params:
            dim = len(oscillators)
            if params['phi_cross_couplings'].shape != (dim, dim):
                raise ValueError('The flux cross coupling matrix'
                        ' is not properly defined.')

        osc_total = 1
        for osc in oscillators:
            osc_total *= osc.num_osc
        num_tot = int(params['num_cpl'])
        if num_tot > osc_total:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        for idx1, osc1 in enumerate(oscillators):
            array = []
            for idx2, osc2 in enumerate(oscillators):
                if idx1 == idx2:
                    array.append(osc1.H())
                else:
                    array.append(osc2.eye())
            if idx1 == 0:
                H = qt.tensor(*array)
            else:
                H += qt.tensor(*array)

        for idx1, osc1 in enumerate(oscillators):
            for idx2, osc2 in enumerate(oscillators):
                if idx2 > idx1:
                    if 'n_cross_couplings' in params:
                        array = []
                        for idx3 in range(len(oscillators)):
                            if idx3 == idx1:
                                array.append(osc1.b() - osc1.b().dag())
                            elif idx3 == idx2:
                                array.append(osc2.b() - osc2.b().dag())
                            else:
                                array.append(oscillators[idx3].eye())
                        H += (params['n_cross_couplings'][idx1, idx2] *
                                qt.tensor(*array))

                    if 'phi_cross_couplings' in params:
                        array = []
                        for idx3 in range(len(oscillators)):
                            if idx3 == idx1:
                                array.append(osc1.b() + osc1.b().dag())
                            elif idx3 == idx2:
                                array.append(osc2.b() + osc2.b().dag())
                            else:
                                array.append(oscillators[idx3].eye())
                        H += (params['phi_cross_couplings'][idx1,idx2] *
                                qt.tensor(*array))

        self._hamiltonian = H
        self._oscillators = oscillators
        self._reset_cache()
        
    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        Model._reset_cache(self)
        self._b_ops = None

    def b(self):
        """Annihilation operators in the multi-oscillator eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            List of annihilation operators.
        """
        if self._b_ops is None:
            self._b_ops = []
            evecs = self.states()
            num_tot = self.num_tot
            for idx1, osc1 in enumerate(self._oscillators):
                array = []
                for idx2, osc2 in enumerate(self._oscillators):
                    if idx1 == idx2:
                        array.append(osc1.b())
                    else:
                        array.append(osc2.eye())
                b = qt.tensor(*array)
                
                evecs_padded = np.hstack([evecs,
                        [0. * evecs[0]] * (b.shape[0] - num_tot)])
                b_op = b.transform(evecs_padded)[:num_tot,:num_tot]
                self._b_ops.append(qt.Qobj(b_op))
        return self._b_ops


class FluxoniumCoupledToOscillators(Model):
    """Fluxonium Hamiltonian coupled to harmonic modes."""
    def __init__(self, fluxonium, oscillators, params):
        if ('n_couplings' in params and
                len(params['frequencies']) != len(params['n_couplings'])):
            raise ValueError('The number of oscillators should be equal'
                    ' to the number of charge couplings.')
        
        if ('phi_couplings' in params and
                len(params['frequencies']) != len(params['phi_couplings'])):
            raise ValueError('The number of oscillators should be equal'
                    ' to the number of flux couplings.')
        
        num_tot = int(params['num_tot'])
        if num_tot > fluxonium.num_qbt * oscillators.num_tot:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor([fluxonium.H(), oscillators.eye()])
        H += qt.tensor([fluxonium.eye(), oscillators.H()])

        for idx, b in enumerate(oscillators.b()):
            if 'n_couplings' in params:
                op = -1.j * params['n_couplings'][idx] * (b - b.dag())
                H += qt.tensor([fluxonium.n(), op])
            if 'phi_couplings' in params:
                op = params['phi_couplings'][idx] * (b + b.dag())
                H += qt.tensor([fluxonium.phi(), op])

        self._hamiltonian = (H + H.dag()) / 2.
        self._fluxonium = fluxonium
        self._coupled_oscillators = oscillators
        self._reset_cache()
        
    def n_ij(self, level1, level2):
        """The charge matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        complex
            Matrix element of the charge operator.
        """
        self._check_level(level1)
        self._check_level(level2)
        evecs = self._spectrum()[1]
        
        op = qt.tensor([self._fluxonium.n(),
                        self._coupled_oscillators.eye()])

        return op.matrix_element(evecs[level1].dag(), evecs[level2])
        
    def phi_ij(self, level1, level2):
        """The flux matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        complex
            Matrix element of the flux operator.
        """
        self._check_level(level1)
        self._check_level(level2)
        evecs = self._spectrum()[1]
        
        op = qt.tensor([self._fluxonium.phi(),
                        self._coupled_oscillators.eye()])

        return op.matrix_element(evecs[level1].dag(), evecs[level2])
        
    def sin_half_phi_ij(self, level1, level2):
        """The sine of half flux matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        complex
            Matrix element of the sine of half flux operator.
        """
        self._check_level(level1)
        self._check_level(level2)
        evecs = self._spectrum()[1]

        op = qt.tensor([(.5 * self._fluxonium.phi()
                  + np.pi * self._fluxonium.phi_ext).sinm(),
                  self._coupled_oscillators.eye()])

        return op.matrix_element(evecs[level1].dag(), evecs[level2])

    def b_ij(self, level1, level2, mode=0):
        """The bosonic operator matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.
            
        mode: int
            Index of the harmonic mode.

        Returns
        -------
        complex
            Matrix element of the bosonic operator.
        """
        self._check_level(level1)
        self._check_level(level2)
        evecs = self._spectrum()[1]
        
        op = qt.tensor([self._fluxonium.eye(),
                        self._coupled_oscillators.b()[mode]])

        return op.matrix_element(evecs[level1].dag(), evecs[level2])


class FluxoniumChargeCoupledToResonator(Model):
    """Fluxonium Hamiltonian charge coupled to a resonator mode."""
    def __init__(self, fluxonium, params):
        resonator = HarmonicOscillator(frequency=params['f_r'],
                                       num_osc=params['num_res'])
        
        osc_tot = fluxonium.num_qbt * params['num_res']
        num_tot = int(params['num_tot'])
        if num_tot > osc_tot:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor(fluxonium.H(), resonator.eye())
        H += qt.tensor(fluxonium.eye(), resonator.H())

        H += (1.j * params['g_r_J']
                  * qt.tensor(fluxonium.n(),
                              resonator.b() - resonator.b().dag()))
        self._hamiltonian = H
        self._fluxonium = fluxonium
        self._reset_cache()


class FluxoniumChargeCoupledToOneModeAndResonator(Model):
    """
    Fluxonium Hamiltonian charge-coupled to a chain and a resonator
    modes.
    """
    def __init__(self, fluxonium, params):
        mode = HarmonicOscillator(frequency=params['f_m'],
                            num_osc=params['num_chn'])
        resonator = HarmonicOscillator(frequency=params['f_r'],
                                 num_osc=params['num_res'])
        
        osc_tot = fluxonium.num_qbt * params['num_chn'] * params['num_res']
        num_tot = int(params['num_tot'])
        if num_tot > osc_tot:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor(fluxonium.H(), mode.eye(), resonator.eye())
        H += qt.tensor(fluxonium.eye(), mode.H(), resonator.eye())
        H += qt.tensor(fluxonium.eye(), mode.eye(), resonator.H())

        mode_n = mode.b() - mode.b().dag()
        resonator_n = resonator.b() - resonator.b().dag()
        H += 1.j * params['g_m_J'] * qt.tensor(fluxonium.n(),
                                               mode_n,
                                               resonator.eye())
        H += 1.j * params['g_r_J'] * qt.tensor(fluxonium.n(),
                                               mode.eye(),
                                               resonator_n)
        H += params['g_m_r'] * qt.tensor(fluxonium.eye(),
                                         mode_n,
                                         resonator_n)
        self._hamiltonian = H
        self._fluxonium = fluxonium
        self._reset_cache()

        
class FluxoniumFluxCoupledToOneModeAndResonator(Model):
    """
    Fluxonium Hamiltonian flux-coupled to a chain and a resonator
    modes.
    """
    def __init__(self, fluxonium, params):
        mode = HarmonicOscillator(frequency=params['f_m'],
                            num_osc=params['num_chn'])
        resonator = HarmonicOscillator(frequency=params['f_r'],
                                 num_osc=params['num_res'])
        
        osc_tot = fluxonium.num_qbt * params['num_chn'] * params['num_res']
        num_tot = int(params['num_tot'])
        if num_tot > osc_tot:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor(fluxonium.H(), mode.eye(), resonator.eye())
        H += qt.tensor(fluxonium.eye(), mode.H(), resonator.eye())
        H += qt.tensor(fluxonium.eye(), mode.eye(), resonator.H())

        mode_phi = mode.b() + mode.b().dag()
        resonator_phi = resonator.b() + resonator.b().dag()
        H += params['g_m_J'] * qt.tensor(fluxonium.phi(),
                                         mode_phi,
                                         resonator.eye())
        H += params['g_r_J'] * qt.tensor(fluxonium.phi(),
                                         mode.eye(),
                                         resonator_phi)
        H += params['g_m_r'] * qt.tensor(fluxonium.eye(),
                                         mode_phi,
                                         resonator_phi)
        self._hamiltonian = H
        self._fluxonium = fluxonium
        self._reset_cache()


class FluxoniumFluxCoupledToTwoModesAndResonator(Model):
    """
    Fluxonium Hamiltonian flux-coupled to two chain and a resonator
    modes.
    """
    def __init__(self, fluxonium, params):
        mode1 = HarmonicOscillator(frequency=params['f_m1'],
                            num_osc=params['num_chn1'])
        mode2 = HarmonicOscillator(frequency=params['f_m2'],
                            num_osc=params['num_chn2'])
        resonator = HarmonicOscillator(frequency=params['f_r'],
                                 num_osc=params['num_res'])
        
        osc_tot = (fluxonium.num_qbt * params['num_chn1']
                   * params['num_chn2'] * params['num_res'])
        num_tot = int(params['num_tot'])
        if num_tot > osc_tot:
            raise ValueError('The number of levels is too high.')
        self._num_tot = num_tot

        H = qt.tensor(fluxonium.H(),
                mode1.eye(), mode2.eye(), resonator.eye())
        H += qt.tensor(fluxonium.eye(),
                mode1.H(), mode2.eye(), resonator.eye())
        H += qt.tensor(fluxonium.eye(),
                mode1.eye(), mode2.H(), resonator.eye())
        H += qt.tensor(fluxonium.eye(),
                mode1.eye(), mode2.eye(), resonator.H())

        mode1_phi = mode1.b() + mode1.b().dag()
        mode2_phi = mode2.b() + mode2.b().dag()
        resonator_phi = resonator.b() + resonator.b().dag()
        H += params['g_m1_J'] * qt.tensor(fluxonium.phi(),
                                          mode1_phi,
                                          mode2.eye(),
                                          resonator.eye())
        H += params['g_m2_J'] * qt.tensor(fluxonium.phi(),
                                          mode1.eye(),
                                          mode2_phi,
                                          resonator.eye())
        H += params['g_r_J'] * qt.tensor(fluxonium.phi(),
                                         mode1.eye(),
                                         mode2.eye(),
                                         resonator_phi)
        H += params['g_m1_r'] * qt.tensor(fluxonium.eye(),
                                          mode1_phi,
                                          mode2.eye(),
                                          resonator_phi)
        H += params['g_m2_r'] * qt.tensor(fluxonium.eye(),
                                          mode1.eye(),
                                          mode2_phi,
                                          resonator_phi)
        H += params['g_m1_m2'] * qt.tensor(fluxonium.eye(),
                                           mode1_phi,
                                           mode2_phi,
                                           resonator.eye())
        self._hamiltonian = H
        self._fluxonium = fluxonium
        self._reset_cache()
