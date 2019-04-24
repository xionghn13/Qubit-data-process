import numpy as np
import qutip as qt
from scipy.special import eval_hermite as hpoly


class HarmonicOscillator(object):
    """Harmonic oscillator Hamiltonian."""
    def __init__(self, frequency, num_osc=20):
        self._frequency = 0.
        self._num_osc = 0
        self.frequency = frequency
        self.num_osc = num_osc

    def __str__(self):
        units = 'GHz'
        return ('Harmonic oscillator f = %.4f %s.'
                % (self.frequency, units))

    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, value):
        if value != self._frequency:
            self._frequency = value
            self._reset_cache()

    @property
    def num_osc(self):
        return self._num_osc

    @num_osc.setter
    def num_osc(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of levels must be positive.')
        if value != self._num_osc:
            self._num_osc = value
            self._reset_cache()

    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        self._H_op = None
        self._b_op = None
        self._eye_op = None

    def eye(self):
        """
        Identity operator.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Identity operator.
        """
        if self._eye_op is None:
            self._eye_op = qt.qeye(self._num_osc)
        return self._eye_op

    def b(self):
        """Annihilation operator."""
        if self._b_op is None:
            self._b_op = qt.destroy(self._num_osc)
        return self._b_op

    def H(self):
        """Hamiltonian."""
        if self._H_op is None:
            self._H_op = self._frequency * self.b().dag() * self.b()
        return self._H_op


class Transmon(object):
    """Transmon qubit."""
    def __init__(self, params):
        for key in ['E_C', 'E_J', 'n_g', 'num_chg', 'num_qbt']:
            if key not in params:
                raise ValueError('%s should be specified.' % key)

        self._E_C = 0.
        self._E_J = 0.
        self._n_g = 0.
        self._num_chg = 0
        self._num_qbt = 0

        # Specify the charging energy.
        self.E_C = params['E_C']
        # Specify the Josephson energy.
        self.E_J = params['E_J']
        # Specify the offset charge defined as a fraction of the single
        # Cooper pair charge, i.e. 2e.
        if np.array(params['n_g']).size > 1:
            n_g = params['n_g'][0]
        else:
            n_g = params['n_g']
        self.n_g = n_g
        # Specify the number of states in the charge basis.
        self.num_chg = params['num_chg']
        # Specify the number of states in the qubit basis.
        self.num_qbt = params['num_qbt']

    def __str__(self):
        units = 'GHz'
        return ('Transmon: E_C = %.4f %s, E_J = %.4f %s, offset charge '
                'n_g / (2e) = %.4f.'
                % (self.E_C, units, self.E_J, units, self.n_g))

    @property
    def E_C(self):
        return self._E_C

    @E_C.setter
    def E_C(self, value):
        if value <= 0:
            raise ValueError('Charging energy must be positive.')
        if value != self._E_C:
            self._E_C = value
            self._reset_cache()

    @property
    def E_J(self):
        return self._E_J

    @E_J.setter
    def E_J(self, value):
        if value <= 0:
            raise ValueError('Josephson energy must be positive.')
        if value != self._E_J:
            self._E_J = value
            self._reset_cache()

    @property
    def n_g(self):
        return self._n_g

    @n_g.setter
    def n_g(self, value):
        # n_g is defined as a fraction of 2e.
        if value != self._n_g:
            self._n_g = value
            self._reset_cache()

    @property
    def num_chg(self):
        return self._num_chg

    @num_chg.setter
    def num_chg(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of charge states '
                    'must be positive.')
        if value != self._num_chg:
            self._num_chg = value
            self._reset_cache()
        
    @property
    def num_qbt(self):
        return self._num_qbt

    @num_qbt.setter
    def num_qbt(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of qubit levels '
                    'must be positive.')
        if value > self._num_chg:
            raise ValueError('The number of qubit levels exceeds '
                    'the number of charge states.')
        if value != self._num_qbt:
            self._num_qbt = value
            self._reset_cache()

    def _reset_cache(self):
        """Reset the cached data that have already been calculated."""
        self._eigvals = None
        self._eigvecs = None
        self._H_op = None
        self._eye_op = None
        self._n_op = None
        self._phi_op = None

    def _phi_chg(self):
        """Flux (phase) operator in the charge basis."""
        raise NotImplementedError()
        
    def _n(self):
        num_chg = self.num_chg
        low = np.ceil(self.n_g - np.double(num_chg) / 2.)
        return np.linspace(low, low + num_chg - 1, num_chg)

    def _n_chg(self):
        """Charge operator in the charge basis."""
        return qt.Obj(np.diag(self._n))

    def _H_chg(self):
        """Qubit Hamiltonian in the charge basis."""
        H = np.diag(4. * self.E_C * (self._n() - self.n_g)**2.)
        off_diag = -.5 * self.E_J * np.ones(self.num_chg - 1)
        H += np.diag(off_diag, 1)
        H += np.diag(off_diag, -1)
        return qt.Qobj(H)

    def _spectrum_chg(self):
        """Eigen-energies and eigenstates in the charge basis."""
        if self._eigvals is None or self._eigvecs is None:
            self._eigvals, self._eigvecs = self._H_chg().eigenstates()
        return self._eigvals, self._eigvecs

    def levels(self):
        """
        Eigen-energies of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenvalues.
        """
        return self._spectrum_chg()[0][:self._num_qbt]
        
    def states(self):
        """
        Eigenstates of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenstates.
        """
        return self._spectrum_chg()[1][:self._num_qbt]
        
    def _check_level(self, level):
        if level < 0 or level >= self._num_qbt:
            raise ValueError('The level index is out of bounds.')

    def level(self, level):
        """
        Energy of a single level of the qubit.

        Parameters
        ----------
        level : int
            Qubit level starting from zero.

        Returns
        -------
        float
            Energy of the level.
        """
        self._check_level(level)
        return self._spectrum_chg()[0][level]

    def frequency(self, level1, level2):
        """
        Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2`
            defined as the difference of energies. Positive
            if `level1` < `level2`.
        """
        self._check_level(level1)
        self._check_level(level2)
        return self.level(level2) - self.level(level1)

    def eye(self):
        """
        Identity operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Identity operator.
        """
        if self._eye_op is None:
            self._eye_op = qt.qeye(self._num_qbt)
        return self._eye_op

    def H(self):
        """
        Qubit Hamiltonian in its eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Hamiltonian operator.
        """
        if self._H_op is None:
            self._H_op = qt.Qobj(np.diag(self.levels()[:self._num_qbt]))
        return self._H_op

    def phi(self):
        """
        Generalized-flux operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Flux operator.
        """
        raise NotImplementedError()

    def n(self):
        """
        Charge operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Charge operator.
        """
        if self._n_op is None:
            evecs = self._spectrum_chg()[1]
            n_op = self._n_chg().transform(evecs)[:self._num_qbt,
                                                  :self._num_qbt]
            self._n_op = qt.Qobj(n_op)
        return self._n_op

    def phi_ij(self, level1, level2):
        """
        Flux matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        complex
            Matrix element of the flux operator.
        """
        raise NotImplementedError()

    def n_ij(self, level1, level2):
        """
        Charge matrix element between two eigenstates.

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
        evecs = self._spectrum_chg()[1]
        return self._n_chg().matrix_element(evecs[level1].dag(),
                                            evecs[level2])


class Blochnium(object):
    """Blochnium qubit."""
    def __init__(self, params):
        for key in ['E_L', 'E_Bs', 'phi_ext', 'num_flx', 'num_qbt']:
            if key not in params:
                raise ValueError('%s should be specified.' % key)

        self._E_L = 0.
        self._E_Bs = np.array([0.])
        self._phi_ext = 0.
        self._num_flx = 0
        self._num_qbt = 0

        # Specify the inductive energy.
        self.E_L = params['E_L']
        # Specify the Bloch band energy.
        self.E_Bs = params['E_Bs']
        # Specify the phi_ext defined as a fraction of Phi_0.
        if np.array(params['phi_ext']).size > 1:
            phi_ext = params['phi_ext'][0]
        else:
            phi_ext = params['phi_ext']
        self.phi_ext = phi_ext
        # Specify the number of states in the flux basis.
        self.num_flx = params['num_flx']
        # Specify the number of states in the qubit basis.
        self.num_qbt = params['num_qbt']

    def __str__(self):
        units = 'GHz'
        return ('Blochnium: E_L = %.4f %s, E_Bs = %s %s, '
                'external flux shift Phi_ext/Phi_0 = %.4f.'
                % (self.E_L, units, str(self.E_Bs), units, self.phi_ext))

    @property
    def E_L(self):
        return self._E_L

    @E_L.setter
    def E_L(self, value):
        if value <= 0:
            raise ValueError('Inductive energy must be positive.')
        if value != self._E_L:
            self._E_L = value
            self._reset_cache()

    @property
    def E_Bs(self):
        return self._E_Bs

    @E_Bs.setter
    def E_Bs(self, value):
        if type(self._E_Bs) != type(value) or \
                np.array(self._E_Bs).size != np.array(value).size or \
                np.any(np.array(self._E_Bs) != np.array(value)):
            self._E_Bs = value
            self._reset_cache()

    @property
    def phi_ext(self):
        return self._phi_ext

    @phi_ext.setter
    def phi_ext(self, value):
        # phi_ext is defined as a fraction of Phi_0.
        if value != self._phi_ext:
            self._phi_ext = value
            self._reset_cache()

    @property
    def num_flx(self):
        return self._num_flx

    @num_flx.setter
    def num_flx(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of flux states '
                    'must be positive.')
        if value != self._num_flx:
            self._num_flx = value
            self._reset_cache()
        
    @property
    def num_qbt(self):
        return self._num_qbt

    @num_qbt.setter
    def num_qbt(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of qubit levels '
                    'must be positive.')
        if value > self._num_flx:
            raise ValueError('The number of qubit levels exceeds '
                    'the number of flux states.')
        if value != self._num_qbt:
            self._num_qbt = value
            self._reset_cache()

    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        self._eigvals = None
        self._eigvecs = None
        self._H_op = None
        self._eye_op = None
        self._n_op = None
        self._phi_op = None
    
    def _phi(self):
        num_flx = self.num_flx
        low = np.ceil(self.phi_ext - np.double(num_flx) / 2.)
        return 2. * np.pi * np.linspace(low, low + num_flx - 1, num_flx)

    def _phi_flx(self):
        """Flux (phase) operator in the flux basis."""
        return qt.Obj(np.diag(self._phi()))

    def _n_flx(self):
        """Charge operator in the flux basis."""
        raise NotImplementedError()

    def _H_flx(self):
        """Qubit Hamiltonian in the flux basis."""
        H = np.diag(.5 * self.E_L * (self._phi() - 2. * np.pi * self.phi_ext)**2.)
        for k in range(len(self.E_Bs)):
            off_diag = .5 * self.E_Bs[k] * np.ones(self.num_flx - k - 1)
            H += np.diag(off_diag,  k + 1)
            H += np.diag(off_diag, -k - 1)
        return qt.Qobj(H)

    def _spectrum_flx(self):
        """Eigen-energies and eigenstates in the flux basis."""
        if self._eigvals is None or self._eigvecs is None:
            self._eigvals, self._eigvecs = self._H_flx().eigenstates()
        return self._eigvals, self._eigvecs

    def levels(self):
        """
        Eigen-energies of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenvalues.
        """
        return self._spectrum_flx()[0][:self._num_qbt]
        
    def states(self):
        """
        Eigenstates of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenstates.
        """
        return self._spectrum_flx()[1][:self._num_qbt]
        
    def _check_level(self, level):
        if level < 0 or level >= self._num_qbt:
            raise ValueError('The level index is out of bounds.')

    def level(self, level):
        """
        Energy of a single level of the qubit.

        Parameters
        ----------
        level : int
            Qubit level starting from zero.

        Returns
        -------
        float
            Energy of the level.
        """
        self._check_level(level)
        return self._spectrum_flx()[0][level]

    def frequency(self, level1, level2):
        """
        Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2`
            defined as the difference of energies. Positive
            if `level1` < `level2`.
        """
        self._check_level(level1)
        self._check_level(level2)
        return self.level(level2) - self.level(level1)

    def eye(self):
        """
        Identity operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Identity operator.
        """
        if self._eye_op is None:
            self._eye_op = qt.qeye(self._num_qbt)
        return self._eye_op

    def H(self):
        """
        Qubit Hamiltonian in its eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Hamiltonian operator.
        """
        if self._H_op is None:
            self._H_op = qt.Qobj(np.diag(self.levels()[:self._num_qbt]))
        return self._H_op

    def phi(self):
        """
        Generalized-flux operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Flux operator.
        """
        if self._phi_op is None:
            evecs = self._spectrum_flx()[1]
            phi_op = self._phi_flx().transform(evecs)[:self._num_qbt,
                                                      :self._num_qbt]
            self._phi_op = qt.Qobj(phi_op)
        return self._phi_op

    def n(self):
        """
        Charge operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Charge operator.
        """
        raise NotImplementedError()

    def phi_ij(self, level1, level2):
        """
        Flux matrix element between two eigenstates.

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
        evecs = self._spectrum_flx()[1]
        return self._phi_flx().matrix_element(evecs[level1].dag(),
                                              evecs[level2])

    def n_ij(self, level1, level2):
        """
        Charge matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        complex
            Matrix element of the charge operator.
        """
        raise NotImplementedError()

    def wavefunc(self, level):
        """
        Flux-dependent wavefunction of an eigenstate.
        
        Parameters
        ----------
        level : int
            Qubit level.
            
        Returns
        -------
        phi : ndarray
            Phi basis (in radians).
        wf : ndarray
            Wavefunction in phi basis.
        """
        self._check_level(level)
        evecs = self._spectrum_flx()[1]
        return (self._phi(), np.array(evecs[level].data.todense())[:,0])


class Fluxonium(object):
    """Fluxonium qubit."""
    def __init__(self, params):
        for key in ['E_L', 'E_C', 'E_J', 'phi_ext', 'num_osc', 'num_qbt']:
            if key not in params:
                raise ValueError('%s should be specified.' % key)

        self._E_L = 0.
        self._E_C = 0.
        self._E_J = 0.
        self._phi_ext = 0.
        self._num_osc = 0
        self._num_qbt = 0

        # Specify the inductive energy.
        self.E_L = params['E_L']
        # Specify the charging energy.
        self.E_C = params['E_C']
        # Specify the Josephson energy.
        self.E_J = params['E_J']
        # Specify the phi_ext defined as a fraction of Phi_0.
        self.phi_ext = np.array(params['phi_ext'])[0]
        # Specify the number of states in the oscillator basis.
        self.num_osc = params['num_osc']
        # Specify the number of states in the qubit basis.
        self.num_qbt = params['num_qbt']

    def __str__(self):
        units = 'GHz'
        return ('Fluxonium: E_L = %.4f %s GHz, E_C = %.4f %s, '
                'and E_J = %.4f %s, external flux shift '
                'Phi_ext/Phi_0 = %.4f.'
                % (self.E_L, units, self.E_C, units,
                   self.E_J, units, self.phi_ext))

    @property
    def E_L(self):
        return self._E_L

    @E_L.setter
    def E_L(self, value):
        if value <= 0:
            raise ValueError('Inductive energy must be positive.')
        if value != self._E_L:
            self._E_L = value
            self._reset_cache()

    @property
    def E_C(self):
        return self._E_C

    @E_C.setter
    def E_C(self, value):
        if value <= 0:
            raise ValueError('Charging energy must be positive.')
        if value != self._E_C:
            self._E_C = value
            self._reset_cache()

    @property
    def E_J(self):
        return self._E_J

    @E_J.setter
    def E_J(self, value):
        if value <= 0:
            raise ValueError('Josephson energy must be positive.')
        if value != self._E_J:
            self._E_J = value
            self._reset_cache()

    @property
    def phi_ext(self):
        return self._phi_ext

    @phi_ext.setter
    def phi_ext(self, value):
        # phi_ext is defined as a fraction of Phi_0.
        if value != self._phi_ext:
            self._phi_ext = value
            self._reset_cache()

    @property
    def num_osc(self):
        return self._num_osc

    @num_osc.setter
    def num_osc(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of oscillator levels '
                    'must be positive.')
        if value != self._num_osc:
            self._num_osc = value
            self._reset_cache()
        
    @property
    def num_qbt(self):
        return self._num_qbt

    @num_qbt.setter
    def num_qbt(self, value):
        value = int(value)
        if value <= 0:
            raise ValueError('The number of qubit levels '
                    'must be positive.')
        if value > self._num_osc:
            raise ValueError('The number of qubit levels exceeds '
                    'the number of oscillator levels.')
        if value != self._num_qbt:
            self._num_qbt = value
            self._reset_cache()

    def _reset_cache(self):
        """Reset cached data that have already been calculated."""
        self._eigvals = None
        self._eigvecs = None
        self._H_op = None
        self._eye_op = None
        self._n_op = None
        self._phi_op = None

    def _phi_osc(self):
        """Flux (phase) operator in the oscillator basis."""
        return (8. * self.E_C / self.E_L)**.25 * qt.position(self.num_osc)

    def _n_osc(self):
        """Charge operator in the oscillator basis."""
        return (self.E_L / (8. * self.E_C))**.25 * qt.momentum(self.num_osc)

    def _H_osc(self):
        """Qubit Hamiltonian in the oscillator basis."""
        E_C = self.E_C
        E_L = self.E_L
        E_J = self.E_J
        phi = self._phi_osc()
        n = self._n_osc()
        delta_phi = phi - 2. * np.pi * self.phi_ext
        return 4. * E_C * n**2. + .5 * E_L * phi**2. - E_J * delta_phi.cosm()

    def _spectrum_osc(self):
        """Eigen-energies and eigenstates in the oscillator basis."""
        if self._eigvals is None or self._eigvecs is None:
            self._eigvals, self._eigvecs = self._H_osc().eigenstates()
        return self._eigvals, self._eigvecs

    def levels(self):
        """
        Eigen-energies of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenvalues.
        """
        return self._spectrum_osc()[0][:self._num_qbt]
        
    def states(self):
        """
        Eigenstates of the qubit.

        Parameters
        ----------
        None.

        Returns
        -------
        numpy.ndarray
            Array of eigenstates.
        """
        return self._spectrum_osc()[1][:self._num_qbt]
        
    def _check_level(self, level):
        if level < 0 or level >= self._num_qbt:
            raise ValueError('The level index is out of bounds.')

    def level(self, level):
        """
        Energy of a single level of the qubit.

        Parameters
        ----------
        level : int
            Qubit level, starting from zero.

        Returns
        -------
        float
            Energy of the qubit level.
        """
        self._check_level(level)
        return self._spectrum_osc()[0][level]

    def frequency(self, level1, level2):
        """
        Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2`
            defined as the difference of energies. Positive
            if `level1` < `level2`.
        """
        self._check_level(level1)
        self._check_level(level2)
        return self.level(level2) - self.level(level1)

    def eye(self):
        """
        Identity operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Identity operator.
        """
        if self._eye_op is None:
            self._eye_op = qt.qeye(self._num_qbt)
        return self._eye_op

    def H(self):
        """
        Qubit Hamiltonian in its eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Hamiltonian operator.
        """
        if self._H_op is None:
            self._H_op = qt.Qobj(np.diag(self.levels()[:self._num_qbt]))
        return self._H_op

    def phi(self):
        """
        Generalized-flux operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Flux operator.
        """
        if self._phi_op is None:
            evecs = self._spectrum_osc()[1]
            phi_op = self._phi_osc().transform(evecs)[:self._num_qbt,
                                                      :self._num_qbt]
            self._phi_op = qt.Qobj(phi_op)
        return self._phi_op

    def n(self):
        """Charge operator in the qubit eigenbasis.

        Parameters
        ----------
        None.

        Returns
        -------
        :class:`qutip.Qobj`
            Charge operator.
        """
        if self._n_op is None:
            evecs = self._spectrum_osc()[1]
            n_op = self._n_osc().transform(evecs)[:self._num_qbt,
                                                  :self._num_qbt]
            self._n_op = qt.Qobj(n_op)
        return self._n_op

    def phi_ij(self, level1, level2):
        """
        The flux matrix element between two eigenstates.

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
        evecs = self._spectrum_osc()[1]
        return self._phi_osc().matrix_element(evecs[level1].dag(),
                                              evecs[level2])

    def n_ij(self, level1, level2):
        """
        Charge matrix element between two eigenstates.

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
        evecs = self._spectrum_osc()[1]
        return self._n_osc().matrix_element(evecs[level1].dag(),
                                            evecs[level2])

    def potential(self, phi):
        """Returns the fluxonium potential as a function of flux.
        
        Parameters
        ----------
        phi : float or ndarray
            Value(s) of the flux variable (in radians).
            
        Returns
        -------
        float or ndarray
            Potential energy.
        """
        return (.5 * self.E_L * phi**2 -
                self.E_J * np.cos(phi - 2. * np.pi * self.phi_ext))

    def wavefunc(self, level, phi_n, basis='flux'):
        """
        Wavefunction of an eigenstate in either flux or charge basis.
        
        Parameters
        ----------
        level : int
            Qubit level.
        phi_n : float or ndarray.
            Value(s) of the flux variable (in radians) or value(s) of
            the charge variable (in Cooper pair charge units).
        basis : str
            Use 'flux' for flux representation (defualt).
            Use 'charge' for charge representation.
            
        Returns
        -------
        ndarray
            Wavefunction in flux or charge basis.
            
        Notes
        -----
        The following equations could be obtained from
        http://wiki.physics.fsu.edu/wiki/index.php/Harmonic_Oscillator_Spectrum_and_Eigenstates
        assuming hbar -> 1 and m * omega -> sqrt(E_L / (8 * E_C)).
        """
        self._check_level(level)
        
        def ho_wf(xp, lvl, ratio):
            # Calcultes wave functions of the harmonic oscillator.
            coeff = (2.**lvl * np.math.factorial(lvl)
                     * np.sqrt(np.pi) * ratio)**(-.5)
            coeff *= np.exp(-.5 * (xp / ratio)**2.) 
            return coeff * hpoly(lvl, xp / ratio)

        wf = np.zeros_like(phi_n, dtype=np.complex)
        evecs = self._spectrum_osc()[1]
        if basis == 'flux':
            ratio = (8. * self.E_C / self.E_L)**.25
            for lvl_idx in range(np.min([150, self.num_osc])):
                coeff = evecs[level][lvl_idx,0]
                wf += coeff * ho_wf(phi_n, lvl_idx, ratio)
        elif basis == 'charge':
            ratio = (self.E_L / (8. * self.E_C))**.25
            for lvl_idx in range(np.min([150, self.num_osc])):
                coeff = evecs[level][lvl_idx,0]
                wf += (-1.j)**lvl_idx * coeff * ho_wf(phi_n, lvl_idx, ratio)
        else:
            raise NotImplementedError("Unknown basis '%s'." % basis)
        return wf

    def overlap_function(self, level1, level2, span=5., N=5000):
        """The overlap function.

        Parameters
        ----------
        level1, level2 : int
            Qubit levels.
            
        span : float
            Integration span: the overlap function is computed from
            -2pi*span to 2pi*span.
        N : int
            Number of points for numerical integration.

        Returns
        -------
        complex
            Matrix element of the flux operator.
            
        References
        ----------
        https://journals.aps.org/prb/pdf/10.1103/PhysRevB.85.024521
        [see Eq. (7)].
        """
        self._check_level(level1)
        self._check_level(level2)

        phi_max = 2. * np.pi * span
        phi = np.linspace(-phi_max, phi_max, N)
        phi += self.phi_ext
        dphi = np.median(np.diff(phi))
        
        Phi1_0 = self.wavefunc(level1, phi, basis='flux')
        Phi1_2pi = self.wavefunc(level1, phi - 2. * np.pi, basis='flux')
        Phi2_0 = self.wavefunc(level2, phi, basis='flux')
        Phi2_2pi = self.wavefunc(level2, phi - 2. * np.pi, basis='flux')
        
        return np.abs(np.sum(Phi1_0 * Phi1_2pi - Phi2_0 * Phi2_2pi) * dphi)
