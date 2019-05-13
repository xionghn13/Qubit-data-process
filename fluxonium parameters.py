from QubitDataProcessPackages import *
import scipy.constants as ct
phi0 = ct.hbar / 2 / ct.e
EJ = 5 * ct.h * 1e9
LJ = phi0 ** 2 / EJ
print('LJ=%.3GnH' % (LJ * 1e9))
