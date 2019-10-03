from QubitDataProcessPackages import *
import scipy.constants as ct
phi0 = ct.hbar / 2 / ct.e
EJ = 5 * ct.h * 1e9
LJ = phi0 ** 2 / EJ
print('LJ=%.3GnH' % (LJ * 1e9))
print('----------------area to EJ----------------')
EJ_last = 14.858
area_last = 0.59 * 0.11
area_now = 0.544 * 0.1267
EJ_now = EJ_last / area_last * area_now
print('EJ now is %.3GGHz.' % EJ_now)

