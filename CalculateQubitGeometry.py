from QubitDataProcessPackages import *


def twoPointInterp(p1, p2, x3):
    slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
    y3 = slope * (x3 - p2[0]) + p2[1]
    return y3


# EL, real sjj area
q1 = [6.5, 0.37 * 0.12]
q2 = [14.858, 0.59 * 0.11]
EL_target = 13
area_target = twoPointInterp(q1, q2, EL_target)
side1 = 0.11
side2 = area_target / side1
print('area should be %.3G * %.3G' % (side2, side1))

