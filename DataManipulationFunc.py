import numpy as np


def AutoRotate(complex):
    complex_out = complex
    disp_complex = complex_out - complex_out[0]
    for i in range(3):
        max_disp_ind = np.argmax(np.abs(disp_complex))
        disp_complex = complex_out - complex_out[max_disp_ind]
    max_disp_ind = np.argmax(np.abs(disp_complex))
    max_disp = disp_complex[max_disp_ind]
    normalized_max_disp = max_disp / np.abs(max_disp)
    # print(normalized_max_disp)
    if normalized_max_disp.real < 0:
        normalized_max_disp = -normalized_max_disp
    CM = np.mean(complex_out) * 0  # origin
    complex_out -= CM
    complex_out /= normalized_max_disp
    complex_out += CM
    return complex_out
