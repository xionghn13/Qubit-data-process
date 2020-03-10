import numpy as np


def ComplexToIQ(data_set):
    return np.transpose([np.real(data_set), np.imag(data_set)])


def IQtoComplex(data_set):
    return np.squeeze(data_set[:, 0] + 1j * data_set[:, 1])


def FindAngle(data_set):
    CovarianceMatrix = np.cov(np.real(data_set), np.imag(data_set))
    Eigenvalues, EigenVectors = np.linalg.eig(CovarianceMatrix)
    ìndex = list(Eigenvalues).index(min(Eigenvalues))
    CorrectEigenVector = EigenVectors[ìndex]
    return np.angle(CorrectEigenVector[0] + 1j * CorrectEigenVector[1])


def rotate(data_set, method='covariance'):
    if method == 'covariance':
        theta = FindAngle(data_set)
        return data_set * np.exp(1j * theta)
    else:
        complex_out = data_set
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
       