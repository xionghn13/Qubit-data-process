import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: height * np.exp(
        -(((center_x - x) / width_x) ** 2 + ((center_y - y) / width_y) ** 2) / 2)

def multi_gaussian(X, Y, param):
    length = len(param)
    n_gaussain = int(length / 5)
    param_mat = np.reshape(param, (n_gaussain, 5))
    print(param_mat.shape)
    Z = 0
    for i in range(n_gaussain):
        print(param_mat[i, :])
        Z += gaussian(*param_mat[i, :])(X, Y)
    return Z

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y


def fitgaussian(X, Y, data):
    """Returns the gaussian parameters of a 2D distribution found by a fit"""
    # params = moments(data)
    # print(np.max(data))
    height_guess = np.max(data)
    width_x_guess = (np.max(X) - np.min(X)) / 5
    width_y_guess = width_x_guess
    center_guess = [[40, 40],
                    [40, 160],
                    [160, 40],
                    [160, 160]]
    n_pts = len(center_guess)
    param_list = []
    for pt in center_guess:
        param_list += [[height_guess] + pt + [width_x_guess, width_y_guess]]
    param_mat = np.array(param_list)
    print(param_mat.shape)
    def error_func(param):
        data_sim = multi_gaussian(X, Y, param)
        return np.ravel(data_sim)
    p, success = optimize.leastsq(error_func, param_mat)
    return p


# Create the gaussian data
Xin, Yin = np.mgrid[0:201, 0:201]
data = gaussian(3, 50, 50, 10, 10)(Xin, Yin) \
       + gaussian(3, 50, 150, 10, 10)(Xin, Yin) \
       + gaussian(3, 150, 50, 10, 10)(Xin, Yin) \
       + gaussian(3, 150, 150, 10, 10)(Xin, Yin) \
       + np.random.random(Xin.shape)

plt.matshow(data, cmap=plt.cm.gist_earth_r)

params = fitgaussian(Xin, Yin, data)
print('-----')
print(params)
fit = multi_gaussian(Xin, Yin, params)

plt.contour(fit, cmap=plt.cm.copper)
ax = plt.gca()
(height, x, y, width_x, width_y) = params[0: 5]

plt.text(0.95, 0.05, """
x : %.1f
y : %.1f
width_x : %.1f
width_y : %.1f""" % (x, y, width_x, width_y),
         fontsize=16, horizontalalignment='right',
         verticalalignment='bottom', transform=ax.transAxes)
plt.show()
