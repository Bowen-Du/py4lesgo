import numpy as np


#TODO using np.gradient and scipy.interpolate.RegularGridInterpolator replace
#these functions
def partialx(p, x):
    y = np.zeros_like(x)

    i = 0
    y[i, :, :] = (x[i + 1, :, :] - x[p.nx - 1, :, :]) / 2 / p.dx
    y[1:p.nx - 1, :, :] = (x[2:p.nx, :, :] - x[0:p.nx - 2, :, :]) / 2 / p.dx
    i = p.nx - 1
    y[i, :, :] = (x[0, :, :] - x[i - 1, :, :]) / 2 / p.dx

    return y


def partialy(p, x):
    y = np.zeros_like(x)

    j = 0
    y[:, j, :] = (x[:, j + 1, :] - x[:, p.ny - 1, :]) / 2 / p.dy
    y[:, 1: p.ny - 1, :] = (x[:, 2: p.ny, :] - x[:, 0: p.ny - 2, :]) / 2 / p.dy
    j = p.ny - 1
    y[:, j, :] = (x[:, 0, :] - x[:, j - 1, :]) / 2 / p.dy

    return y


def partialz_uv_w(p, x):
    y = np.zeros_like(x)

    k = 0
    y[:, :, k] = x[:, :, k] * 2 / p.dz  #see it in the MKE.py, need log law
    y[:, :, 1:np.shape(x)[2] - 1] = (x[:, :, 1: np.shape(x)[2] - 1] - x[:, :, 0: np.shape(x)[2] - 2]) / p.dz
    k = np.shape(x)[2] - 1
    y[:, :, k] = 0

    return y


def partialz_w_uv(p, x):
    y = np.zeros_like(x)

    y[:, :, 0:np.shape(x)[2] - 1] = (x[:, :, 1:np.shape(x)[2]] - x[:, :, 0:np.shape(x)[2] - 1]) / p.dz
    k = np.shape(x)[2] - 1
    y[:, :, k] = 0

    return y


def interp_w_uv(p, x):
    y = np.zeros_like(x)
    
    y[:, :, 0:np.shape(x)[2] - 1] = (x[:, :, 1:np.shape(x)[2]] + x[:, :, 0:np.shape(x)[2] - 1]) / 2
    k = np.shape(x)[2] - 1
    y[:, :, k] = x[:, :, k]

    return y


def interp_uv_w(p, x):
    y = np.zeros_like(x)

    k = 0
    y[:, :, k] = x[:, :, k]
    y[:, :, 1:np.shape(x)[2]] = (x[:, :, 1:np.shape(x)[2]] + x[:, :, 0:np.shape(x)[2]-1]) / 2
    
    return y
