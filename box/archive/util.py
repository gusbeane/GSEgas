import numpy as np
from numba import njit
import sys

# @njit

# x has shape (N, Nbin) or (Nbin,)
# x0 have shape (N,)
# if x has shape (Nbin,), then x is cast to (N, Nbin)

def _integrate_reshape_arrays(x, y, xmin, xmax):
    assert xmin.shape == xmax.shape
    if len(xmin.shape) == 1:
        N = xmin.shape[0]
    else:
        raise ValueError('xmin must have shape (N,)')

    if len(x.shape) == 1:
        Nbin = x.shape[0]
        x = x.reshape((1, Nbin))
        x = np.repeat(x, N, axis=0)
    elif len(x.shape) == 2:
        Nbin = x.shape[1]
        if x.shape[0] != N:
            raise ValueError('x must have shape (Nbin,) or (N, Nbin)')        
    else:
        raise ValueError('x must have shape (Nbin,) or (N, Nbin)')

    if len(y.shape) == 1:
        assert y.shape[0] == Nbin
        y = y.reshape((1, Nbin))
        y = np.repeat(y, N, axis=0)
    elif len(y.shape) == 2:
        assert y.shape[1] == Nbin
        if y.shape[0] != N:
            raise ValueError('y must have shape (Nbin,) or (N, Nbin)')
    else:
        raise ValueError('y must have shape (Nbin,) or (N, Nbin)')
    
    # double check
    assert x.shape == y.shape

    return x, y, xmin, xmax, N, Nbin

def get_bin_indices(x, x0):
    N = x.shape[0]
    Nbin = x.shape[1]

    i = (x0[:,None] >= x).sum(axis=1) - 1
    i = np.clip(i, 0, Nbin-2)

    x_i = x[np.arange(N),i]
    x_ip1 = x[np.arange(N),i+1]

    s = (x0 - x_i) / (x_ip1 - x_i)
    s = np.clip(s, 0, 1)

    return i, s

def trap_integrate(x, y, xmin, xmax, bin_width):
    x, y, xmin, xmax, N, Nbin = _integrate_reshape_arrays(x, y, xmin, xmax)
    # now x, y have shape (N, Nbin)
    # xmin, xmax have shape (N,)

    # assert bin_width constraint
    assert np.allclose(np.diff(x), bin_width)
    
    result = np.zeros(N)
    
    ilow, slow = get_bin_indices(x, xmin)
    ihigh, shigh = get_bin_indices(x, xmax)
    ihigh += 1

    mask_low = np.arange(Nbin) < ilow[:, np.newaxis]
    mask_high = np.arange(Nbin) > ihigh[:, np.newaxis]

    result = np.copy(y)
    result[mask_low | mask_high] = 0
    result = np.sum(result, axis=1)

    # now fix partial bins. collect y at ilow,ihigh
    y_ilow = y[np.arange(N), ilow]
    y_ilowp1 = y[np.arange(N), ilow+1]
    y_ihigh = y[np.arange(N), ihigh]
    y_ihighm1 = y[np.arange(N), ihigh-1]

    # subtract off partial bin at lower end
    result -= 0.5 * (y_ilow + y_ihigh)

    corr_low = np.zeros_like(result)
    corr_low[slow < 0.5] = slow[slow < 0.5] * y_ilow[slow < 0.5]
    corr_low[slow >= 0.5] = 0.5 * y_ilow[slow >= 0.5] + (slow[slow >= 0.5] - 0.5) * y_ilowp1[slow >= 0.5]
    result -= corr_low

    corr_high = np.zeros_like(result)
    corr_high[shigh < 0.5] = 0.5 * y_ihigh[shigh < 0.5] + (0.5 - shigh[shigh < 0.5]) * y_ihighm1[shigh < 0.5]
    corr_high[shigh >= 0.5] = (1-shigh[shigh >= 0.5]) * y_ihigh[shigh >= 0.5]
    result -= corr_high

    result *= bin_width

    return result

def auto_convert_to_ndarray(func):
    def wrapper(*args, **kwargs):
        # ignore first arg if we are a class method
        arg0 = args[1] if isinstance(args[0], object) else args[0]
        args_in = [np.array([arg]) if isinstance(arg, (float, int)) 
                   else np.asarray(arg) if isinstance(arg, list)
                   else arg for arg in args]
        result = func(*args_in, **kwargs)
        if isinstance(arg0, (float, int)):
            if len(result) == 1:
                result = float(result[0])
            else:
                result = tuple([float(res) for res in result])
        return result
    return wrapper
