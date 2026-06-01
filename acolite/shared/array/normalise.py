## def normalise
## normalise array to value according to linear interpolation across given axis
## written by Quinten Vanhellemont, RBINS
## 2026-06-01
## modifications:

def linear_weight(x, x_val):
    import numpy as np
    idx = np.interp(x_val, x, np.arange(len(x)))
    idx_0 = int(np.floor(idx))
    w = 1 - (idx - idx_0)
    return(idx_0, w, idx_0+1, (1-w))

def normalise(x, y, x_val, axis = 0):
    import numpy as np

    #idx = np.interp(y_val, y, np.arange(len(y)))
    #idx_0 = int(np.floor(idx))
    #w = 1 - (idx - idx_0)

    idx_0, w0, idx_1, w1 = linear_weight(x, x_val)

    if axis == 0:
        #scale_x = x[idx_0, :, :] * w + (1-w) * x[idx_0+1, :, :]
        scale_y = y[idx_0, :, :] * w0 + y[idx_1, :, :] * w1
    else:
        print('Axis = {} not configured'.format(axis))
    return(y / scale_y)
