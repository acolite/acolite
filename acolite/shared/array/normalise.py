## def normalise
## normalise array to value according to linear interpolation across given axis
## written by Quinten Vanhellemont, RBINS
## 2026-06-01
## modifications: 2026-06-03 (QV) separated linear weight
##                2026-06-09 (QV) use take_along_axis

def normalise(x, y_, x_val, axis = 0):
    import acolite as ac
    import numpy as np

    idx_0, w0, idx_1, w1 = ac.shared.array.linear_weight(x, x_val)

    y = np.asarray(y_)
    dims = y.shape
    axis_ = [i for i in range(len(dims))]

    y0 = np.take_along_axis(y, np.expand_dims(idx_0, axis = axis_), axis = axis)
    y1 = np.take_along_axis(y, np.expand_dims(idx_1, axis = axis_), axis = axis)
    scale_y = y0 * w0 + y1 * w1

    return(y / scale_y)
