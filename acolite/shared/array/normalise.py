## def normalise
## normalise array to value according to linear interpolation across given axis
## written by Quinten Vanhellemont, RBINS
## 2026-06-01
## modifications: 2026-06-03 (QV) separated linear weight

def normalise(x, y, x_val, axis = 0):
    import acolite as ac
    idx_0, w0, idx_1, w1 = ac.shared.array.linear_weight(x, x_val)

    if axis == 0:
        scale_y = y[idx_0, :, :] * w0 + y[idx_1, :, :] * w1
    else:
        print('Axis = {} not configured'.format(axis))
    return(y / scale_y)
