## def tiles_interp
## interpolates tiled dataset to full scene extent
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-11
## modifications: 2020-11-17 (QV) added target mask for limiting interpolation extent and processing time
##                2020-11-18 (QV) added dtype to convert from griddata float64, by default float32
##                                this improves peak memory use when several datasets are kept in memory
##                2021-02-11 (QV) added smooth keyword,  default to nearest
##                2024-04-18 (QV) new version using interpn


def tiles_interp(data, xnew, ynew, smooth = False, kern_size=2, method='nearest', mask=None,
                 target_mask=None, target_mask_full=False, fill_nan = True, dtype='float32'):

    import numpy as np
    from scipy.interpolate import interpn
    from scipy.ndimage import uniform_filter,percentile_filter, distance_transform_edt
    import acolite as ac

    if mask is not None: data[mask] = np.nan

    ## fill nans with closest value
    if fill_nan:
        cur_data = ac.shared.fillnan(data)
    else:
        cur_data = data*1.0

    dim = cur_data.shape
    if smooth: cur_data = uniform_filter(cur_data, size = kern_size)

    ## grid positions
    x = np.arange(0., dim[1], 1)
    y = np.arange(0., dim[0], 1)

    ## interpolate
    if target_mask is None:
        znew = interpn((y,x), cur_data, (ynew[:,None], xnew[None, :]), method = method, bounds_error = False)
    else:
        ## limit to target mask
        vd = np.where(target_mask)
        if target_mask_full:
            znew = np.zeros((len(ynew), len(xnew))).astype(dtype)+np.nan
            znew[vd] = interpn((y,x), cur_data, (ynew[vd[0]], xnew[vd[1]]), method = method, bounds_error = False)
        else:
            znew = interpn((y,x), cur_data, (ynew[vd[0]], xnew[vd[1]]), method = method, bounds_error = False)

    ## to convert data type - scipy always returns float64
    if dtype is not None: znew = znew.astype(np.dtype(dtype))
    return(znew)
