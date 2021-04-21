## def intercept
## sorts data and compute intercept using first # pixels
## written by Quinten Vanhellemont, RBINS
## 2021-03-01
## modifications:
##

def intercept(data, pixels):
    import acolite as ac
    import numpy as np
    import scipy.stats

    tmp = data[np.where(np.isfinite(data) & (data >0))]
    idx = np.argsort(tmp)
    npix = min(int(pixels), len(idx))
    pidx = np.arange(npix)
    #my, by, ry, smy, sby = ac.shared.lsqfity(pidx, tmp[idx[0:npix]])
    if len(pidx) == 0: return(np.asarray(0.0))

    slope, intercept, r, p, se = scipy.stats.linregress(pidx, tmp[idx[0:npix]])
    return(np.asarray(intercept))
