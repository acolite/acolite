## def stats
## compute a set of statistics for given x and y
## written by Quinten Vanhellemont, RBINS
## 2026-02-03
## modifications: 2026-02-03 (QV) rewrote old version for integration in acolite.shared

def stats(x_, y_):
    import numpy as np
    import scipy.stats

    ## convert and flatten inputs and check length
    x = np.asarray(x_).flatten()
    y = np.asarray(y_).flatten()
    if len(x) != len(x):
        print('Different number of elements in x ({}) and y ({})'.format(len(x), len(y)))
        return

    ## find finites
    finite = np.where(np.isfinite(x) * np.isfinite(y))
    if len(finite[0]) == 0:
        print('No finite values in both x and y')
        return

    ## subset finites
    x = x[finite]
    y = y[finite]

    ## output dict
    s = {}
    s['n'] = len(x)
    s['x_mean'] = np.mean(x)
    s['y_mean'] = np.mean(y)
    s['x_std'] = np.std(x)
    s['y_std'] = np.std(y)
    s['x_max'] = np.max(x)
    s['y_max'] = np.max(y)
    s['x_min'] = np.min(x)
    s['y_min'] = np.min(y)
    s['x_median'] = np.median(x)
    s['y_median'] = np.median(y)

    ## Pearson correlation
    s['r'] = scipy.stats.pearsonr(x, y)[0]
    s['Rsq'] = s['r'] ** 2

    ## RMA regression
    s['RMA_m'] = np.sign(s['r']) * np.std(y, ddof = 1) / np.std(x, ddof = 1)
    s['RMA_b'] = s['y_mean'] - s['RMA_m'] * s['x_mean']

    ## OLS regression
    m, b, r, p, se = scipy.stats.linregress(x, y)
    s['OLS_m'] = m
    s['OLS_b'] = b

    ## difference
    diff = y - x

    ## mean and median difference
    s['md'] = np.mean(diff)
    s['medd'] = np.median(diff)

    ## mean and median absolute difference
    adiff = np.abs(diff)
    s['mad'] = np.mean(adiff)
    s['medad'] = np.median(adiff)
    del adiff

    ## root mean squared difference
    s['rmsd'] = np.mean(diff**2) ** 0.5
    ## unbiased rmsd RMSD^2 = MD^2 + unb-RMSD^2.
    s['unb-rmsd'] = ((s['rmsd']**2) - (s['mad']**2))**0.5

    ## mean value of x and y
    mn = (y + x)/2

    ## mean relative difference, relative bias
    s['mrd'] = np.mean(diff/mn)

    ydiff = y - s['y_median']
    xdiff = x - s['x_median']
    unb_dev = np.abs(ydiff - xdiff)
    s['mean_unb_dev'] = np.mean(unb_dev)
    s['med_unb_dev'] = np.median(unb_dev)
    del unb_dev

    ## mean unbiased absolute relative difference
    ard = np.abs(diff / mn)
    s['murd'] = np.mean(ard)
    s['med_urd'] = np.median(ard)
    del ard

    ## mean absolute relative difference ref y
    ard = np.abs(diff / y)
    s['mard_y'] = np.nanmean(ard)
    s['med_ard_y'] = np.median(ard)
    del ard

    ## absolute percentage difference ref x
    ard = np.abs(diff / x)
    s['mard_x'] = np.nanmean(ard)
    s['med_ard_x'] = np.median(ard)
    del ard
    del diff, mn

    ## ratio
    ratio = y / x

    ## mean ratio
    s['mr'] = np.mean(ratio)
    s['med_r'] = np.median(ratio)
    del ratio

    return(s)
