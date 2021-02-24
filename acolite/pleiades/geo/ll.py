## def ll
## returns lon and lat arrays for Pleiades image (crop supported)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications:
##

def ll(metadata, sub=None):
    import numpy as np
    from acolite.pleiades import geo

    ## get interpolator for lon / lat
    (zlon, zlat), _ = geo.init(metadata)

    if sub is None:
        x0=1
        y0=1
        ns = float(metadata['NCOLS'])
        nl = float(metadata['NROWS'])
    else:
        x0, y0, ns, nl = sub

    x = np.arange(x0, x0+ns, 1)
    y = np.arange(y0, y0+nl, 1)

    lat = zlat(x, y)
    lon = zlon(x, y)

    return lon, lat
