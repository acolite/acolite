## def crop
## finds crop position (x0, y0, ns, nl) for Pleiades image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications: 2025-01-28 (QV) removed 0 indexig after switch to LinearNDInterpolator
##

def crop(metadata, limit):
    import numpy as np
    from acolite.pleiades import geo

    ## get interpolator for row / col
    _, (zcol, zrow) = geo.init(metadata)

    south = limit[0]
    east = limit[1]
    north = limit[2]
    west = limit[3]

    se_x = np.floor(zcol(east,south))
    se_y = np.ceil(zrow(east,south))

    nw_x = np.ceil(zcol(west,north))
    nw_y = np.floor(zrow(west,north))

    ns = nw_x - se_x
    nl = se_y - nw_y

    sub = [se_x, nw_y, ns, nl]
    sub = [int(s) for s in sub]
    return(sub)
