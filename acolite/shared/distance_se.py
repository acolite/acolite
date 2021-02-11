## def distance_se
## computes sun - earth distance from day of year
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications:

def distance_se(doy):
    import numpy as np
    doy=float(doy)
    se = 1.00014-0.01671*np.cos(np.pi*(0.9856002831*doy-3.4532868)/180.)-0.00014*np.cos(2*np.pi*(0.9856002831*doy-3.4532868)/180.)
    return(se)
