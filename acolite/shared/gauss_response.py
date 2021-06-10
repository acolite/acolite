## def gauss_response
## compute gaussian rsr for center wave and fwhm
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-01-30
## modifications: 2021-06-07 (QV) changed numpy import added to generic acolite

def gauss_response(center, fwhm, step=1):
    import numpy as np
    wrange = (center - 1.5*fwhm, center + 1.5*fwhm)
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    x = np.linspace(wrange[0], wrange[1], int(1+(wrange[1]-wrange[0])/step))
    y = np.exp(-((x-center)/sigma)**2 )
    return(x,y)
