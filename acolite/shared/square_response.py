## def square_response
## compute square rsr for center wave and fwhm
## written by Quinten Vanhellemont, RBINS
## 2024-07-03
## modifications:

def square_response(center, fwhm, step = 1, factor = 0.5):
    import numpy as np
    wrange = (center - factor*fwhm, center + factor*fwhm)
    x = np.linspace(wrange[0], wrange[1], int(1+(wrange[1]-wrange[0])/step))
    y = np.zeros(len(x)) + 1
    return(x,y)
