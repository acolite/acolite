## def malh
## (modified) Astoreca Line Height algorithm
## see https://doi.org/10.1093/plankt/fbn116
## and https://doi.org/10.1016/j.rse.2022.113270
##
## algorithm wavelengths in increasing order
## with support for 3D arrays (wavelength, x, y)
##
## written by Quinten Vanhellemont, RBINS
##
## 2026-06-01
## modifications:

def malh(wavelength, reflectance, wave = [470, 482.5, 490, 700], awNIR = 0.57):
    import acolite as ac

    ## get weighting factor
    w = (wave[1] - wave[0])/(wave[2]-wave[0])

    ## extract reflectance interpolated to required wavelengths
    r_ = []
    for wl in wave:
        idx_0, w0, idx_1, w1 = ac.shared.array.linear_weight(wavelength, wl)
        if len(reflectance.shape) == 1:
            r_.append(reflectance[idx_0] * w0 + reflectance[idx_1] * w1)
        elif len(reflectance.shape) == 2:
            r_.append(reflectance[idx_0, :] * w0 + reflectance[idx_1, :] * w1)
        elif len(reflectance.shape) == 3:
            r_.append(reflectance[idx_0, :, :] * w0 + reflectance[idx_1, :, :] * w1)

    ## compute modified astoreca
    ac3 = (1/r_[1] - (1/r_[0]**(1-w)) * (1/r_[2]**w)) * awNIR * r_[3]

    ## delete reference reflectance dataset
    del r_
    
    return(ac3)
