## def read_band
## generic function to read and convert landsat toa data
##
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-09 (QV) added warp_to option

def read_toa(fm, mus = 1, sub=None, warp_to=None, usgs_reflectance = True, usgs_radiance = False, usgs_bt=True):
    import numpy as np
    import acolite as ac

    ## read data
    data = ac.shared.read_band(fm['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)

    ## mask data
    data[data<fm['QUANTIZE_CAL_MIN']] = np.nan
    data[data>fm['QUANTIZE_CAL_MAX']] = np.nan

    ## check if thermal
    thermal = False
    if ('K1_CONSTANT' in fm) & ('K2_CONSTANT' in fm):
        usgs_reflectance = False
        usgs_radiance = True
        thermal = True
        
    if ('REFLECTANCE_MULT' not in fm) or ('REFLECTANCE_ADD' not in fm):
        usgs_reflectance = False

    ## convert to reflectance
    if usgs_reflectance:
        slope = float(fm['REFLECTANCE_MULT'])
        offset = float(fm['REFLECTANCE_ADD'])
        data *= slope
        data += offset
        ## normalise to sun zenith angle
        if len(np.atleast_1d(mus))>1:
            pad = mus.shape[0]-data.shape[0], mus.shape[1]-data.shape[1]
            if pad[0]<0 or pad[1]<0:
                print('read_toa padding error')
                print(pad)
                print(mus.shape)
            data /= mus[0:mus.shape[0]-pad[0], 0:mus.shape[1]-pad[1]]
        else:
            data /= mus
    ## convert to radiance
    else:
        slope = float(fm['RADIANCE_MULT'])
        offset = float(fm['RADIANCE_ADD'])
        data *= slope
        data += offset
        ## or not
        if not usgs_radiance:
            d = fm['se_distance']
            f0 = fm['f0']
            data *= (np.pi * d * d) / (f0 * mus)
        if thermal & usgs_bt:
            data = fm['K2_CONSTANT'] / np.log((fm['K1_CONSTANT']/data)+1.)

    return(data)
