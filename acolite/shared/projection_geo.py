## def projection_geo
## gets geolocation for generic image/projection dict
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications:

def projection_geo(dct, xy=False):
    import numpy as np
    xdim = np.linspace(dct['xrange'][0],dct['xrange'][1],dct['xdim']).reshape(1,dct['xdim'])
    ydim = np.linspace(dct['yrange'][0],dct['yrange'][1],dct['ydim']).reshape(dct['ydim'],1)
    xdim = np.tile(xdim, (dct['ydim'],1))
    ydim = np.tile(ydim, (1,dct['xdim']))

    if xy:
        return(xdim, ydim)
    else:
        lon,lat = dct['p'](xdim,ydim,inverse=True)
        return(lon, lat)
    return(sub_dct)
