## def projection_geo
## gets geolocation for generic image/projection dict
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-11 (QV) added half pixel offset option

def projection_geo(dct, xy=False, add_half_pixel=False):
    import numpy as np

    if not add_half_pixel:
        xdim = np.linspace(dct['xrange'][0],dct['xrange'][1]-dct['pixel_size'][0],dct['xdim']).reshape(1,dct['xdim'])
        ydim = np.linspace(dct['yrange'][0],dct['yrange'][1]-dct['pixel_size'][1],dct['ydim']).reshape(dct['ydim'],1)
    else:
        xdim = np.linspace(dct['xrange'][0],dct['xrange'][1]-dct['pixel_size'][0],dct['xdim']).reshape(1,dct['xdim']) + dct['pixel_size'][0]/2
        ydim = np.linspace(dct['yrange'][0],dct['yrange'][1]-dct['pixel_size'][1],dct['ydim']).reshape(dct['ydim'],1) + dct['pixel_size'][1]/2

    xdim = np.tile(xdim, (dct['ydim'],1))
    ydim = np.tile(ydim, (1,dct['xdim']))

    if xy:
        return(xdim, ydim)
    else:
        lon,lat = dct['p'](xdim,ydim,inverse=True)
        return(lon, lat)
    return(sub_dct)
