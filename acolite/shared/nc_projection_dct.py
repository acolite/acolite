## def nc_projection_dct
## set up projection based on nc_projection
## written by Quinten Vanhellemont, RBINS
## 2024-01-26
## modifications:

def nc_projection_dct(nc_projection):
    from pyproj import Proj
    import numpy as np

    ## x and y data ranges
    xrange = nc_projection['x']['data'][0], nc_projection['x']['data'][-1]
    yrange = nc_projection['y']['data'][0], nc_projection['y']['data'][-1]
    xdim = len(nc_projection['x']['data'])
    ydim = len(nc_projection['y']['data'])

    ## x and y pixel resolution
    xres = (nc_projection['x']['data'][-1]-nc_projection['x']['data'][0])/(len(nc_projection['x']['data'])-1)
    yres = (nc_projection['y']['data'][-1]-nc_projection['y']['data'][0])/(len(nc_projection['y']['data'])-1)

    ## get projection key and wkt
    pk = [k for k in list(nc_projection.keys()) if k not in ['x', 'y']][0]
    wkt = nc_projection[pk]['attributes']['crs_wkt']

    ## set up projection
    p = Proj(wkt)

    ## make acolite generic projection dict
    dct = {'p': p, 'epsg': p.crs.to_epsg(),
                   'Wkt': wkt,  'proj4_string': p.crs.to_proj4(),
                   'xrange': xrange, 'yrange': yrange,
                   'xdim':xdim, 'ydim': ydim,
                   'dimensions':(xdim, ydim),
                   'pixel_size': (xres, yres)}
    dct['projection'] = 'EPSG:{}'.format(dct['epsg'])
    return(dct)
