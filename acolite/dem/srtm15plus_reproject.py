## def srtm15
## get srtm15 z for given dct - too slow for large pixel counts
## written by Quinten Vanhellemont, RBINS
## 2022-01-07
## modifications: 2022-01-09 (QV) renamed to reproject

def srtm15plus_reproject(dct = None, gatts = None, sea_level = 0):
    import os
    from pyproj import Proj
    import pyproj

    import acolite as ac
    import numpy as np
    from pyresample.bilinear import NumpyBilinearResampler
    from pyresample import geometry

    if dct is None:
        if gatts is None:
            print('No dct or gatts given.')
            return()
        dct = {}
        dct['xrange'] = gatts['xrange']
        dct['yrange'] = gatts['yrange']
        dct['pixel_size'] = gatts['pixel_size']
        dct['xdim'] = (dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0]
        dct['ydim'] =(dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1]
        dct['p'] = pyproj.Proj(gatts['proj4_string'])

    file = '/Volumes/SSD/Data/Bathymetry/SRTM15_V2.3.nc'

    #dct = {'xrange': xrange, 'yrange': yrange, 'p': p,
    #       'pixel_size': target_pixel_size, 'xdim': nx, 'ydim': ny}
    #nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)

    xrange = dct['xrange']
    yrange = dct['yrange']
    nx = dct['xdim']
    ny = dct['ydim']
    projection = dct['p'].srs

    l_ = dct['p'](xrange, yrange, inverse=True)
    limit = l_[1][1], l_[0][0], l_[1][0], l_[0][1]

    ## set up target definition
    target_definition = geometry.AreaDefinition('area_id', 'description', 'proj_id',
                                          projection, nx, ny, [xrange[0],yrange[1],xrange[1],yrange[0]])

    ## read lat/lon
    lat = ac.shared.nc_data(file, 'lat')
    lon = ac.shared.nc_data(file, 'lon')

    lonoff = 0.25
    latoff = 0.25
    sublon = np.where((lon >= limit[1]-lonoff) & (lon <= limit[3]+lonoff))
    sublat = np.where((lat >= limit[0]-latoff) & (lat <= limit[2]+latoff))

    sub = [sublon[0][0], sublat[0][0], sublon[0][-1]-sublon[0][0]+1, sublat[0][-1]-sublat[0][0]+1]

    yi = lat[sublat].shape[0]
    xi = lon[sublon].shape[0]

    lons = np.tile(lon[sublon], yi).reshape(yi, xi)
    lats = np.rot90(np.tile(lat[sublat], xi).reshape(xi, yi))

    ## set up source definition
    source_definition = geometry.SwathDefinition(lons=lons, lats=lats)

    ## set up resampler
    resampler = NumpyBilinearResampler(source_definition, target_definition, 30e3)

    zin, zatt = ac.shared.nc_data(file, 'z', attributes=True, sub = sub)
    zin = np.flipud(zin)

    zout = resampler.resample(zin)
    if sea_level is not None:
        zout[zout<=sea_level] = 0

    return(zout)
