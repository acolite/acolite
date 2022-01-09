## def srtm15plus_lonlat
## returns srtm15plus data for given lon, lat grid
##
## SRTM15+ is provided by Tozer et al. 2019 https://doi.org/10.1029/2019EA000658
## link to dataset: https://topex.ucsd.edu/WWW_html/srtm15_plus.html
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-09
## modifications:

def srtm15plus_lonlat(lon1, lat1, path=None, sea_level=0):

    import os
    from pyproj import Proj
    import pyproj

    import acolite as ac
    import numpy as np
    from pyresample.bilinear import NumpyBilinearResampler
    from pyresample import geometry

    if path is None:
        file = ac.dem.srtm15plus(path=None \
                                    if 'srtm15plus_path' not in ac.config \
                                    else ac.config['srtm15plus_path'])
    else:
        file = '{}'.format(path)

    limit = np.nanmin(lat1), np.nanmin(lon1), np.nanmax(lat1), np.nanmax(lon1)

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

    lon0 = np.tile(lon[sublon], yi).reshape(yi, xi)
    lat0 = np.rot90(np.tile(lat[sublat], xi).reshape(xi, yi))
    lon = None
    lat = None

    ## read z
    zin, zatt = ac.shared.nc_data(file, 'z', attributes=True, sub = sub)
    zin = np.flipud(zin)

    result = ac.shared.reproject2(zin, lon0, lat0, lon1, lat1, nearest=False, radius_of_influence=5000)

    if sea_level is not None:
        result[result<sea_level] = sea_level

    return(result)
