## def dem_lonlat
## returms dem for given lon lat arrays
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-09
## modifications: 2022-07-06 (QV) added Copernicus DEM
##                2022-07-07 (QV) added SRTM1 DEM
##                2026-07-02 (QV) added SRTMGL3S, moved srtm to dem.srtm

def dem_lonlat(lon, lat, source='copernicus30', default='copernicus30'):
    import acolite as ac
    import os
    import numpy as np

    if source not in ['srtm', 'srtm15plus', \
                              'srtmgl1', 'srtmgl3', 'srtmgl3s',
                              'copernicus30', 'copernicus90',
                              'COP-DEM_GLO-30-DGED__2021_1', 'COP-DEM_GLO-30-DGED__2022_1',
                              'COP-DEM_GLO-90-DGED__2021_1', 'COP-DEM_GLO-90-DGED__2022_1']:
        print('DEM {} not recognised, using {}.'.format(source, default))
        source = '{}'.format(default)

    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)

    dem = None
    if source.lower().startswith('srtm'):
        if source.lower() in ['srtm', 'srtmgl1', 'srtmgl3', 'srtmgl3s']:
            dem = ac.dem.srtm.hgt_lonlat(lon, lat, source = source)
        elif source.lower() == 'srtm15plus':
            dem = ac.dem.srtm15plus_lonlat(lon, lat)
        else:
            print('dem_source={} not configured.'.format(source))

    if source in ['copernicus30', 'copernicus90',
                          'COP-DEM_GLO-30-DGED__2021_1', 'COP-DEM_GLO-30-DGED__2022_1',
                          'COP-DEM_GLO-90-DGED__2021_1', 'COP-DEM_GLO-90-DGED__2022_1']:
        dem = ac.dem.copernicus_dem_lonlat(lon, lat, source=source)

    if dem is not None:
        if dem.shape == (1,): dem = dem[0]
    return(dem)
