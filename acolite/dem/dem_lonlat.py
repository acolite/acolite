## def dem_lonlat
## returms dem for given lon lat arrays
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-09
## modifications: 2022-07-06 (QV) added Copernicus DEM
##                2022-07-07 (QV) added SRTM1 DEM

def dem_lonlat(lon, lat, source='copernicus30', default='copernicus30'):
    import acolite as ac
    import os

    if source.lower() not in ['srtm', 'srtm15plus', \
                              'srtmgl1', 'srtmgl3',
                              'copernicus30', 'copernicus90']:
        print('DEM {} not recognised, using {}.'.format(source, default))
        source = 'copernicus30'

    if source.lower() == 'srtm':
        dem = ac.dem.hgt_lonlat(lon, lat)
    if source.lower() in ['srtmgl1', 'srtmgl3']:
        dem = ac.dem.hgt_lonlat(lon, lat, source=source)

    if source.lower() == 'srtm15plus':
        dem = ac.dem.srtm15plus_lonlat(lon, lat)
    if source.lower() in ['copernicus30', 'copernicus90']:
        dem = ac.dem.copernicus_dem_lonlat(lon, lat, source=source.lower())

    return(dem)
