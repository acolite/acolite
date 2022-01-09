## def dem_lonlat
## returms dem for given lon lat arrays
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-09
## modifications:

def dem_lonlat(lon, lat, source='srtm'):
    import acolite as ac
    import os

    if source.lower() not in ['srtm', 'srtm15plus']:
        print('DEM {} not recognised, using srtm.'.format(source))
        source = 'srtm'

    if source.lower() == 'srtm':
        dem = ac.dem.hgt_lonlat(lon, lat)
    if source.lower() == 'srtm15plus':
        dem = ac.dem.srtm15plus_lonlat(lon, lat)

    return(dem)
