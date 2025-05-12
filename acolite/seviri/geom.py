## def geom
## read geolocation and view geometry for MSG/SEVIRI MSG located at subsatellite point lon_0
## datasets will be computed and saved to NetCDF first runtime
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-18
## modifications: 2025-05-12 (QV) added sensor and ssd keywords, FCI support

def geom(lon_0 = 0.0, sub = None, geolocation = True, geometry = True, sensor = 'SEVIRI', ssd = 0.5):
    import os
    import numpy as np
    import acolite as ac

    if (not geolocation) & (not geometry): return

    ## geolocation and geometry datasets
    if sensor == 'SEVIRI':
        geol = ac.config['data_dir'] + '/GEO/{}/lonlat_{}.nc'.format(sensor, lon_0)
        geom = ac.config['data_dir'] + '/GEO/{}/vaavza_{}.nc'.format(sensor, lon_0)
    else:
        geol = ac.config['data_dir'] + '/GEO/{}/lonlat_{}_{}km.nc'.format(sensor, lon_0, ssd)
        geom = ac.config['data_dir'] + '/GEO/{}/vaavza_{}_{}km.nc'.format(sensor, lon_0, ssd)

    ## compute geolocation
    if (geolocation) & (not os.path.exists(geol)):
        print('Computing {} geolocation for sub centre point {}'.format(sensor, lon_0))
        lon, lat = ac.seviri.lonlat(lon_0 = lon_0, ssd = ssd, sensor = sensor)
        print('Saving geolocation to {}'.format(geol))
        gemo = ac.gem.gem(geol, new = True)
        gemo.write('lon', lon)
        gemo.write('lat', lat)
        gemo.close()
        gemo = None
        print('Saved geolocation to {}'.format(geol))
    ## end compute geolocation

    ## compute geometry
    if (geometry) & (not os.path.exists(geom)):
        ## read geolocation
        print('Reading {} geolocation from {}'.format(sensor, geol))
        gem = ac.gem.gem(geol)
        lon = gem.data('lon')
        lat = gem.data('lat')
        gem.close()
        gem = None
        print('Computing {} geometry for sub centre point {}'.format(sensor, lon_0))
        vaa, vza = ac.seviri.vaavza(lon, lat, lon_0 = lon_0)
        print('Saving geometry to {}'.format(geom))
        gemo = ac.gem.gem(geom, new = True)
        gemo.write('vaa', vaa)
        gemo.write('vza', vza)
        gemo.close()
        gemo = None
        print('Saved geometry to {}'.format(geom))
    ## end compute geolocation

    ## read geolocation dataset
    if geolocation:
        print('Reading {} geolocation from {}'.format(sensor, geol))
        gem = ac.gem.gem(geol)
        lon = gem.data('lon', sub = sub)
        lat = gem.data('lat', sub = sub)
        gem.close()
        gem = None
        if not geometry: return(lon, lat)

    ## read geometry dataset
    if geometry:
        print('Reading {} geometry from {}'.format(sensor, geom))
        gem = ac.gem.gem(geom)
        vaa = gem.data('vaa', sub = sub)
        vza = gem.data('vza', sub = sub)
        gem.close()
        gem = None
        if not geolocation: return(vaa, vza)

    return(lon, lat, vaa, vza)
