## def ged_lonlat
## returns ASTER GED data for given lon, lat grid
##
## ASTER Band 10 - 8.3 μm Band 11 - 8.6 μm Band 12 - 9.1 μm Band 13 - 10.6 μm Band 14 - 11.3 μm
## idx0 8.3 μm / idx1 8.6 μm / idx2 9.1 μm / idx3 10.6 μm / idx4 11.3 μm
##
## ASTER GED data is obtained from USGS OpenDAP/hyrax
## https://opendap.cr.usgs.gov/opendap/hyrax/ASTER/ASTT/AG100.003/2000.01.01/contents.html
##
## function written by Quinten Vanhellemont, RBINS
## 2022-08-04
## modifications:

def ged_lonlat(lon1, lat1, bands = [10, 11, 12, 13, 14], nearest = False, fill = True):

    import os
    import acolite as ac
    import numpy as np

    ## get limit based on lat lon and find and download GED tiles
    limit = np.nanmin(lat1), np.nanmin(lon1), np.nanmax(lat1), np.nanmax(lon1)
    ged_tiles = ac.ged.ged_find(limit)

    ##
    bands_all = [10, 11, 12, 13, 14]
    bands = [b for b in bands if b in bands_all]

    ged = None
    ## run through dem files and reproject data to target lat,lon
    for i, ged_tile in enumerate(ged_tiles):
        if len(ged_tile) < 21: continue
        if not os.path.exists(ged_tile):
            print('{} does not exist.'.format(ged_tile))
            continue
        print(ged_tile)
        lon0 = ac.shared.nc_data(ged_tile, 'Geolocation_Longitude')
        lat0 = ac.shared.nc_data(ged_tile, 'Geolocation_Latitude')

        em0 = ac.shared.nc_data(ged_tile, 'Emissivity_Mean')
        em0 = em0.astype(float)/1000
        #gdem0 = ac.shared.nc_data(ged_tile, 'ASTER_GDEM_ASTGDEM')

        cstack = None
        for j in range(em0.shape[0]):
            if bands_all[j] not in bands: continue ## skip bands

            result = ac.shared.reproject2(em0[j,:,:], lon0, lat0, lon1, lat1, nearest=nearest)
            result[result<0] = np.nan
            result[result>1] = 1

            if cstack is None:
                cstack = result * 1.0
            else:
                cstack = np.dstack((cstack, result))

        if ged is None:
            ged = cstack
        else:
            ged[cstack>0] = cstack[cstack>0]

    ## fill holes
    ## nans are dilated a few iterations, then filled with closest values
    ## some holes in the middle of water have low emissivity pixels on their edges
    if fill:
        import scipy.ndimage
        struct = scipy.ndimage.generate_binary_structure(2, 2)
        if ged is not None:
            if len(ged.shape) == 3:
                for i in range(ged.shape[2]):
                    if len(np.where(np.isnan(ged[:,:,i]))[0]) > 0:
                        print('Filling GED holes: B{}'.format(bands[i]))
                        tmp = scipy.ndimage.binary_dilation(np.isnan(ged[:,:,i]), structure=struct, iterations=5)
                        ged[:,:,i][tmp] = np.nan
                        ged[:,:,i] = ac.shared.fillnan(ged[:,:,i])
            else:
                if len(np.where(np.isnan(ged))[0]) > 0:
                    print('Filling GED holes: B{}'.format(bands[0]))
                    tmp = scipy.ndimage.binary_dilation(np.isnan(ged), structure=struct, iterations=5)
                    ged[tmp] = np.nan
                    ged = ac.shared.fillnan(ged)

    return(ged)
