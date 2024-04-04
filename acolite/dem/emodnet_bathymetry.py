## def emodnet_bathymetry
## returns EMODNET bathymetry for given lon, lat grid
##
## function written by Quinten Vanhellemont, RBINS
## 2024-04-04
## modifications:

def emodnet_bathymetry(lon1, lat1, output = None, filename = 'bathymetry', override = True, dataset = 'elevation', nearest=False):
    import os
    import acolite as ac
    import numpy as np

    ## emodnet data
    lat_range = [15.000520833333333, 89.99947916660017]
    lon_range = [-35.99947916666667, 42.99947916663591]
    url_base = 'https://erddap.emodnet.eu/erddap/griddap/bathymetry_2022.nc?elevation%5B({}):1:({})%5D%5B({}):1:({})%5D'

    ## get limit based on lat lon and find and download NC
    limit = np.nanmin(lat1), np.nanmin(lon1), np.nanmax(lat1), np.nanmax(lon1)

    ## test limits
    if limit[0] < lat_range[0]:
        print('Southern limit {} out of dataset range {}'.format(limit[0], lat_range))
        return
    if limit[2] > lat_range[1]:
        print('Northern limit {} out of dataset range {}'.format(limit[2], lat_range))
        return
    if limit[1] < lon_range[0]:
        print('Western limit {} out of dataset range {}'.format(limit[1], lon_range))
        return
    if limit[3] > lon_range[1]:
        print('Eastern limit {} out of dataset range {}'.format(limit[3], lon_range))
        return

    url = url_base.format(limit[0], limit[2], limit[1], limit[3])

    if output is None: output = ac.settings['run']['output']
    file = '{}/{}.nc'.format(output, filename)
    if (not os.path.exists(file)) | (override):
        print('Downloading bathymetry dataset to {}'.format(file))
        ac.shared.download_file(url, file)

    ## read data
    print('Reading bathymetry from {}'.format(dataset))
    bath, att = ac.shared.nc_read(file, dataset)
    lat0, lata = ac.shared.nc_read(file, 'latitude')
    lon0, lona = ac.shared.nc_read(file, 'longitude')

    ## make 2D lon lat arrays
    lon00 = np.tile(lon0, len(lat0)).reshape(bath.shape)
    lat00 = np.rot90(np.tile(lat0, len(lon0)).reshape(bath.shape[1], bath.shape[0]), k=3)

    ## reproject
    result = ac.shared.reproject2(bath, lon00, lat00, lon1, lat1, nearest=nearest)
    result[result.mask] = np.nan ## hard mask

    return(result)
