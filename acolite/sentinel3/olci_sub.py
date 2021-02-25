## def olci_sub
## finds crop in OLCI image
## written by Quinten Vanhellemont, RBINS 2019-04-02
##              QV 2020-11-28 added more accurate limit subsetting

def olci_sub(bundle, limit, use_tpg=True):
    import acolite as ac
    from numpy import where

    if use_tpg:
        file = '{}/{}.nc'.format(bundle, 'geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        data_shape = lat.shape
        lat = None

        file = '{}/{}.nc'.format(bundle, 'tie_geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        lon = ac.shared.nc_data(file, 'longitude')
        tpg_shape = lat.shape
    else:
        file = '{}/{}.nc'.format(bundle, 'geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        lon = ac.shared.nc_data(file, 'longitude')
        data_shape = lat.shape

    ## new version
    sub = ac.shared.geolocation_sub(lat, lon, limit)
    lat = None
    lon = None

    if sub is not None:
        if use_tpg:
            ys = int((data_shape[1]-1)/(tpg_shape[1]-1))
            y0s=max(0,sub[0]-1)*ys
            y1s=min(tpg_shape[1],sub[0]+sub[2]+1)*ys
            sub[0] = int(y0s)
            sub[2] = int(y1s-y0s)

        ## limit to scene
        if sub[0]<0: sub[0] = 0
        if sub[1]<0: sub[1] = 0
        if sub[1]+sub[3] > data_shape[0]:
            sub[3] = data_shape[0]-sub[1]
        if sub[0]+sub[2] > data_shape[1]:
            sub[2] = data_shape[1]-sub[0]

    return(sub, data_shape)
