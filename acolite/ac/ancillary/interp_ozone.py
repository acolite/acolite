## ancillary_interp_ozone
## interpolates NRT ozone data from to given lon, lat
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2017-10-18 (QV) fixed latitude indexing, renamed from ancillary_interp_toast
##                2017-10-24 (QV) added option to use nearest neighbour (kind from scipy= ‘linear’, ‘cubic’, ‘quintic’)
##                2018-03-12 (QV) added file closing
##                2021-03-01 (QV) simplified for acg renamed from ancillary_interp_ozone
##                2022-11-17 (QV) added 2D interpolation
##                2025-01-28 (QV) switched to LinearNDInterpolator

def interp_ozone(file, lon, lat, dataset='ozone', kind='linear'):
    import acolite as ac
    from pyhdf.SD import SD, SDC
    import numpy as np
    from scipy import interpolate

    ## input geolocation dimensions
    dim = np.atleast_1d(lon).shape
    onedim = (len(dim) == 1) & (dim[0] == 1)

    ## open ozone file
    f = SD(file, SDC.READ)
    datasets_dic = f.datasets()
    meta = f.attributes()
    sds_obj = f.select(dataset)
    data = sds_obj.get()
    f.end()
    f = None

    ## make lons and lats for this file
    lons = np.linspace(meta["Westernmost Longitude"], meta["Easternmost Longitude"],
                        num = meta['Number of Columns'])
    lats = np.linspace(meta["Northernmost Latitude"], meta["Southernmost Latitude"],
                        num = meta['Number of Rows'])

    ## make lons/lats 2D for reproject2
    if not onedim:
        lons = np.repeat(np.broadcast_to(lons, (1, data.shape[1])), data.shape[0], axis=0)
        lats = np.repeat(np.expand_dims(lats, axis=1), data.shape[1], axis=1)

    ## old 1D interp
    if onedim:
        ## do interpolation in space
        if kind == 'nearest':
            xi,xret = min(enumerate(lons), key=lambda x: abs(x[1]-float(lon)))
            yi,yret = min(enumerate(lats), key=lambda x: abs(x[1]-float(lat)))
            uoz = data[yi,xi]/1000.
        else:
            interp = interpolate.LinearNDInterpolator(lons, lats, data)
            uoz = (interp(lon, lat))[0]
    ## 2D interp
    else:
        uoz = ac.shared.reproject2(data, lons, lats, lon, lat,
                                   nearest=kind == 'nearest',
                                   radius_of_influence=10e5)

    anc_ozone = {'ozone':{'interp':uoz}}
    return(anc_ozone)
