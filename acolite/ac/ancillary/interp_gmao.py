## interp_gmao
## interpolates GMAO GEOS data from hourly files to given lon, lat and time (float)
##
## written by Quinten Vanhellemont, RBINS
## 2023-10-16
## modifications:

def interp_gmao(files, lon, lat, isodate, datasets = ['PS', 'TO3', 'TQV', 'U10M', 'V10M'], method = 'linear'):
    import acolite as ac
    import os, sys
    import numpy as np
    import datetime, dateutil.parser
    from scipy import interpolate

    ## requested date/time
    dt = dateutil.parser.parse(isodate)
    ftime = dt.hour + dt.minute/60 + dt.second/3600

    ## input geolocation dimensions
    dim = np.atleast_1d(lon).shape
    onedim = ((len(dim) == 1) & (dim[0] == 1))

    ## run through files
    interp_data = {ds:[] for ds in datasets}
    ftimes = []
    jdates = []
    for file in files:
        #ds = ac.shared.nc_datasets(file)

        ## read lon/lat
        lons = ac.shared.nc_data(file, 'lon')
        lats = ac.shared.nc_data(file, 'lat')

        ## make lons/lats 2D for reproject2
        if not onedim:
            shape = lats.shape[0], lons.shape[0]
            lons = np.repeat(np.broadcast_to(lons, (1, shape[1])), shape[0], axis=0)
            lats = np.repeat(np.expand_dims(lats, axis=1), shape[1], axis=1)

        ## get modelled time
        gatts = ac.shared.nc_gatts(file)
        file_dt = dateutil.parser.parse(gatts['time_coverage_start'])
        file_ftime = file_dt.hour + file_dt.minute/60 + file_dt.second/3600

        ftimes.append(file_ftime)
        jdates.append(int(file_dt.strftime("%j")))

        for dataset in datasets:
            data = ac.shared.nc_data(file, dataset)
            ## interpolation in space
            if onedim:
                #if method == 'nearest':
                #    xi,xret = min(enumerate(lons), key=lambda x: abs(x[1]-float(lon)))
                #    yi,yret = min(enumerate(lats), key=lambda x: abs(x[1]-float(lat)))
                #    interp_data[dataset].append(data[yi,xi])
                #else:
                interp = interpolate.RegularGridInterpolator((lats, lons), data, method=method)
                idata = interp((lat, lon))
                interp_data[dataset].append(idata)
            ## add QC?
            else:
                interp_data[dataset].append(ac.shared.reproject2(data, lons, lats, lon, lat,
                                            nearest=kind == 'nearest',  radius_of_influence=10e5))

    ## add check for year for files[-1]?
    if (ftimes[-1] == 0.) & \
        ((jdates[-1] == jdates[0]+1) | (jdates[0] >= 365 & jdates[-1] == 1)): ftimes[-1] = 24.0

    ## do interpolation in time
    anc_data = {}
    if (ftime >= ftimes[0]) & (ftime <= ftimes[-1]):
        ## linear interpolation weigths
        ip = np.interp(ftime, ftimes, np.arange(len(ftimes)))
        i0 = int(np.floor(ip))
        i1 = i0+1
        w0 = 1 - (ip-i0)
        w1 = 1 - w0
        for dataset in datasets:
            ti = (w0 * interp_data[dataset][i0]) + (w1 * interp_data[dataset][i1])
            anc_data[dataset] = {"interp":ti, "series":interp_data[dataset]}
    return(anc_data)
