## interp_gmao
## interpolates GMAO GEOS data from hourly files to given lon, lat and time (float)
##
## written by Quinten Vanhellemont, RBINS
## 2023-10-16
## modifications: 2025-10-06 (QV) do all interpolations with rgi

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
    #onedim = ((len(dim) == 1) & (dim[0] == 1))

    ## run through files
    interp_data = {ds:[] for ds in datasets}
    ftimes = []
    jdates = []
    for file in files:
        ## read lon/lat
        lons = ac.shared.nc_data(file, 'lon')
        lats = ac.shared.nc_data(file, 'lat')

        ## get modelled time
        gatts = ac.shared.nc_gatts(file)
        file_dt = dateutil.parser.parse(gatts['time_coverage_start'])
        file_ftime = file_dt.hour + file_dt.minute/60 + file_dt.second/3600

        ftimes.append(file_ftime)
        jdates.append(int(file_dt.strftime("%j")))

        ## read data and interpolate
        for dataset in datasets:
            data = ac.shared.nc_data(file, dataset)
            interp = interpolate.RegularGridInterpolator((lats, lons), data, method = method)
            idata = interp((lat, lon))
            interp_data[dataset].append(idata)

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
