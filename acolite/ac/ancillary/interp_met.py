## ancillary_interp_met
## interpolates NCEP MET data from 6 hourly files to given lon, lat and time (float)
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2017-10-18 (QV) fixed latitude indexing
##                           (QV) added bz2 support
##                2017-10-24 (QV) added option to use nearest neighbour (kind from scipy= ‘linear’, ‘cubic’, ‘quintic’)
##                2018-03-05 (QV) fixed end of year rollover
##                2018-03-12 (QV) added file closing to enable file deletion for Windows
##                2021-03-01 (QV) simplified for acg renamed from ancillary_interp_met
##                2022-11-17 (QV) added 2D interpolation
##                2025-01-28 (QV) switched to LinearNDInterpolator

def interp_met(files, lon, lat, time, datasets=['z_wind','m_wind','press','rel_hum','p_water'], kind='linear'):
    import acolite as ac
    import os, sys, bz2
    from pyhdf.SD import SD, SDC
    import numpy as np
    from scipy import interpolate

    ## input geolocation dimensions
    dim = np.atleast_1d(lon).shape
    onedim = ((len(dim) == 1) & (dim[0] == 1))

    interp_data = {ds:[] for ds in datasets}
    ftimes = []
    jdates = []
    for file in files:
        zipped = False
        # uncompress bz2 files
        if file[-4:len(file)] == '.bz2':
            try:
                zipped=True
                file = file.strip('.bz2')
                file_zipped = file + '.bz2'
                with bz2.open(file_zipped, 'rb') as f: data = f.read()
                with open(file,'wb') as f: f.write(data)
            except:
                print("Error extracting file {}, probably incomplete download".format(file_zipped))
                continue

        f = SD(file, SDC.READ)
        datasets_dic = f.datasets()
        meta = f.attributes()

        ftime = meta['Start Millisec'] / 3600000.
        ftimes.append(ftime)
        jdates.append(meta['Start Day'])

        ## make lons and lats for this file
        lons = np.linspace(meta["Westernmost Longitude"], meta["Easternmost Longitude"],
                        num = meta['Number of Columns'])
        lats = np.linspace(meta["Northernmost Latitude"], meta["Southernmost Latitude"],
                        num = meta['Number of Rows'])

        ## make lons/lats 2D for reproject2
        if not onedim:
            shape = lats.shape[0], lons.shape[0]
            lons = np.repeat(np.broadcast_to(lons, (1, shape[1])), shape[0], axis=0)
            lats = np.repeat(np.expand_dims(lats, axis=1), shape[1], axis=1)

        for dataset in datasets:
            sds_obj = f.select(dataset)
            data = sds_obj.get()
            ## old 1D interp
            if onedim:
                ## do interpolation in space
                if kind == 'nearest':
                    xi,xret = min(enumerate(lons), key=lambda x: abs(x[1]-float(lon)))
                    yi,yret = min(enumerate(lats), key=lambda x: abs(x[1]-float(lat)))
                    interp_data[dataset].append(data[yi,xi])
                else:
                    interp = interpolate.LinearNDInterpolator(lons, lats, data)
                    idata = interp(lon, lat)
                    interp_data[dataset].append(idata[0])
                ## add QC?
            else:
                interp_data[dataset].append(ac.shared.reproject2(data, lons, lats, lon, lat,
                                                   nearest=kind == 'nearest',
                                                   radius_of_influence=10e5))
        f.end()
        f = None

        ### delete unzipped file
        if (zipped): os.remove(file)

    ## add check for year for files[-1]?
    if (ftimes[-1] == 0.) & \
        ((jdates[-1] == jdates[0]+1) | (jdates[0] >= 365 & jdates[-1] == 1)): ftimes[-1] = 24.0

    ## do interpolation in time
    anc_data = {}
    if (time >= ftimes[0]) & (time <= ftimes[-1]):
        ## linear interpolation weigths
        ip = np.interp(time, ftimes, np.arange(len(ftimes)))
        i0 = int(np.floor(ip))
        i1 = i0+1
        w0 = 1 - (ip-i0)
        w1 = 1 - w0
        for dataset in datasets:
            ## old 1D interp
            #tinp = interpolate.interp1d(ftimes,interp_data[dataset])
            #ti = tinp(time).flatten()[0]
            ti = (w0 * interp_data[dataset][i0]) + (w1 * interp_data[dataset][i1])
            anc_data[dataset] = {"interp":ti, "series":interp_data[dataset]}

    return(anc_data)
