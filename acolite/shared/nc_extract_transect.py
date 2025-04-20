## def nc_extract_transect
## extracts transect from an ACOLITE NetCDF file for a line between given lon1, lat1 and lon2, lat2
## written by Quinten Vanhellemont, RBINS
## 2025-04-14
## modifications: 2025-04-20 (QV) added datasets keyword

def nc_extract_transect(file, st_lon1, st_lat1, st_lon2, st_lat2,
                        datasets = None,
                        transect_length = None, pixel_coordinates = False,
                        lat_par = 'lat', lon_par = 'lon', verbosity = 0):
    import os
    import numpy as np
    import acolite as ac

    if not os.path.exists(file):
        print('File does not exist: {}'.format(file))
        return

    ## read netcdf attributes and datasets
    gem = ac.gem.gem(file)
    if gem.gatts is None:
        print('Could not read attributes from file {}'.format(file))
        gem.close()
        return

    ## test datasets
    if datasets is None:
        datasets = gem.datasets
    else:
        for ds in datasets:
            if ds not in gem.datasets:
                datasets.remove(ds)

    for ds in ['transverse_mercator', 'x', 'y']:
        if ds in datasets:
            datasets.remove(ds)

    if len(datasets) == 0:
        print('No datasets in file {}'.format(file))
        gem.close()
        return


    ## read file info
    #datasets = ac.shared.nc_datasets(file)
    #gatts = ac.shared.nc_gatts(file)

    ## read lat lon
    if not pixel_coordinates:
        lat = gem.data(lat_par)
        lon = gem.data(lon_par)

        if not (np.isfinite(lon).any() & np.isfinite(lat).any()):
            print('No finite lat/lon in scene {}'.format(file))
            gem.close()
            return

        ## are requested points in this scene
        latrange = np.nanpercentile(lat, (0,100))
        lonrange = np.nanpercentile(lon, (0,100))
        if (st_lat1 < latrange[0]) | (st_lat1 > latrange[1]) | (st_lon1 < lonrange[0]) | (st_lon1 > lonrange[1]):
            print('Point {}N {}E not in scene {}'.format(st_lat1, st_lon1, file))
            gem.close()
            return
        if (st_lat2 < latrange[0]) | (st_lat2 > latrange[1]) | (st_lon2 < lonrange[0]) | (st_lon2 > lonrange[1]):
            print('Point {}N {}E not in scene {}'.format(st_lat2, st_lon2, file))
            gem.close()
            return

        ## find pixel1
        tmp = ((lon - st_lon1)**2 + (lat - st_lat1)**2)**0.5
        i1, j1 = np.where(tmp == np.nanmin(tmp))
        i1 = i1[0]
        j1 = j1[0]

        ## find pixel2
        tmp = ((lon - st_lon2)**2 + (lat - st_lat2)**2)**0.5
        i2, j2 = np.where(tmp == np.nanmin(tmp))
        i2 = i2[0]
        j2 = j2[0]
        del lat, lon
    else:
        i1 = int(st_lon1)
        j1 = int(st_lat1)
        i2 = int(st_lon2)
        j2 = int(st_lat2)

        xdim = gem.nc.dimensions['x'].size
        ydim = gem.nc.dimensions['y'].size
        if (i1 < 0) | (i2 < 0) | (j1 < 0) | (j2 < 0):
            print('Negative position given for point 1: ({}, {}) or 2 ({}, {})'.format(i1, j1, i2, j2))
            gem.close()
            return
        if (i1 > xdim-1) | (i2 > xdim-1) | (j1 > ydim-1) | (j2 > ydim-1):
            print('Position outside image dimensions ({}, {}) for point 1: ({}, {}) or 2 ({}, {})'.format(xdim, ydim, i1, j1, i2, j2))
            gem.close()
            return

    ## create transect coordinates (NN)
    ## if transect length is None, then extract just the covered pixels
    if transect_length is None: transect_length = int(np.hypot(i2-i1, j2-j1))
    x = np.linspace(i1, i2, transect_length).astype(int)
    y = np.linspace(j1, j2, transect_length).astype(int)
    if verbosity > 1: print('Transect length {} from ({}, {}) to ({}, {})'.format(transect_length, i1, j1, i2, j2))

    ## run through datasets and extract transect
    transect = {}
    for ds in datasets:
        if verbosity > 1: print('Extracting {}'.format(ds))
        transect[ds] = gem.data(ds)[x, y]
    gem.close()

    return(transect)
