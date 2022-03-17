## def extract_point
## extracts data from an ACOLITE NetCDF file for a given lon, lat position
## written by Quinten Vanhellemont, RBINS
## 2022-03-02
## modifications:

def nc_extract_point(ncf, st_lon, st_lat, extract_datasets = None, box_size = 1, shift_edge = False):
    import acolite as ac
    import numpy as np

    if (box_size & 1) == 0:
        print('Box size has to be odd.')
        return()

    ## read netcdf attributes and datasets
    gatts = ac.shared.nc_gatts(ncf)
    datasets = ac.shared.nc_datasets(ncf)
    for ds in ['transverse_mercator', 'x', 'y']:
        if ds in datasets:
            datasets.remove(ds)

    ## find datasets to extract
    dataset_list = []
    if extract_datasets == None:
        dataset_list += datasets
    else:
        if type(extract_datasets) is not list:
            extract_datasets = [extract_datasets]
        if ('rhot' in extract_datasets) or ('rhot_*' in extract_datasets):
            dataset_list += [ds for ds in datasets if 'rhot_' in ds]
        if ('rhos' in extract_datasets) or ('rhos_*' in extract_datasets):
            dataset_list += [ds for ds in datasets if 'rhos_' in ds]
        if ('rhow' in extract_datasets) or ('rhow_*' in extract_datasets):
            dataset_list += [ds for ds in datasets if 'rhow_' in ds]
        if ('Rrs' in extract_datasets) or ('Rrs_*' in extract_datasets):
            dataset_list += [ds for ds in datasets if 'Rrs_' in ds]
        for ds in extract_datasets:
            if ds in dataset_list: continue
            if ds not in datasets: continue
            dataset_list.append(ds)
    if len(dataset_list) == 0: return()

    ## read lat lon
    lon = ac.shared.nc_data(ncf, 'lon')
    lat = ac.shared.nc_data(ncf, 'lat')

    ## is requested point in this scene?
    latrange = np.nanpercentile(lat, (0,100))
    lonrange = np.nanpercentile(lon, (0,100))
    if (st_lat < latrange[0]) | (st_lat > latrange[1]) | (st_lon < lonrange[0]) | (st_lon > lonrange[1]):
        print('Point {}N {}E not in scene {}'.format(st_lat, st_lon, ncf))
        return()

    ## find pixel
    tmp = ((lon - st_lon)**2 + (lat - st_lat)**2)**0.5
    i, j = np.where(tmp == np.nanmin(tmp))
    i = i[0]
    j = j[0]

    ## extract data
    if box_size == 1:
        sub = {ds: ac.shared.nc_data(ncf, ds, sub=[j, i, 1, 1])[0,0] for ds in dataset_list}
    else:
        hbox = int(box_size/2)
        i0 = i - hbox
        if i0 < 0:
            if shift_edge:
                i0 = 0
                print('Point at the edge of scene, setting i0 to {} for extracting box'.format(i0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size, box_size))
                return()

        if (i0 + box_size) > lat.shape[0]-1:
            if shift_edge:
                i0 = lat.shape[0]-1 - box_size
                print('Point at the edge of scene, setting i0 to {} for extracting box'.format(i0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size, box_size))
                return()

        j0 = j - hbox
        if j0 < 0:
            if shift_edge:
                j0 = 0
                print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size, box_size))
                return()

        if (j0 + box_size) > lat.shape[1]-1:
            if shift_edge:
                j0 = lat.shape[1]-1 - box_size
                print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size, box_size))
                return()

        ## extract data
        sub = {ds: ac.shared.nc_data(ncf, ds, sub=[j0, i0, box_size, box_size]) for ds in dataset_list}

    ## create return dict
    dct = {}
    dct['gatts'] = gatts
    dct['data'] = sub

    ## store common spectral datasets and centre wavelengths
    dct['datasets'] = list(dct['data'].keys())
    dct['rhot_datasets'] = [ds for ds in dct['datasets'] if 'rhot_' in ds]
    dct['rhot_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhot_datasets']]
    dct['rhos_datasets'] = [ds for ds in dct['datasets'] if 'rhos_' in ds]
    dct['rhos_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhos_datasets']]
    dct['rhow_datasets'] = [ds for ds in dct['datasets'] if 'rhow_' in ds]
    dct['rhow_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhow_datasets']]
    dct['Rrs_datasets'] = [ds for ds in dct['datasets'] if 'Rrs_' in ds]
    dct['Rrs_wave'] = [int(ds.split('_')[-1]) for ds in dct['Rrs_datasets']]

    if box_size > 1:
        dct['mean'] = {ds: np.nanmean(dct['data'][ds]) for ds in dct['data']}
        dct['std'] = {ds: np.nanstd(dct['data'][ds]) for ds in dct['data']}
        dct['median'] = {ds: np.nanmedian(dct['data'][ds]) for ds in dct['data']}
        dct['n'] = {ds: np.count_nonzero(~np.isnan(dct['data'][ds])) for ds in dct['data']}
    return(dct)
