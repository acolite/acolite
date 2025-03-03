## def extract_point
## extracts data from an ACOLITE NetCDF file for a given lon, lat position
## written by Quinten Vanhellemont, RBINS
## 2022-03-02
## modifications: 2024-03-21 (QV) added circular extraction, and (p)ixel and (m)etre options for box and radius units
##                2024-03-25 (QV) optional external mask, added output of image subset position
##                2024-04-25 (QV) use gem for file reading
##                2024-06-05 (QV) added lat and lon parameter names as keyword
##                2025-03-03 (QV) check for finite lat/lon

def nc_extract_point(ncf, st_lon, st_lat, extract_datasets = None,
                     box_size = 1, box_size_units = 'p', shift_edge = False,
                     lat_par = 'lat', lon_par = 'lon',
                     extract_circle = False, extract_circle_radius = 1, extract_cicle_units = 'p',
                     external_mask = None):
    import acolite as ac
    import numpy as np

    ## read netcdf attributes and datasets
    gem = ac.gem.gem(ncf)
    gatts = gem.gatts # ac.shared.nc_gatts(ncf)
    datasets = gem.datasets # ac.shared.nc_datasets(ncf)
    for ds in ['transverse_mercator', 'x', 'y']:
        if ds in datasets:
            datasets.remove(ds)

    ## get image resolution
    if 'scene_pixel_size' in gatts:
        resolution = gatts['scene_pixel_size'][0]
    else:
        print('Could not determine resolution from file attributes')

    ## if extracting circle determine radius in pixels
    if extract_circle:
        if extract_cicle_units[0] == 'm':
            radius_pix = extract_circle_radius/resolution
        else:
            radius_pix = 1 * extract_circle_radius

        if (radius_pix <= 1):
            print('Radius size in pixels has to > 1')
            gem.close()
            return
    ## else determine box size in pixels
    else:
        if box_size_units[0] == 'm':
            box_size_pix = box_size/resolution
        else:
            box_size_pix = 1 * box_size
        box_size_pix = np.round(box_size_pix).astype(int)
        if (box_size_pix & 1) == 0:
            print('Box size in pixels has to be odd.')
            gem.close()
            return

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
    lon = gem.data(lon_par) # ac.shared.nc_data(ncf, 'lon')
    lat = gem.data(lat_par) # ac.shared.nc_data(ncf, 'lat')

    if not (np.isfinite(lon).any() & np.isfinite(lat).any()):
        print('No finite lat/lon in scene {}'.format(ncf))
        gem.close()
        return

    ## is requested point in this scene?
    latrange = np.nanpercentile(lat, (0,100))
    lonrange = np.nanpercentile(lon, (0,100))
    if (st_lat < latrange[0]) | (st_lat > latrange[1]) | (st_lon < lonrange[0]) | (st_lon > lonrange[1]):
        print('Point {}N {}E not in scene {}'.format(st_lat, st_lon, ncf))
        gem.close()
        return

    ## find pixel
    tmp = ((lon - st_lon)**2 + (lat - st_lat)**2)**0.5
    i, j = np.where(tmp == np.nanmin(tmp))
    i = i[0]
    j = j[0]

    ## extract data
    if extract_circle:
        ## currently does not test whether circle overlaps with the scene edge
        ## pixel coordinates
        y = np.arange(0, lon.shape[0])
        x = np.arange(0, lon.shape[1])
        ## circular mask
        mask = (x[None, :] - j) ** 2 + (y[:, None] - i) ** 2 < radius_pix**2
        mask_idx = np.where(mask)
        ## subset in pixel grid
        gsub = [np.min(mask_idx[1]), np.min(mask_idx[0]),
                np.max(mask_idx[1]) - np.min(mask_idx[1]) + 1,
                np.max(mask_idx[0]) - np.min(mask_idx[0]) + 1]
        ## mask  in subset
        mask_sub = mask[np.min(mask_idx[0]):np.max(mask_idx[0]) + 1,\
                        np.min(mask_idx[1]):np.max(mask_idx[1]) + 1]
    else:
        if box_size_pix == 1:
            gsub = [j, i, 1, 1]
        else:
            hbox = int(box_size_pix/2)
            i0 = i - hbox
            if i0 < 0:
                if shift_edge:
                    i0 = 0
                    print('Point at the edge of scene, setting i0 to {} for extracting box'.format(i0))
                else:
                    print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                    gem.close()
                    return

            if (i0 + box_size_pix) > lat.shape[0]-1:
                if shift_edge:
                    i0 = lat.shape[0]-1 - box_size_pix
                    print('Point at the edge of scene, setting i0 to {} for extracting box'.format(i0))
                else:
                    print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                    gem.close()
                    return

            j0 = j - hbox
            if j0 < 0:
                if shift_edge:
                    j0 = 0
                    print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
                else:
                    print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                    gem.close()
                    return

            if (j0 + box_size_pix) > lat.shape[1]-1:
                if shift_edge:
                    j0 = lat.shape[1]-1 - box_size_pix
                    print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
                else:
                    print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                    gem.close()
                    return

            ##  set up subset
            gsub = [j0, i0, box_size_pix, box_size_pix]

    ## extract data
    #sub = {ds: ac.shared.nc_data(ncf, ds, sub=gsub) for ds in dataset_list}
    #gem = ac.gem.gem(ncf)
    sub = {ds: gem.data(ds, sub=gsub) for ds in dataset_list}
    gem.close()
    gem = None

    ## get single element for box size 1
    if gsub[2:] == [1,1]: sub = {ds:sub[ds][0,0] for ds in sub}

    ## mask circle edges
    if extract_circle: sub = {ds:sub[ds] * mask_sub for ds in sub}

    ## apply additional external mask
    external_mask_sub = None
    if external_mask is not None:
        if external_mask.shape != lon.shape:
            print('Provided external_mask is the wrong shape.')
            print('Image: {},  external_mask: {}'.format(lon.shape, external_mask.shape))
            return
        else:
            ## subset to crop position
            external_mask_sub = external_mask[gsub[1]:gsub[1]+gsub[3],
                                              gsub[0]:gsub[0]+gsub[2]]
            ## apply mask to extracted data
            sub = {ds:sub[ds] * external_mask_sub for ds in sub}
    ## end external mask

    ## create return dict
    dct = {}
    dct['gatts'] = gatts
    dct['data'] = sub
    dct['sub'] = gsub ## store image subset location
    if extract_circle: dct['mask_sub'] = mask_sub ## store circular mask
    if external_mask_sub is not None: dct['external_mask_sub'] = external_mask_sub ## store external mask

    ## store common spectral datasets and centre wavelengths
    dct['datasets'] = list(dct['data'].keys())
    dct['rhot_datasets'] = [ds for ds in dct['datasets'] if 'rhot_' in ds]
    dct['rhot_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhot_datasets']]
    dct['rhotc_datasets'] = [ds for ds in dct['datasets'] if 'rhotc_' in ds]
    dct['rhotc_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhotc_datasets']]
    dct['rhos_datasets'] = [ds for ds in dct['datasets'] if 'rhos_' in ds]
    dct['rhos_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhos_datasets']]
    dct['rhosu_datasets'] = [ds for ds in dct['datasets'] if 'rhosu_' in ds]
    dct['rhosu_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhosu_datasets']]
    dct['rhoe_datasets'] = [ds for ds in dct['datasets'] if 'rhoe_' in ds]
    dct['rhoe_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhoe_datasets']]
    dct['rhow_datasets'] = [ds for ds in dct['datasets'] if 'rhow_' in ds]
    dct['rhow_wave'] = [int(ds.split('_')[-1]) for ds in dct['rhow_datasets']]
    dct['Rrs_datasets'] = [ds for ds in dct['datasets'] if 'Rrs_' in ds]
    dct['Rrs_wave'] = [int(ds.split('_')[-1]) for ds in dct['Rrs_datasets']]

    ## compute means when more than 1x1 pixel is extracted
    if gsub[2:] != [1,1]:
        dct['mean'] = {ds: np.nanmean(dct['data'][ds]) for ds in dct['data']}
        dct['std'] = {ds: np.nanstd(dct['data'][ds]) for ds in dct['data']}
        dct['median'] = {ds: np.nanmedian(dct['data'][ds]) for ds in dct['data']}
        dct['n'] = {ds: np.count_nonzero(~np.isnan(dct['data'][ds])) for ds in dct['data']}
    return(dct)
