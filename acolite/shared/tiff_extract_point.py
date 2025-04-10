## def tiff_extract_point
## extracts data from tiff file for a given lon, lat position
## useable but work in progress
## written by Quinten Vanhellemont, RBINS
## 2025-03-31
## modifications: 2025-04-09 (QV) added rhor

def tiff_extract_point(file, st_lon, st_lat, box_size = 1, shift_edge = False):
    import os, re
    import numpy as np

    from osgeo import gdal
    gdal.UseExceptions()

    import acolite as ac

    if not os.path.exists(file):
        print('File {} does not exist'.format(file))
        return

    ## get dataset names
    datasets = []
    with gdal.Open(file) as ds:
        for bi in range(ds.RasterCount):
            bds = ds.GetRasterBand(bi+1)
            dataset_name = bds.GetDescription()
            datasets.append(dataset_name)

    ## read file projection and determine lat/lon
    dct = ac.shared.projection_read(file)
    lon, lat = ac.shared.projection_geo(dct, add_half_pixel = False)

    ## find pixel
    tmp = ((lon - st_lon)**2 + (lat - st_lat)**2)**0.5
    i, j = np.where(tmp == np.nanmin(tmp))
    i = i[0]
    j = j[0]

    box_size_pix = box_size * 1

    ## determine subset
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
                return

        if (i0 + box_size_pix) > lat.shape[0]-1:
                if shift_edge:
                    i0 = lat.shape[0]-1 - box_size_pix
                    print('Point at the edge of scene, setting i0 to {} for extracting box'.format(i0))
                else:
                    print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                    return

        j0 = j - hbox
        if j0 < 0:
            if shift_edge:
                j0 = 0
                print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                return

        if (j0 + box_size_pix) > lat.shape[1]-1:
            if shift_edge:
                j0 = lat.shape[1]-1 - box_size_pix
                print('Point at the edge of scene, setting j0 to {} for extracting box'.format(j0))
            else:
                print('Point at the edge of scene, cannot extract {}x{} box'.format(box_size_pix, box_size_pix))
                return

        ##  set up subset
        gsub = [j0, i0, box_size_pix, box_size_pix]

    ## read tiff data
    meta, stack = ac.shared.read_band(file, gdal_meta = True)
    if len(datasets) == stack.shape[0]:
        dataset_names = [ds for ds in datasets]
    else:
        dataset_names = ['band_{}'.format(i+1) for i in range(stack.shape[0])]

    #print(lon[gsub[1]:gsub[1]+gsub[3], gsub[0]:gsub[0]+gsub[2]])
    #print(lat[gsub[1]:gsub[1]+gsub[3], gsub[0]:gsub[0]+gsub[2]])

    dct = {}
    dct['data']  = {ds: stack[di, gsub[1]:gsub[1]+gsub[3], gsub[0]:gsub[0]+gsub[2]] for di, ds in enumerate(dataset_names)}
    del stack
    dct['data']['lon'] = lon[gsub[1]:gsub[1]+gsub[3], gsub[0]:gsub[0]+gsub[2]]
    dct['data']['lat'] = lat[gsub[1]:gsub[1]+gsub[3], gsub[0]:gsub[0]+gsub[2]]
    del lon, lat

    dct['datasets'] = list(dct['data'].keys())
    dct['mean'] = {ds: np.nanmean(dct['data'][ds]) for ds in dct['datasets']}
    dct['std'] = {ds: np.nanstd(dct['data'][ds]) for ds in dct['datasets']}

    dct['rhow_datasets'] = [ds for ds in dct['datasets'] if 'rhow' in ds]
    #dct['rhow_wave'] = [int(re.findall(r'\b\d+\b', ds)[0]) for ds in dct['rhow_datasets']]
    dct['rhow_wave'] = [int(re.findall(r'\d+', ds)[0]) for ds in dct['rhow_datasets']]

    dct['rhor_datasets'] = [ds for ds in dct['datasets'] if 'rhor' in ds]
    dct['rhor_wave'] = [int(re.findall(r'\d+', ds)[0]) for ds in dct['rhor_datasets']]

    dct['band_datasets'] = [ds for ds in dct['datasets'] if 'band_' in ds]
    dct['band_idx'] = [int(re.findall(r'\d+', ds)[0]) for ds in dct['band_datasets']]

    return(dct)
