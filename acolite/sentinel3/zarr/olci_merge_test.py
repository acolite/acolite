## def zarr.olci_merge_test
## finds out whether OLCI zarr images can be merged, and how to subset if the limit is specified
## sub is None if images just have to be merged, 4 elements if the final dataset is a subset of the merge
## also returned are the shape of the merged images, the order in which to merge,
## and the crop position for each of the inputs and the output crop position
## written by Quinten Vanhellemont, RBINS 
## 2025-11-05
## modifications:

def olci_merge_test(bundles, limit = None, use_tpg = True, max_time_diff_sec = 1, max_orbit_diff = 1):
    import acolite as ac
    import numpy as np
    import dateutil.parser
    import zarr

    ## read bundles attributes
    bundles_meta = []
    for bundle in bundles:
        z = zarr.open(bundle, mode='r')
        md = ac.zarr.meta(z)
        meta = ac.zarr.meta_parse(z)
        bundles_meta.append(meta)

    ## find out info about bundles
    orbits = np.asarray([m['sat:absolute_orbit'] for m in bundles_meta])
    start_times = np.asarray([dateutil.parser.parse(m['start_datetime']) for m in bundles_meta])
    end_times = np.asarray([dateutil.parser.parse(m['end_datetime']) for m in bundles_meta])
    products = [m['product:type'] for m in bundles_meta]

    ## test if we can merge
    merge = True
    if len(set(products)) != 1:
        print('Products can not be merged, different product types:')
        print('Products: '+', '.join(products))
        print('Bundles: '+', '.join(bundles))
        merge = False

    ## sort bundles based on start_times
    sort_bundles = np.argsort(start_times)

    ## check orbit numbers
    if np.max(np.diff(orbits[sort_bundles])) > max_orbit_diff:
        print('Products can not be merged, orbit numbers differ by >{}:'.format(max_orbit_diff))
        print('Orbits: '+', '.join(orbits))
        print('Bundles: '+', '.join(bundles))
        merge = False

    ## compute time difference between each bundle
    diff_sec = np.asarray([(end_times[sort_bundles][c]-start_times[sort_bundles][c+1]).total_seconds() for c in range(len(bundles)-1)])
    if any(np.abs(diff_sec)>max_time_diff_sec):
        print('Products can not be merged, time difference greater than {:.1f} sec:'.format(max_time_diff_sec))
        print('Time difference: '+', '.join([str(v) for v in diff_sec]))
        print('Bundles: '+', '.join(bundles))
        merge = False

    ## do merging
    if merge:
        bundles = [bundles[bi] for bi in sort_bundles]
        data_shapes = []
        ## iterate over bundles to construct lat/lon
        for bi, bundle in enumerate(bundles):
            ## open zarr
            z = zarr.open(bundle, mode='r')
            ## get tpg and data shape
            tpg_shape = z['conditions']['geometry']['latitude'].shape
            data_shape = z['conditions']['image']['latitude'].shape
            data_shapes.append(data_shape[0])

            if use_tpg:
                lat = z['conditions']['geometry']['latitude'][:]
                lon = z['conditions']['geometry']['longitude'][:]
            else:
                lat = z['conditions']['image']['latitude'][:]
                lon = z['conditions']['image']['longitude'][:]

            if bi == 0:
                scene_offsets = [0]
                scene_index_merged = np.zeros(data_shape, dtype=int)
                lat_merged = lat * 1.0
                lon_merged = lon * 1.0
                data_shape_merged = [data_shape[0], data_shape[1]]
            else:
                scene_offsets += [data_shapes[-2]]
                scene_index_merged = np.vstack((scene_index_merged, np.zeros(data_shape, dtype=int)+bi))
                lat_merged = np.vstack((lat_merged, lat))
                lon_merged = np.vstack((lon_merged, lon))
                data_shape_merged[0] += data_shape[0]

        ## do subsetting if limit is given
        if limit:
            sub = ac.shared.geolocation_sub(lat_merged, lon_merged, limit)

            if sub is not None:
                ## convert from tpg to per pixel coordinates
                if use_tpg:
                    ys = int((data_shape[1]-1)/(tpg_shape[1]-1))
                    y0s=max(0,sub[0]-1)*ys
                    y1s=min(tpg_shape[1],sub[0]+sub[2]+1)*ys
                    sub[0] = int(y0s)
                    sub[2] = int(y1s-y0s)

                ## limit to merged data extents
                if sub[0]<0: sub[0] = 0
                if sub[1]<0: sub[1] = 0
                if sub[1]+sub[3] > data_shape_merged[0]:
                    sub[3] = data_shape_merged[0]-sub[1]
                if sub[0]+sub[2] > data_shape[1]:
                    sub[2] = data_shape_merged[1]-sub[0]

                ## subset in merged scene
                scene_index_merged_sub = np.zeros(data_shape_merged)*np.nan
                scene_index_merged_sub[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]] = scene_index_merged[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                scene_index_merged = scene_index_merged_sub * 1.0
                scene_index_merged_sub = None
            else:
                print('Limit not in scenes: ', limit)
                print('Bundles: '+', '.join(bundles))
                return
        else:
            sub = None

        ## get the subset, and target in the new array for each scene
        crop_in = []
        crop_out = []
        sort_bundles_out = []
        for bi, bundle in enumerate(bundles):
            ## index range in input bundle
            si = np.where(scene_index_merged == bi)
            if len(si[0]) == 0: continue ## if bundle does not cover limit
            sort_bundles_out.append(sort_bundles[bi])

            cropi = si[1][0], si[1][-1]+1, si[0][0]-scene_offsets[bi], si[0][-1]-scene_offsets[bi]+1 ## for nc_data
            crop_in.append(cropi)

            ## index range in output region
            if sub is None:
                so = np.where(scene_index_merged  == bi)
            else:
                so = np.where(scene_index_merged[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]] == bi)
            so = so[1][0], so[1][-1]+1, so[0][0], so[0][-1]+1,  ## for array subsetting
            crop_out.append(so)

        lat_merged = None
        lon_merged = None
        return(sub, data_shape_merged, sort_bundles_out, crop_in, crop_out)

    else:
        return
