## def warp_and_merge
## warp and merge a list of tiles to common projection (e.g. WV or Pléiades tiles)
## written by Quinten Vanhellemont, RBINS
## 2022-07-22
## modifications: 2022-07-23 (QV) added dem
##                2022-08-13 (QV) added estimate of image extent and resolution

def warp_and_merge(tiles, output = None, limit = None,
                   use_tile_projection = True,
                   delete_warped_tiles = True,
                   rpc_dem = None, find_dem = True,
                   decimals = 0, dct = None, resolution = None, utm = True):
    import os
    import numpy as np
    import acolite as ac

    from osgeo import gdal, gdalconst
    from packaging import version
    if version.parse(gdal.__version__) < version.parse('3.3'):
        from osgeo.utils import gdal_merge
    else:
        from osgeo_utils import gdal_merge
    if output is None: output = '{}/'.format(ac.config['scratch_dir'])
    if type(tiles) is not list: tiles = [tiles]

    ## find tile extent in geographic coordinates
    tile_extents = []
    tile_resolutions = []
    for tile in tiles:
        try:
            tile_limit, (xdim, ydim), corners = ac.shared.image_extent(tile)
            tile_extents.append(tile_limit)
            res = ac.shared.image_resolution((tile_limit, (xdim, ydim), corners), decimals=decimals)
            tile_resolutions.append(res)
        except:
            pass
    image_limit = None
    image_resolution = None
    if len(tile_extents) > 0:
        for ti, t in enumerate(tile_extents):
            if image_limit is None:
                image_limit = t
            else:
                if t[0] < image_limit[0]: image_limit[0] = t[0]
                if t[1] < image_limit[1]: image_limit[1] = t[1]
                if t[2] > image_limit[2]: image_limit[2] = t[2]
                if t[3] > image_limit[3]: image_limit[3] = t[3]
    if len(tile_resolutions) > 0:
        image_resolution = np.mean(tile_resolutions)

    ## find if tiles are projected
    dct_tiles = None
    tile_projections = []
    for tile in tiles:
        try:
            dct_tile = ac.shared.projection_read(tile)
            tile_projections.append(dct_tile)
        except:
            pass

    ## find roi
    ## if no projection is given determine dct from projected tiles if possible
    if len(tile_projections)>0:
        for ri, r in enumerate(tile_projections):
            if ri == 0:
                dct_tiles = {k:r[k] for k in r}
                continue
            for k in r:
                if dct_tiles[k] == r[k]: continue
                if k not in ['xrange', 'yrange', 'xdim', 'ydim', 'dimensions']:
                    print('Crucial key differs: {}'.format(k))
                    continue
                if (k == 'xrange') & (dct_tiles['pixel_size'][0]>0):
                    xr = [min(dct_tiles[k][0], r[k][0]),max(dct_tiles[k][1], r[k][1])]
                    dct_tiles[k] = xr
                elif (k == 'yrange') & (dct_tiles['pixel_size'][1]<0):
                    yr = [max(dct_tiles[k][0], r[k][0]),min(dct_tiles[k][1], r[k][1])]
                    dct_tiles[k] = yr
                else:
                    continue
                    print(k, dct_tiles[k], r[k])
        dct_tiles['xdim'] = int(np.round((dct_tiles['xrange'][1]-dct_tiles['xrange'][0])/dct_tiles['pixel_size'][0]))
        dct_tiles['ydim'] = int(np.round((dct_tiles['yrange'][1]-dct_tiles['yrange'][0])/dct_tiles['pixel_size'][1]))        #stop

    ## set image limit if needed
    if (limit is None) & (image_limit is not None):
        print('Using image limit {}'.format(image_limit))
        limit = image_limit

    ## set image resolution if needed
    if (resolution is None) & (image_resolution is not None):
        print('Using image resolution {}'.format(image_resolution))
        resolution = image_resolution

    ## determine final projection dct
    dct_limit = None
    if (limit is not None):
        if (dct is not None):
            dct_limit = ac.shared.projection_sub(dct, limit)
            print('Used provided output projection with limit {}.'.format(limit))
        elif (dct_tiles is not None) & (use_tile_projection):
            dct_limit = ac.shared.projection_sub(dct_tiles, limit)
            print('Output projection was determined from input files with limit {}.'.format(limit))
        elif (resolution is not None):
            dct_limit, nc_projection, warp_to = ac.shared.projection_setup(limit, resolution, utm=utm)
            print('Computed new UTM output projection at {} m resolution with limit {}.'.format(resolution, limit))
        else:
            print('Could not determine subset for limit.')
    elif dct is not None:
        dct_limit = {k:dct[k] for k in dct}
        print('Used provided output projection.')
    elif (dct_tiles is not None) & (use_tile_projection):
        dct_limit = {dct_tiles[k] for k in dct_tiles}
        print('Output projection was determined from input files.')

    if dct_limit is not None:
        test = [dct_limit['xdim'],dct_limit['ydim']]
        if 'sub' in dct_limit: test += dct_limit['sub']
        if any([s<0 for s in test]):
            print('Negative dimensions found for warping.')
            print(dct_limit)
            return()
        print(dct_limit)
    else:
        print('Output projection will be determined by GDAL.')
        if limit is not None:
            print('User defined limit {} will not be used.'.format(','.join([str(l) for l in limit])))

    ## find DEM
    if (find_dem) & (rpc_dem is None) & (dct_limit is not None):
        pos = dct_limit['p']((dct_limit['xrange'][0],dct_limit['xrange'][0],\
                              dct_limit['xrange'][1],dct_limit['xrange'][1]),\
                             (dct_limit['yrange'][0],dct_limit['yrange'][1],\
                              dct_limit['yrange'][0],dct_limit['yrange'][1]), inverse=True)
        pos_limit = [min(pos[1]), min(pos[0]), max(pos[1]), max(pos[0])]
        dem_files = ac.dem.copernicus_dem_find(pos_limit)

        if len(dem_files) == 1:
            rpc_dem = dem_files[0]
        elif len(dem_files) > 1:
            rpc_dem ='{}/dem_merged.tif'.format(output)
            if os.path.exists(rpc_dem): os.remove(rpc_dem)
            print('Merging {} tiles to {}'.format(len(dem_files), rpc_dem))
            gdal_merge.main(['', '-o', rpc_dem, '-n', '0']+dem_files)

        print(pos_limit)
        print(dem_files)
        print(rpc_dem)

    merge = False
    if len(tiles) > 1: merge = True

    ## warp tiles
    warped_tiles=[]
    for tile in tiles:
        bn, ex = os.path.splitext((os.path.basename(tile)))
        warped_file='{}/{}_warped.tif'.format(output, bn)
        if os.path.exists(warped_file): os.remove(warped_file)
        warped_tile, dim = ac.shared.warp_inputfile(tile, target=warped_file, rpc_dem=rpc_dem, dct=dct_limit)
        print(warped_tile, dim)
        warped_tiles.append(warped_tile)

    ## merge tiles
    if merge:
        merged_file='{}/{}_merged.tif'.format(output, os.path.splitext(os.path.basename(tiles[0]))[0])
        if os.path.exists(merged_file): os.remove(merged_file)
        print('Merging {} tiles to {}'.format(len(warped_tiles), merged_file))
        gdal_merge.main(['', '-o', merged_file, '-n', '0']+warped_tiles)
        dct_merged = ac.shared.projection_read(merged_file)

    if (merge) & (delete_warped_tiles):
        for warped_tile in warped_tiles: os.remove(warped_tile)
        return(merged_file, dct_merged, dct_limit)
    elif (~merge):
        dct_merged = ac.shared.projection_read(warped_tiles[0])
        return(warped_tiles[0], dct_merged, dct_limit)
    else:
        return(warped_tiles, merged_file, dct_merged, dct_limit)
