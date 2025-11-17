## def zarr.projection
## get projection dict from Sentinel-2 zarr dataset
## written by Quinten Vanhellemont, RBINS
## 2025-11-06
## modifications: 2025-11-12 (QV) made into function

def projection(z, meta = None, s2_target_res = 10):
    import os, zarr, pyproj
    import acolite as ac

    ## get meta
    if type(z) is str: z = zarr.open(z, mode = 'r')
    if meta is None: meta = ac.zarr.meta_parse(z)

    ## get scene projection and extent
    dct = None
    p = pyproj.Proj(meta['horizontal_CRS_code'])
    grid_key = 'r{:.0f}m'.format(s2_target_res)
    if grid_key in z['measurements']['reflectance']:
        dct = {'p': p, 'epsg': p.crs.to_epsg(), 'projection': meta['horizontal_CRS_code']}
        xgrid = z['measurements']['reflectance'][grid_key]['x'][:]
        ygrid = z['measurements']['reflectance'][grid_key]['y'][:]
        dct['dimensions'] = len(xgrid), len(ygrid)
        dct['pixel_size'] = xgrid[1] - xgrid[0], ygrid[1] - ygrid[0]
        dct['xrange'] = [xgrid[0], xgrid[-1]+dct['pixel_size'][0]]
        dct['yrange'] = [ygrid[0], ygrid[-1]+dct['pixel_size'][1]]
        dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0]) # + 1
        dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1]) #+ 1
    ## end set up projection

    return(dct)
