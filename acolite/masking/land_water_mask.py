## def land_water_mask
## rasterizes a land/water mask optionally including lakes
## written by Quinten Vanhellemont, RBINS
## 2024-02-29

def land_water_mask(ncf, add_lakes = True, extend = False, extend_km = 20, poly_lakes = None, poly_land = None):
    import os, time
    import numpy as np
    from osgeo import ogr,osr,gdal
    import acolite as ac

    ## get lake dataset
    if (add_lakes) & (poly_lakes is None):
        poly_lakes = ac.shared.polylakes('worldlakes')
        if ac.settings['run']['verbosity'] > 3: print('Using lake dataset {}'.format(poly_lakes))

    ## get land mask dataset
    if poly_land is None:
        poly_land = ac.shared.polylakes('gshhg')
        if ac.settings['run']['verbosity'] > 3: print('Using land dataset {}'.format(poly_land))

    ## read datasets and gatts
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## get nc projection
    nc_projection = ac.shared.nc_projection_read(ncf)
    if (nc_projection is None):
        if ac.settings['run']['verbosity'] > 2: print('Could not read nc_projection, not extending grid.')

    ## get pixel size
    if nc_projection is not None: ## from nc projection
        dct = ac.shared.nc_projection_dct(nc_projection)
        pixel_size = 1.0 * dct['pixel_size'][0]

    ## in case the image needs to be extended
    if extend:
        ## amount of pixels to extend the image
        extend_metres = extend_km * 1000
        xext = int(np.round(extend_metres/dct['pixel_size'][0]))
        yext = int(np.round(extend_metres/dct['pixel_size'][1]))

        ## get extended ranges
        xrange_ = [dct['xrange'][0], dct['xrange'][1]]
        yrange_ = [dct['yrange'][0], dct['yrange'][1]]

        ## extend right
        xrange_[1] = dct['xrange'][1] + dct['pixel_size'][0] * xext
        ## extend_left
        xrange_[0] = dct['xrange'][0] - dct['pixel_size'][0] * xext
        ## extend_top
        yrange_[0] = dct['yrange'][0] + dct['pixel_size'][1] * yext
        ## extend_bottom
        yrange_[1] = dct['yrange'][1] - dct['pixel_size'][1] * yext

        ## extended dimensions
        xdim_ = int((xrange_[1]-xrange_[0]+dct['pixel_size'][0])/dct['pixel_size'][0])
        ydim_ = int((yrange_[1]-yrange_[0]+dct['pixel_size'][1])/dct['pixel_size'][1])

        ## update dct
        dct_old = {k:dct[k] for k in dct}
        dct['xrange'] = xrange_
        dct['yrange'] = yrange_
        dct['xdim'] = xdim_
        dct['ydim'] = ydim_
        dct['dimensions'] = (xdim_, ydim_)

        ## subsetting for extended array
        y0 = int(np.abs((dct_old['yrange'][0]-yrange_[0])/dct_old['pixel_size'][1]))
        y1 = int(ydim_ - np.abs((dct_old['yrange'][1]-yrange_[1])/dct_old['pixel_size'][1]))
        x0 = int(np.abs((dct_old['xrange'][0]-xrange_[0])/dct_old['pixel_size'][0]))
        x1 = int(xdim_ - np.abs((dct_old['xrange'][1]-xrange_[1])/dct_old['pixel_size'][0]))

    ## set up target dataset for rasterisation
    xSrc = dct['xdim']
    ySrc = dct['ydim']
    gt = dct['xrange'][0], dct['pixel_size'][0], 0.0,\
         dct['yrange'][0], 0.0, dct['pixel_size'][1]
    if 'Wkt' in dct:
        wkt = dct['Wkt']
    else:
        epsg = dct['epsg']

    srs = osr.SpatialReference()
    if wkt is not None:
        srs.ImportFromWkt(wkt)
    elif epsg is not None:
        srs.ImportFromEPSG(epsg)
        wkt = srs.ExportToWkt()

    ## set partial reprojection needed for GSHHG land mask
    pr = gdal.GetConfigOption('OGR_ENABLE_PARTIAL_REPROJECTION')
    gdal.SetConfigOption('OGR_ENABLE_PARTIAL_REPROJECTION', 'TRUE')

    ## in memory target dataset
    drv = gdal.GetDriverByName('MEM')
    target_ds = drv.Create('', xSrc, ySrc, 1,  gdal.GDT_Byte)
    target_ds.SetGeoTransform(gt)
    target_ds.SetProjection(wkt)

    ## rasterize land water mask
    t0 = time.time()
    opt_land = gdal.RasterizeOptions(allTouched=True, burnValues=[1])
    err = gdal.Rasterize(target_ds, poly_land, options=opt_land)
    t1 = time.time()
    if ac.settings['run']['verbosity'] > 4: print('Rasterizing land mask took {:.1f} seconds'.format(t1-t0))

    ## rasterize lakes
    if add_lakes:
        opt_lakes = gdal.RasterizeOptions(allTouched=True, burnValues=[0]) #
        err = gdal.Rasterize(target_ds, poly_lakes, options=opt_lakes)
        t2 = time.time()
        if ac.settings['run']['verbosity'] > 4: print('Rasterizing lakes took {:.1f} seconds'.format(t2-t1))

    ## set back to original
    gdal.SetConfigOption('OGR_ENABLE_PARTIAL_REPROJECTION', pr)

    ## read dataset
    data = target_ds.ReadAsArray()
    target_ds = None

    if not extend:
        return(data)
    else:
        return(data, (y0,y1, x0,x1), dct)
