## def warp_from_source
## warps dataset given source to projection dict
## written by Quinten Vanhellemont, RBINS
## 2021-02-23
## modifications: 2022-09-22 (QV) added rpc_dem option
##                2023-07-25 (QV) removed proj4 string

def warp_from_source(source, dct, data, warp_to = None, rpc_dem=None):
    import os
    from osgeo import ogr,osr,gdal

    ## if target is a file copy information from there
    if type(source) is str:
        if os.path.exists(source):
            g = gdal.Open(source)
            xSrc = g.RasterXSize
            ySrc = g.RasterYSize
            gt = g.GetGeoTransform()
            wkt = g.GetProjection()
            RPCs = g.GetMetadata('RPC')
            g = None
        else:
            print('Could not access {}'.format(source))
            return()
    elif type(source) is dict:
        xSrc = source['xdim']
        ySrc = source['ydim']
        gt = source['xrange'][0], source['pixel_size'][0], 0.0,\
             source['yrange'][0], 0.0, source['pixel_size'][1]
        if 'Wkt' in source:
            wkt = source['Wkt']
        else:
            epsg = source['epsg']

    srs = osr.SpatialReference()
    if wkt is not None:
        srs.ImportFromWkt(wkt)
    elif epsg is not None:
        srs.ImportFromEPSG(epsg)
        wkt = srs.ExportToWkt()
    else:
        print('Failed to determine projection.')
        return

    ## in memory source dataset based chosen band
    drv = gdal.GetDriverByName('MEM')
    source_ds = drv.Create('', xSrc, ySrc, 1,  gdal.GDT_Float32)
    source_ds.SetGeoTransform(gt)
    source_ds.SetProjection(wkt)

    ## put data in source_ds
    source_ds.GetRasterBand(1).WriteArray(data)
    source_ds.FlushCache()

    if warp_to is None:
        ## set up the warp
        xyr = [min(dct['xrange']),
               min(dct['yrange'])+dct['pixel_size'][1],
               max(dct['xrange'])+dct['pixel_size'][0],
               max(dct['yrange']),
               dct['proj4_string']]
        warp_to_region = (dct['proj4_string'], xyr,
                          dct['pixel_size'][0], dct['pixel_size'][1],'average')
    else:
        warp_to_region = warp_to

    ## warp the data
    if True:
        dstSRS = warp_to_region[0] ## target projection
        ## target bounds in projected space
        if len(warp_to_region[1]) == 5:
            outputBounds = warp_to_region[1][0:4]
            outputBoundsSRS = warp_to_region[1][4]
        else:
            outputBounds = warp_to_region[1]
            outputBoundsSRS = dstSRS
        #targetAlignedPixels = True
        targetAlignedPixels = False

        target_res = None
        ## if we don't know target resolution, figure out from the outputBounds
        if target_res is None:
            xRes = None
            yRes = None
        else:
            if type(target_res) in (int, float):
                xRes = target_res * 1
                yRes = target_res * 1
            else:
                xRes = target_res[0]
                yRes = target_res[1]

        ## if given use target resolution
        if len(warp_to_region) >= 4:
            xRes = warp_to_region[2]
            yRes = warp_to_region[3]
        #if (xRes is None) or (yRes is None): targetAlignedPixels = False

        ## use given warp algorithm
        if len(warp_to_region) >= 5:
            warp_alg = warp_to_region[4]

        ## add transformeroptions
        rpc = False
        if len(RPCs) > 0: rpc = True
        transformerOptions = []
        if rpc_dem is not None: transformerOptions+=['RPC_DEM={}'.format(rpc_dem)]

        ## warp in memory and read dataset to array
        ## https://gdal.org/python/osgeo.gdal-module.html
        #ds = gdal.Warp('', source_ds,
        #                xRes = xRes, yRes = yRes,
        #                outputBounds = outputBounds, outputBoundsSRS = outputBoundsSRS,
        #                dstSRS=dstSRS, targetAlignedPixels = targetAlignedPixels,
        #                format='VRT', resampleAlg=warp_alg)

        ds = gdal.Warp('', source_ds,
                        xRes = xRes, yRes = yRes,
                        outputBounds = outputBounds, outputBoundsSRS = outputBoundsSRS,
                        dstSRS=dstSRS, targetAlignedPixels = targetAlignedPixels,
                        rpc = rpc, transformerOptions = transformerOptions,
                        format='VRT', resampleAlg=warp_alg)

        data = ds.ReadAsArray()
        ds = None
    source_ds = None
    return(data)
