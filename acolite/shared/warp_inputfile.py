## def warp_inputfile
## function to warp inputfile before processing (e.g. Pléiades or Worldview)
## written by Quinten Vanhellemont, RBINS
## 2022-03-30
## modifications: 2022-07-12 (QV) added RPC dem option

def warp_inputfile(file, target=None, dct=None, rpc_dem=None):
    from osgeo import gdal
    import os

    bn = os.path.basename(file)
    bn, ext = os.path.splitext(bn)
    dn = os.path.dirname(file)

    dsi = gdal.Open(file)
    dimxi, dimyi = dsi.RasterXSize, dsi.RasterYSize
    ## get projection info from file
    transform = dsi.GetGeoTransform()
    projection = dsi.GetProjection()
    ## get RPC data
    RPCs = dsi.GetMetadata('RPC')
    ## get GCP data
    GCPs = dsi.GetGCPs()
    GCPProjection = dsi.GetGCPProjection()
    dsi = None

    if False:
        print(dimxi, dimyi)
        print(transform, projection)
        print(RPCs)
        print(GCPs, GCPProjection)

    ## some defaults
    rpc = True if len(RPCs) > 0 else False
    errorThreshold = 0
    targetAlignedPixels = False
    resampleAlg = 'average'

    ## add transformeroptions
    transformerOptions = []
    if rpc_dem is not None: transformerOptions+=['RPC_DEM={}'.format(rpc_dem)]

    ## if target projection is given
    if dct is not None:
        xRes = dct['pixel_size'][0]
        yRes = dct['pixel_size'][1]
        outputBounds = (min(dct['xrange']), min(dct['yrange'])+dct['pixel_size'][1],
                        max(dct['xrange'])+dct['pixel_size'][0],max(dct['yrange']))
        if 'epsg' in dct:
            outputBoundsSRS = dct['epsg']
            dstSRS = dct['epsg']
        else:
            outputBoundsSRS = dct['proj4_string']
            dstSRS = dct['proj4_string']

        wopt = gdal.WarpOptions(rpc=rpc, errorThreshold=errorThreshold,
                                xRes = xRes, yRes = yRes,
                                outputBounds = outputBounds, outputBoundsSRS = outputBoundsSRS,
                                dstSRS=dstSRS, targetAlignedPixels = targetAlignedPixels,
                                resampleAlg=resampleAlg, transformerOptions=transformerOptions)
    else:
        wopt = gdal.WarpOptions(rpc=rpc, errorThreshold=errorThreshold, resampleAlg=resampleAlg,
                                transformerOptions=transformerOptions)

    if target is None:
        ofile = '{}/{}_warped{}'.format(dn, bn, ext)
    else:
        ofile = '{}'.format(target)

    print('Reprojecting {} to {}'.format(file, ofile))
    if os.path.exists(ofile): os.remove(ofile)

    dso = gdal.Warp(ofile, file, options=wopt) ## warp to warp options
    dimxo, dimyo = dso.RasterXSize, dso.RasterYSize
    dso = None

    return((ofile, (dimxo, dimyo)))
