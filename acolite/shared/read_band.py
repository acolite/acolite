## read data from generic band
##
## written by Quinten Vanhellemont, RBINS
##
## 2020-01-25
## modifications:  2021-02-08 (QV) renamed from generic read, integrated in acolite-gen
##                                 added in col and row diff check from landsat reader

def read_band(file, idx = None, warp_to=None, warp_alg = 'near', # 'cubic', 'bilinear'
                 target_res=None, sub=None, gdal_meta = False):

    import os, sys, fnmatch
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(file)
    nrows=ds.RasterYSize
    ncols=ds.RasterXSize

    if gdal_meta:
        try:
            md = ds.GetMetadata_Dict()
        except:
            md = {}

    if idx is not None:
        ds = None
        tmp = gdal.Open(file)
        ds = tmp.GetRasterBand(idx)

    if warp_to is None:
        if sub is None:
            data = ds.ReadAsArray()
        else:
            cdiff = ncols - (sub[0]+sub[2])
            if cdiff < 0:
                sub[2]+=cdiff
            rdiff = nrows-(sub[1]+sub[3])
            if rdiff < 0:
                sub[3]+=rdiff
            data = ds.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None
    else:
        if len(warp_to) >= 2:
            dstSRS = warp_to[0] ## target projection

            ## target bounds in projected space
            if len(warp_to[1]) == 5:
                outputBounds = warp_to[1][0:4]
                outputBoundsSRS = warp_to[1][4]
            else:
                outputBounds = warp_to[1]
                outputBoundsSRS = dstSRS

            #targetAlignedPixels = True
            targetAlignedPixels = False

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
            if len(warp_to) >= 4:
                xRes = warp_to[2]
                yRes = warp_to[3]
            if (xRes is None) or (yRes is None): targetAlignedPixels = False

            ## use given warp algorithm
            if len(warp_to) >= 5:
                warp_alg = warp_to[4]

            ## warp in memory and read dataset to array
            ## https://gdal.org/python/osgeo.gdal-module.html
            ds = gdal.Warp('', file,
                            xRes = xRes, yRes = yRes,
                            outputBounds = outputBounds, outputBoundsSRS = outputBoundsSRS,
                            dstSRS=dstSRS, targetAlignedPixels = targetAlignedPixels,
                            format='VRT', resampleAlg=warp_alg)
            if idx is not None:
                tmp = ds.GetRasterBand(idx)
                data = tmp.ReadAsArray()
                tmp = None
            else:
                data = ds.ReadAsArray()
            ds = None

    if gdal_meta:
        return(md, data)
    else:
        return(data)
