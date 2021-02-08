## read data from generic band
##
## written by Quinten Vanhellemont, RBINS
##
## 2020-01-25
## modifications:  2021-02-08 (QV) renamed from generic read, integrated in acolite-gen
##                                 added in col and row diff check from landsat reader

def read_band(file, warp_to=None, warp_alg = 'near', # 'cubic', 'bilinear'
                 target_res=None, sub=None):

    import os, sys, fnmatch
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(file)
    nrows=ds.RasterYSize
    ncols=ds.RasterXSize

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
        if len(warp_to) == 2:
            dstSRS = warp_to[0]
            outputBounds = warp_to[1]

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

            ds = gdal.Warp('', file,
                            xRes = xRes, yRes = yRes,
                            outputBounds = outputBounds, outputBoundsSRS=dstSRS,
                            dstSRS=dstSRS, format='VRT', resampleAlg=warp_alg)
            data = ds.ReadAsArray()
            ds = None

    return(data)
