## def read_nat
## get TOA data from nat file
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-11
## modifications:

def read_nat(file, b, sub = None, radiance = False, return_meta = False):
    import os
    from osgeo import gdal
    gdal.UseExceptions()

    if not os.path.exists(file): return

    ## normal resolution bands
    if (b > 0) & (b < 12):
        if radiance:
            ds = gdal.Open('RAD:{}'.format(file))
        else:
            ds = gdal.Open('{}'.format(file))
        band_ds = ds.GetRasterBand(b)
    ## HRV
    elif b == 12:
        ds = gdal.Open('HRV:{}'.format(file))
        band_ds = ds.GetRasterBand(1)
    else:
        print('Band {} not configured for SEVIRI data {}'.format(b, file))
        return

    ## dimensions
    nrows = ds.RasterYSize
    ncols = ds.RasterXSize
    transform = ds.GetGeoTransform()
    meta = ds.GetMetadata_Dict()

    if sub is None:
        data = band_ds.ReadAsArray()
    else:
        data = band_ds.ReadAsArray(sub[0], sub[1], sub[2], sub[3])

    ## calibrate to radiance
    if (not radiance) & (b != 12):
        offset, slope = [float(v) for v in meta['ch{}_cal'.format('{}'.format(b).zfill(2))].split()]
        data = offset + slope * data

    ## close datasets
    band_ds = None
    ds = None

    if return_meta:
        return(data, meta)
    else:
        return(data)
