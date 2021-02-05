## def read_band
## simple image reading for landsat files
## sub keyword (xoff, yoff, xcount, ycount)
##
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: 2018-04-19 (QV) added check for sub size (especially for the L8 pan)

def read_band(file, sub=None):
    import os, sys, fnmatch
    if not os.path.isfile(file):
        print('File '+file+' not found.')
        sys.exit()

    if fnmatch.fnmatch(file,'*.TIF'):
        from osgeo import gdal
        gdal.UseExceptions()
        band = gdal.Open(file)
        nrows=band.RasterYSize
        ncols=band.RasterXSize
        if sub is None:
            data = band.ReadAsArray()
        else:
            cdiff = ncols - (sub[0]+sub[2])
            if cdiff < 0:
                sub[2]+=cdiff

            rdiff = nrows-(sub[1]+sub[3])
            if rdiff < 0:
                sub[3]+=rdiff

            data = band.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None

    return(data)
