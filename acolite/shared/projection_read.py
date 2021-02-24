## def projection_read
## gets projection dict from target image file
## written by Quinten Vanhellemont, RBINS
## 2021-02-24
## modifications:

def projection_read(file):
    from pyproj import Proj
    from osgeo import gdal,osr

    ## open file
    ds = gdal.Open(file)
    transform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    dimx, dimy = ds.RasterXSize, ds.RasterYSize
    ds = None

    ## get projection info
    src = osr.SpatialReference()
    src.ImportFromWkt(projection)
    proj4_string = src.ExportToProj4()
    p = Proj(proj4_string)

    ## derive projection extent
    x0 = transform[0]
    dx = transform[1]
    y0 = transform[3]
    dy = transform[5]
    pixel_size = (dx, dy)
    xrange = (x0,x0+dimx*dx)
    yrange = (y0,y0+dimy*dy)

    ## make acolite generic dict
    dct = {'p': p, 'epsg': p.crs.to_epsg(),
           'xrange': xrange, 'yrange': yrange,
           'xdim':dimx, 'ydim': dimy,
           'proj4_string':proj4_string,
           'dimensions':(dimx, dimy),
           'pixel_size': pixel_size}
    return(dct)
