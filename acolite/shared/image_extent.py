## def image_extent
## gets geographic extent of image using gdal
## written by Quinten Vanhellemont, RBINS
## 2022-08-14
## modifications:

def image_extent(image_file):
    from osgeo import gdal

    ds = gdal.Open(image_file)
    xdim = ds.RasterXSize
    ydim = ds.RasterYSize
    tr = gdal.Transformer(ds, None, [])

    corners = {'UL':{'pixel':[1,1,0]},
               'UR':{'pixel':[xdim,1,0]},
               'LL':{'pixel':[1,ydim,0]},
               'LR':{'pixel':[xdim,ydim,0]}}
    lons = []
    lats = []
    for corner in corners:
        col, row, z = corners[corner]['pixel']
        s, crd = tr.TransformPoint(0, col, row, z)
        corners[corner]['lon'] = crd[0]
        corners[corner]['lat'] = crd[1]
        lons.append(crd[0])
        lats.append(crd[1])
    limit = [min(lats), min(lons), max(lats), max(lons)]
    ds = None

    return(limit, (xdim, ydim), corners)
