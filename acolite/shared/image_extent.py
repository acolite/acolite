## def image_extent
## gets geographic extent of image using gdal
## written by Quinten Vanhellemont, RBINS
## 2022-08-14
## modifications: 2022-09-21 (QV) added support to pass open gdal dataset
##                2022-10-05 (QV) added tr_opt

def image_extent(image_file, tr_opt = ['METHOD=RPC', 'RPC_PIXEL_ERROR_THRESHOLD=0.05']):
    from osgeo import gdal

    close = False
    if type(image_file) == str:
        ds = gdal.Open(image_file)
        close = True
    elif type(image_file) == gdal.Dataset:
        ds = image_file
    else:
        print('{} not recognised'.format(image_file))
        return()

    xdim = ds.RasterXSize
    ydim = ds.RasterYSize
    tr = gdal.Transformer(ds, None, tr_opt)

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
    if close: ds = None ## close if we open file here

    return(limit, (xdim, ydim), corners)
