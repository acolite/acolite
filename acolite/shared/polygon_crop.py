## def polygon_crop
## crops to given polygon file (geojson or shapefile tested, other files supported my ogr may work too)
## written by Quinten Vanhellemont, RBINS
## 2021-02-23
## modifications: 2021-02-23 (QV) renamed from crop_to_polygon

def polygon_crop(source, poly, return_sub = False):
    import os
    import numpy as np
    from osgeo import ogr,osr,gdal

    ## if target is a file copy information from there
    if type(source) is str:
        if os.path.exists(source):
            g = gdal.Open(source)
            xSrc = g.RasterXSize
            ySrc = g.RasterYSize
            gt = g.GetGeoTransform()
            pr = osr.SpatialReference(wkt=g.GetProjection()).ExportToProj4()
            g = None
    elif type(source) is dict:
        xSrc = source['xdim']
        ySrc = source['ydim']
        gt = source['xrange'][0], source['pixel_size'][0], 0.0,\
             source['yrange'][0], 0.0, source['pixel_size'][1]
        pr = source['proj4_string']

    srs = osr.SpatialReference()
    srs.ImportFromProj4(pr)
    wkt = srs.ExportToWkt()

    ## in memory source dataset based chosen band
    drv = gdal.GetDriverByName('MEM')
    target_ds = drv.Create('', xSrc, ySrc, 1,  gdal.GDT_Byte)
    target_ds.SetGeoTransform(gt)
    target_ds.SetProjection(wkt)

    ## open vector layer
    vector_ds = ogr.Open(poly)
    lyr = vector_ds.GetLayer()

    err = gdal.RasterizeLayer(target_ds, [1], lyr, options=['ALL_TOUCHED=True'])
    data = target_ds.ReadAsArray()

    target_ds = None
    vector_ds = None

    if not return_sub:
        return(data)
    else:
        s = np.where(data)
        if len(s[0]) == 0:
            print('Polygon not in target dataset.')
            return(None, None)
        xmin, xmax = np.min(s[0]), np.max(s[0])
        ymin, ymax = np.min(s[1]), np.max(s[1])
        sub = [int(ymin), int(xmin), int(ymax-ymin), int(xmax-xmin)]
        mask = data[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]].astype(bool)
        mask = mask == False
        return(sub, mask)
