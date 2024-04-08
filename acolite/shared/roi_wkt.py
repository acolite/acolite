## def roi_wkt
## converts a given roi (polygon file, 4 element limit, wkt) into (bounding) wkt
## written by Quinten Vanhellemont, RBINS
## 2023-09-20
## modifications: 2024-04-08 (QV) added geojson

def roi_wkt(roi):
    import os
    import acolite as ac
    from osgeo import ogr,osr,gdal

    wkt = None
    limit = None
    if type(roi) is str:
        ## try to read from polygon
        if os.path.exists(roi):
            try:
                wkt = ac.shared.polygon_wkt(roi)
            except:
                pass

        ## try if roi is valid wkt
        if wkt is None:
            try:
                geom = ogr.CreateGeometryFromWkt(roi)
                if geom is not None: wkt  = '{}'.format(roi)
            except:
                pass

        ## try if roi is valid geojson
        if wkt is None:
            try:
                geom = ogr.CreateGeometryFromJson(roi)
                if geom is not None: wkt  =  geom.ExportToWkt()
            except:
                pass

        ## maybe its a limit string
        if wkt is None:
            if len(roi.split(',')) == 4: limit = [float(r.strip()) for r in roi.split(',')]

    ## try if roi is limit list
    if wkt is None:
        if type(roi) is list:
            if len(roi) == 4: limit = [r for r in roi]
        if limit is not None:
            try:
                wkt = ac.shared.limit_wkt(limit)
            except:
                pass
    return(wkt)
