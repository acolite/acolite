## def polygon_limit
## gets [S, E, N, W] extent from polygon file (geojson or shapefile tested, other files supported my ogr may work too)
## written by Quinten Vanhellemont, RBINS
## 2021-02-23
## modifications:
## 2021-10-21 (shundt@usgs.gov): Use GetExtent method for envelope. This works for single and multi-polygon 

def polygon_limit(poly):
    from osgeo import ogr,osr,gdal
    vector_ds = ogr.Open(poly)
    lyr = vector_ds.GetLayer()
    env = lyr.GetExtent()
    vector_ds = None
    limit = [env[2], env[0], env[3], env[1]]
    return(limit)
