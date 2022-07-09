## def polygon_limit
## gets [S, E, N, W] extent from polygon file (geojson or shapefile tested, other files supported my ogr may work too)
## written by Quinten Vanhellemont, RBINS
## 2021-02-23
## modifications:
##                2021-10-21 (shundt@usgs.gov): Use GetExtent method for envelope. This works for single and multi-polygon
##                2022-07-09 (QV) convert srs to WGS84 in degrees

def polygon_limit(poly):
    from osgeo import ogr,osr,gdal
    vector_ds = ogr.Open(poly)
    lyr = vector_ds.GetLayer()
    env = lyr.GetExtent()

    ## convert to WGS84 in degrees
    source_srs = lyr.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)
    transform = osr.CoordinateTransformation(source_srs, target_srs)
    point1 = ogr.CreateGeometryFromWkt("POINT ({} {})".format(env[0], env[2]))
    point1.Transform(transform)
    point2 = ogr.CreateGeometryFromWkt("POINT ({} {})".format(env[1], env[3]))
    point2.Transform(transform)
    limit = [point1.GetX(),point1.GetY(),point2.GetX(),point2.GetY()]

    ## close dataset
    vector_ds = None
    return(limit)
