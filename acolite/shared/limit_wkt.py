## def limit_wkt
## returns wkt string for ACOLITE limit (S, W, N, E)
## written by Quinten Vanhellemont, RBINS
## 2023-09-12
## modifications:

def limit_wkt(limit):
    from osgeo import ogr,osr,gdal

    ## make new polygon
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint_2D(limit[1], limit[0])
    ring.AddPoint_2D(limit[3], limit[0])
    ring.AddPoint_2D(limit[3], limit[2])
    ring.AddPoint_2D(limit[1], limit[2])
    ring.AddPoint_2D(limit[1], limit[0])

    ## set up new polygon
    ply = ogr.Geometry(ogr.wkbPolygon)
    ply.AddGeometry(ring)
    return(ply.ExportToWkt())
