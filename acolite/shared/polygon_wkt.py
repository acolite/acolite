## def polygon_wkt
## returns wkt string for polygon file
## written by Quinten Vanhellemont, RBINS
## 2023-09-12
## modifications:

def polygon_wkt(poly):
    from osgeo import ogr,osr,gdal
    import json

    vector_ds = ogr.Open(poly)
    lyr = vector_ds.GetLayer()
    env = lyr.GetExtent()

    ## find features
    nft = lyr.GetFeatureCount()
    wkt = []
    if nft > 1: print('Multiple features found!')
    for i in range(nft):
        ft = lyr.GetFeature(i)
        ## make json version of the feature
        ft_json = ft.ExportToJson()

        ## load points
        points = json.loads(ft_json)['geometry']['coordinates'][0]
        ## make new polygon
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for pt in points:
            ring.AddPoint_2D(pt[0], pt[1])
        ## set up new polygon
        ply = ogr.Geometry(ogr.wkbPolygon)
        ply.AddGeometry(ring)
        wkt.append(ply.ExportToWkt())
    ## close dataset
    vector_ds = None

    if len(wkt) == 1: wkt = wkt[0]
    return(wkt)
