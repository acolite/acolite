## def detector_footprint
## returns raster dataset containing detector footprints based on the GML files
##
## written by Quinten Vanhellemont, RBINS
## 2021-02-17
## modifications:

def detector_footprint(target_file, gml_file):
    from osgeo import ogr,osr,gdal

    ## get target projection from band one
    #g = gdal.Open(safe_files[granule]['B1']['path'])
    g = gdal.Open(target_file)

    ## in memory target dataset based chosen band
    drv = gdal.GetDriverByName('MEM')
    target_ds = drv.Create('', g.RasterXSize, g.RasterYSize, 1,  gdal.GDT_Byte)
    target_ds.SetGeoTransform(g.GetGeoTransform())

    ## set up projections
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    to_proj = osr.SpatialReference()
    to_proj.ImportFromWkt(g.GetProjectionRef())
    target_ds.SetProjection (to_proj.ExportToWkt())
    tx = osr.CoordinateTransformation(wgs84, to_proj)

    # open gml file and select layer
    vector_ds = ogr.Open(gml_file)
    ly = vector_ds.GetLayer(0)
    nfeat = ly.GetFeatureCount()
    dvals = []
    for nf in range(nfeat):
        ## set up destination layer
        drv = ogr.GetDriverByName('Memory')
        dst_ds = drv.CreateDataSource('out')
        dst_layer = dst_ds.CreateLayer('', srs = to_proj, geom_type=ogr.wkbPolygon )
        feat = ly.GetFeature(nf)
        field = feat.GetFieldDefnRef(0)
        dst_layer.CreateField(field)
        defn = dst_layer.GetLayerDefn()

        bv = nf+1
        fname = feat.GetField(0)
        bv = int(fname.split('-')[-2])
        dvals.append(bv)

        # get geometry and transform
        geom = feat.GetGeometryRef()
        geom.Transform(tx)
        # create output geometry
        out_geom = ogr.Feature(defn)
        out_geom.SetGeometry(geom)
        dst_layer.CreateFeature(out_geom)
        dst_layer.ResetReading()
        # delete features
        out_geom.Destroy
        geom.Destroy

        ## burn the vector to the target dataset
        ret = gdal.RasterizeLayer(target_ds, [1], dst_layer,
                    burn_values=[bv],options=['ALL_TOUCHED=True'])
    dvals.sort()

    ## read detector footprints
    data = target_ds.ReadAsArray()
    target_ds = None
    return(dvals, data)
