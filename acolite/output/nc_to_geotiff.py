## def nc_to_geotiff
## writes exports datasets from ACOLITE NetCDF to GeoTIFF
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-06-06
## modifications:  QV 2018-06-07 fixed S2 offset bug
##                2018-07-18 (QV) changed acolite import name
##                2021-02-05 (QV) adapted for acolite-gen

def nc_to_geotiff(f, skip_geo=True):
    import acolite as ac

    gatts = ac.shared.nc_gatts(f)
    datasets = ac.shared.nc_datasets(f)
    tags = ['xrange', 'yrange', 'pixel_size', 'proj4_string']
    tags += ['oname', 'ofile', 'sensor']

    if all([t in gatts for t in tags]):
        from osgeo import osr, gdal
        out = gatts['ofile'].replace('.nc', '')

        xrange = gatts['xrange']
        yrange = gatts['yrange']
        pixel_size = gatts['pixel_size']

        ## make WKT
        srs = osr.SpatialReference()
        srs.ImportFromProj4(gatts['proj4_string'])
        wkt = srs.ExportToWkt()

        ## make geotransform
        trans = (xrange[0], pixel_size[0], 0.0, \
                 yrange[0], 0.0, pixel_size[1])

        for ds in datasets:
            if (skip_geo) & (ds in ['lat', 'lon', 'x', 'y']): continue

            data = ac.shared.nc_data(f, ds)
            y,x = data.shape

            driver = gdal.GetDriverByName('GTiff')
            outfile = '{}_{}{}'.format(out, ds, '.tif')
            dataset = driver.Create(outfile,x,y,1,gdal.GDT_Float32)

            dataset.SetGeoTransform(trans)
            dataset.SetProjection(wkt)
            dataset.GetRasterBand(1).WriteArray(data)
            dataset.FlushCache()
            print('Wrote {}'.format(outfile))
    else:
        print('File {} not  recognised. Not outputting GeoTIFF files.'.format(f))
