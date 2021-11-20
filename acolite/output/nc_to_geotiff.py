## def nc_to_geotiff
## writes exports datasets from ACOLITE NetCDF to GeoTIFF
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-06-06
## modifications:  QV 2018-06-07 fixed S2 offset bug
##                2018-07-18 (QV) changed acolite import name
##                2021-02-05 (QV) adapted for acolite-gen
##                2021-02-09 (QV) more generic, only depends on presence of xrange, yrange, pixel_size and p4 string tags
##                2021-02-11 (QV) disabled half pixel offset option, need to check Landsat
##                2021-11-20 (QV) added match_file to extract projection from (esp if data is using RPC for geolocation?)

def nc_to_geotiff(f, skip_geo=True, match_file=None):
    import acolite as ac
    import numpy as np
    import os

    gatts = ac.shared.nc_gatts(f)
    datasets = ac.shared.nc_datasets(f)
    tags = ['xrange', 'yrange', 'pixel_size', 'proj4_string']

    if all([t in gatts for t in tags]) or (match_file is not None):
        from osgeo import osr, gdal
        if 'ofile' in gatts:
            out = gatts['ofile'].replace('.nc', '')
        else:
            out = f.replace('.nc', '')

        if match_file is None:
            xrange = gatts['xrange']
            yrange = gatts['yrange']
            pixel_size = gatts['pixel_size']

            ## make WKT
            srs = osr.SpatialReference()
            srs.ImportFromProj4(gatts['proj4_string'])
            wkt = srs.ExportToWkt()

            ## make geotransform, add half a pixel for pixel centers
            trans = (xrange[0]+pixel_size[0]/2, pixel_size[0], 0.0, \
                     yrange[0]+pixel_size[1]/2, 0.0, pixel_size[1])
        else:
            if os.path.exists(match_file):
                ## get projection info from match file
                src_ds = gdal.Open(match_file)
                transform = src_ds.GetGeoTransform()
                projection = src_ds.GetProjection()
                dimx, dimy = src_ds.RasterXSize, src_ds.RasterYSize
                ## get RPC data
                rpcs = src_ds.GetMetadata('RPC')
                src_ds = None
            else:
                print('File {} not fount. Not outputting GeoTIFF files.'.format(match_file))
                return

        for ds in datasets:
            if (skip_geo) & (ds in ['lat', 'lon', 'x', 'y']): continue

            data = ac.shared.nc_data(f, ds)
            y,x = data.shape
            if data.dtype == np.float32:
                dt = gdal.GDT_Float32
            else:
                print(data.dtype)

            outfile = '{}_{}{}'.format(out, ds, '.tif')

            if match_file is None:
                driver = gdal.GetDriverByName('GTiff')
                dataset = driver.Create(outfile, x, y, 1, dt)
                dataset.SetGeoTransform(trans)
                dataset.SetProjection(wkt)
            else:
                ## create output file
                driver = gdal.GetDriverByName('GTiff')
                dataset = driver.Create(outfile, dimx, dimy, 1, dt)
                dataset.SetGeoTransform(transform)
                dataset.SetProjection(projection)
                #dataset = driver.CreateCopy(outfile, src_ds, 0 )
                ## write RPC data
                dataset.SetMetadata(rpcs ,'RPC')
                src_ds = None

            dataset.GetRasterBand(1).WriteArray(data)
            dataset.GetRasterBand(1).SetNoDataValue(np.nan)
            dataset.FlushCache()
            print('Wrote {}'.format(outfile))
    else:
        print('File {} not  recognised. Not outputting GeoTIFF files.'.format(f))
