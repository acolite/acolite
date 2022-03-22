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
##                2021-12-08 (QV) added support for the netcdf projection
##                2022-03-22 (QV) added support for match_file with GCP

def nc_to_geotiff(f, skip_geo=True, match_file=None, datasets=None, cloud_optimized_geotiff=False):
    import acolite as ac
    import numpy as np
    import os
    from osgeo import osr, gdal

    creationOptions = None
    format = 'GTiff'
    if cloud_optimized_geotiff:
        creationOptions = ['COMPRESS=DEFLATE', 'PREDICTOR=2', 'OVERVIEWS=NONE', 'BLOCKSIZE=1024']
        format = 'COG'

    gatts = ac.shared.nc_gatts(f)
    datasets_file = ac.shared.nc_datasets(f)
    if 'ofile' in gatts:
        out = gatts['ofile'].replace('.nc', '')
    else:
        out = f.replace('.nc', '')

    if 'projection_key' in gatts:
        for ds in datasets_file:
            if datasets is not None:
                if ds not in datasets: continue
            if ds in ['x', 'y', gatts['projection_key']]: continue
            if (skip_geo) & (ds in ['lat', 'lon']): continue
            outfile = '{}_{}{}'.format(out, ds, '.tif')
            ## write geotiff
            dt = gdal.Translate(outfile, 'NETCDF:"{}":{}'.format(f, ds),
                                format=format, creationOptions=creationOptions)
            ## set no data value
            if True:
                 data = ac.shared.nc_data(f, ds)
                 dt.GetRasterBand(1).WriteArray(data)
                 dt.GetRasterBand(1).SetNoDataValue(np.nan)

            dt = None
            print('Wrote {}'.format(outfile))
    else:
        tags = ['xrange', 'yrange', 'pixel_size', 'proj4_string']
        if all([t in gatts for t in tags]) or (match_file is not None):
            if match_file is None:
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

            else:
                if os.path.exists(match_file):
                    ## get projection info from match file
                    src_ds = gdal.Open(match_file)
                    transform = src_ds.GetGeoTransform()
                    projection = src_ds.GetProjection()
                    dimx, dimy = src_ds.RasterXSize, src_ds.RasterYSize
                    ## get RPC data
                    RPCs = src_ds.GetMetadata('RPC')
                    ## get GCP data
                    GCPs = src_ds.GetGCPs()
                    GCPProjection = src_ds.GetGCPProjection()
                    src_ds = None
                else:
                    print('File {} not found. Not outputting GeoTIFF files.'.format(match_file))
                    return

            for ds in datasets_file:
                if datasets is not None:
                    if ds not in datasets: continue
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
                    driver = gdal.GetDriverByName('GTiff')
                    dataset = driver.Create(outfile, dimx, dimy, 1, dt)
                    ## write RPC data
                    if len(RPCs) > 0:
                        dataset.SetMetadata(RPCs ,'RPC')
                    ## write GCP data
                    if len(GCPs) > 0:
                        dataset.SetGCPs(GCPs, GCPProjection)
                    else:
                        dataset.SetGeoTransform(transform)
                        dataset.SetProjection(projection)
                    src_ds = None

                dataset.GetRasterBand(1).WriteArray(data)
                dataset.GetRasterBand(1).SetNoDataValue(np.nan)
                dataset.FlushCache()
                print('Wrote {}'.format(outfile))
        else:
            print('File {} not recognised. Not outputting GeoTIFF files.'.format(f))
