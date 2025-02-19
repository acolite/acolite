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
##                2024-02-27 (QV) changed writing of nodata, changed COG options
##                2024-03-14 (QV) update settings handling
##                                removed some keywords
##                2024-04-16 (QV) use gem NetCDF handling
##                2025-02-19 (QV) use settings.merge

def nc_to_geotiff(f, settings = None, datasets = None, use_projection_key = True):
    import acolite as ac
    import numpy as np
    import os
    from osgeo import osr, gdal

    ## combine default and user defined settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

    ## get settings from ac.settings
    match_file = setu['export_geotiff_match_file']
    skip_geo = setu['export_geotiff_coordinates'] is False

    ## set creation options for COG
    creationOptions = None
    format = 'GTiff'
    if setu['export_cloud_optimized_geotiff']:
        format = 'COG'
        creationOptions = setu['export_cloud_optimized_geotiff_options']

    ## track outputfiles
    outfiles = []

    ## open file
    gem = ac.gem.gem(f)

    ## get metadata and datasets
    if 'ofile' in gem.gatts:
        out = gem.gatts['ofile'].replace('.nc', '')
    else:
        out = f.replace('.nc', '')

    ## use projection key in NetCDF metadata
    if ('projection_key' in gem.gatts) & (use_projection_key):
        for ds in gem.datasets:
            if datasets is not None:
                if ds not in datasets: continue
            if ds in ['x', 'y', gem.gatts['projection_key']]: continue
            if (skip_geo) & (ds in ['lat', 'lon']): continue
            outfile = '{}_{}{}'.format(out, ds, '.tif')
            ## write geotiff
            dt = gdal.Translate(outfile, 'NETCDF:"{}":{}'.format(f, ds),
                                noData = np.nan,
                                format = format, creationOptions = creationOptions)
            dt = None
            print('Wrote {}'.format(outfile))
            outfiles.append(outfile)
    else:
        tags = ['xrange', 'yrange', 'pixel_size', 'proj4_string']
        if all([t in gem.gatts for t in tags]) or (match_file is not None):
            if match_file is None:
                xrange = gem.gatts['xrange']
                yrange = gem.gatts['yrange']
                pixel_size = gem.gatts['pixel_size']

                ## make WKT
                srs = osr.SpatialReference()
                srs.ImportFromProj4(gem.gatts['proj4_string'])
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

            for ds in gem.datasets:
                if datasets is not None:
                    if ds not in datasets: continue
                if (skip_geo) & (ds in ['lat', 'lon', 'x', 'y']): continue
                if ('projection_key' in gem.gatts):
                    if (ds == gem.gatts['projection_key']): continue

                data = gem.data(ds)
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
                    if (dimx != x) | (dimy != y):
                        print('Warning: dataset shape ({}x{}) does not correspond to match_file shape ({}x{})'.format(x,y,dimx,dimy))
                        dataset = driver.Create(outfile, x, y, 1, dt)
                    else:
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
                outfiles.append(outfile)
        else:
            print('Unprojected data {}. Not outputting GeoTIFF files.'.format(f))

    gem.close()
