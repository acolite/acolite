## def hillshade_nc
## function to compute hillshade from ACOLITE NetCDF (with isodate and lat/lon) containing DEM data
## written by Quinten Vanhellemont, RBINS
## QV 2021-11-04
## modifications: 2022-07-07 (QV) renamed from nc_dem_hillshade, updated sun position function,
##                                use gdal.Translate with NetCDF for MEM dataset

def hillshade_nc(ncf, write_to_netcdf = False, options = [], # ['-combined']
                     dataset = 'dem', data = None):
    import acolite as ac
    from osgeo import osr, gdal
    import numpy as np

    ## read NetCDF
    datasets = ac.shared.nc_datasets(ncf)
    if data is None:
        if dataset not in datasets:
            print('DEM not in {}'.format(ncf))
            return()
        data = ac.shared.nc_data(ncf, dataset)
    if 'hillshade' in datasets:
        hillshade = ac.shared.nc_data(ncf, 'hillshade')
        return(hillshade)

    ## read NetCDF global attributes
    gatts = ac.shared.nc_gatts(ncf)

    ## compute sun position
    centre_lon = np.nanmedian(ac.shared.nc_data(ncf, 'lon'))
    centre_lat = np.nanmedian(ac.shared.nc_data(ncf, 'lat'))
    sun = ac.shared.sun_position(gatts['isodate'], centre_lon, centre_lat)
    alt = sun['elevation'][0]
    azi = sun['azimuth'][0]

    if False:
        ## get projection info and set up source dataset
        nc_projection = ac.shared.nc_projection_read(ncf)
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

        ## in memory source dataset
        drv = gdal.GetDriverByName('MEM')
        ySrc,xSrc = data.shape
        source_ds = drv.Create('', xSrc, ySrc, 1,  gdal.GDT_Float32)
        source_ds.SetGeoTransform(trans)
        source_ds.SetProjection(wkt)
        ## put data in source_ds
        source_ds.GetRasterBand(1).WriteArray(data)
        source_ds.FlushCache()
        ## remove data
        data = None
    else:
        if 'projection_key' in gatts:
            source_ds = gdal.Translate('', 'NETCDF:"{}":{}'.format(ncf, 'dem'),
                                       format='MEM', creationOptions=None)


    ## compute hillshade
    ds = gdal.DEMProcessing('', source_ds, 'hillshade', options=options, format='MEM',  azimuth=azi, altitude=alt)
    hillshade = ds.ReadAsArray()
    ds = None
    source_ds = None

    if write_to_netcdf:
        ac.output.nc_write(ncf, 'hillshade', hillshade)

    return(hillshade)
