## def reproject_acolite_netcdf
## reprojects projected ACOLITE NetCDF data to a defined projection and extent
## written by Quinten Vanhellemont, RBINS
## 2022-01-10
## modifications:

def reproject_acolite_netcdf(ncf, dct, ncfo=None, output=None, settings = {},
                            targetAlignedPixels = False, warp_alg = 'bilinear'):
    import acolite as ac
    import numpy as np
    import os

    from osgeo import ogr,osr,gdal
    from pyproj import Proj

    #try:
    if True:
        ## target projection
        xRes = dct['pixel_size'][0]
        yRes = dct['pixel_size'][1]
        gt_tar = dct['xrange'][0], dct['pixel_size'][0], 0.0, dct['yrange'][0], 0.0, dct['pixel_size'][1]
        pr_tar = dct['proj4_string']
        outputBounds = min(dct['xrange']), min(dct['yrange']), max(dct['xrange']), max(dct['yrange'])
        if 'p' not in dct: dct['p'] = Proj(pr_tar)
        if 'xdim' not in dct:
            dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])
        if 'ydim' not in dct:
            dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])
        nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=False)
        ## end target projection
    else:
        print('No target projection dct given.')
        return()

    try:
        bn = os.path.basename(ncf)
        dn = os.path.dirname(ncf)
        gatts = ac.shared.nc_gatts(ncf)
        datasets = ac.shared.nc_datasets(ncf)
    except:
        print('Could not access {}'.format(ncf))
        return()

    ## parse settings
    setu = ac.acolite.settings.parse(gatts['sensor'], settings=settings)

    ## output name
    odir = dn if output is None else output
    on = bn.replace('.nc', '_reprojected.nc')
    if ncfo is None: ncfo = '{}/{}'.format(odir, on)
    gatts['ofile'] = os.path.splitext(ncfo)[0]

    ## source projection
    keys = ['xrange', 'yrange', 'pixel_size', 'proj4_string']
    dct_src = {k: gatts[k] for k in keys}
    dct_src['xdim'] = int((dct_src['xrange'][1] - dct_src['xrange'][0]) / dct_src['pixel_size'][0])
    dct_src['ydim'] = int((dct_src['yrange'][1] - dct_src['yrange'][0]) / dct_src['pixel_size'][1])
    xSrc = int(dct_src['xdim'])
    ySrc = int(dct_src['ydim'])
    gt_src = dct_src['xrange'][0], dct_src['pixel_size'][0], 0.0, dct_src['yrange'][0], 0.0, dct_src['pixel_size'][1]
    pr_src = dct_src['proj4_string']
    srs_src = osr.SpatialReference()
    srs_src.ImportFromProj4(pr_src)
    wkt_src = srs_src.ExportToWkt()
    ## end source projection

    new = True
    gatts['projection_key'] = [k for k in nc_projection if k not in ['x', 'y']][0]
    for dsname in datasets:

        if True:
            data_in, att = ac.shared.nc_data(ncf, dsname, attributes=True)
            if len(data_in.shape) != 2: continue
            print('Reprojecting {}'.format(dsname))

            ySrc, xSrc = data_in.shape

            ## in memory source dataset based chosen band
            drv = gdal.GetDriverByName('MEM')
            source_ds = drv.Create('', xSrc, ySrc, 1,  gdal.GDT_Float32)
            source_ds.SetGeoTransform(gt_src)
            source_ds.SetProjection(wkt_src)

            ## put data in source_ds
            source_ds.GetRasterBand(1).WriteArray(data_in)
            source_ds.GetRasterBand(1).SetNoDataValue(np.nan)
            source_ds.FlushCache()

            ds = gdal.Warp("", source_ds,
                                xRes = xRes, yRes = yRes,
                                outputBounds = outputBounds, outputBoundsSRS = pr_tar,
                                dstSRS=pr_tar, targetAlignedPixels = targetAlignedPixels,
                                format='VRT', resampleAlg=warp_alg, srcNodata=np.nan, dstNodata=np.nan)
            source_ds = None
        else:
            if dsname in ['x', 'y', 'transverse_mercator', 'polar_stereographic']: continue
            print('Reprojecting {}'.format(dsname))
            att = ac.shared.nc_atts(ncf, dsname)
            ds = gdal.Warp("", 'NETCDF:"{}":{}'.format(ncf, dsname),
                                xRes = xRes, yRes = yRes,
                                outputBounds = outputBounds, outputBoundsSRS = pr_tar,
                                dstSRS=pr_tar, targetAlignedPixels = targetAlignedPixels,
                                format='VRT', resampleAlg=warp_alg, srcNodata=np.nan, dstNodata=np.nan)
        ds.GetRasterBand(1).SetNoDataValue(np.nan)
        data = ds.ReadAsArray()
        data[data>9e36] = np.nan
        ds = None

        print(dsname, data.shape)

        lsd = None
        if dsname not in ['lat', 'lon', 'vza', 'sza', 'vaa', 'saa', 'raa']:
            lsd = setu['netcdf_compression_least_significant_digit']

        ac.output.nc_write(ncfo, dsname, data, attributes = gatts,
                           netcdf_compression=setu['netcdf_compression'],
                           netcdf_compression_level=setu['netcdf_compression_level'],
                           netcdf_compression_least_significant_digit=lsd,
                           nc_projection = nc_projection,
                           dataset_attributes = att, new = new)
        new = False
    return(ncfo)
