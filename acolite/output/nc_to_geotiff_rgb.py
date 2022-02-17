## def nc_to_geotiff_rgb
## makes a RGB GeoTIFF from an ACOLITE NetCDF
## written by Quinten Vanhellemont, RBINS
## 2022-01-02
## modifications: 2022-02-17 (QV) fix for paths with spaces

def nc_to_geotiff_rgb(f, settings = {}, remove_temp_files=True, oformat='GTiff'):
    import acolite as ac
    import numpy as np
    from osgeo import gdal
    import os, sys, subprocess

    creationOptions = None
    setu =  ac.acolite.settings.parse(None, settings=settings)

    ## get attributes and datasets
    gatts = ac.shared.nc_gatts(f)
    datasets_file = ac.shared.nc_datasets(f)
    if 'ofile' in gatts:
        out = gatts['ofile'].replace('.nc', '')
    else:
        out = f.replace('.nc', '')

    ## find datasets and wavelengths
    rhot_ds = [ds for ds in datasets_file if ds[0:5] == 'rhot_']
    rhot_wave = [int(ds.split('_')[-1]) for ds in rhot_ds]
    rhos_ds = [ds for ds in datasets_file if ds[0:5] == 'rhos_']
    rhos_wave = [int(ds.split('_')[-1]) for ds in rhos_ds]

    for base in ['rhot', 'rhos']:
        if (base == 'rhot'):
            if (setu['rgb_rhot'] is False): continue
            if len(rhot_ds) < 3: continue
            cwaves = [w for w in rhot_wave]

        if (base == 'rhos'):
            if (setu['rgb_rhot'] is False): continue
            if len(rhos_ds) < 3: continue
            cwaves = [w for w in rhos_wave]

        ## make temporary band files
        tempfiles = []
        for ii, wave in enumerate([setu['rgb_red_wl'], setu['rgb_green_wl'],setu['rgb_blue_wl']]):
            wi, w = ac.shared.closest_idx(cwaves, wave)
            ds = '{}_{}'.format(base, w)
            outfile = '{}_{}{}'.format(out, ds, '_temp.tif')
            ifile = '{}_{}{}'.format(out, ds, '.tif')

            ## input scale
            scale_in = setu['rgb_min'][ii], setu['rgb_max'][ii]
            scale_out = [0, 255]

            # standard scaling
            if not setu['rgb_autoscale']:
                scale_cur = [scale_in[0],scale_in[1]]
            else:
                data = ac.shared.nc_data(f, ds)
                scale_cur = np.nanpercentile(data, setu['rgb_autoscale_percentiles'])
                data = None

            ## Byte scaling options
            options = gdal.TranslateOptions(options='-ot Byte -scale {:.2f} {:.2f} {} {}'.format(scale_cur[0], scale_cur[1],
                                                                                                 scale_out[0], scale_out[1]))
            ## warp
            if os.path.exists(ifile): ## use existing tif file
                dt = gdal.Translate(outfile, ifile, format=oformat,
                                    creationOptions=creationOptions, options=options)
            else:
                dt = gdal.Translate(outfile, 'NETCDF:"{}":{}'.format(f, ds), format=oformat,
                                creationOptions=creationOptions, options=options)
            ## set no data value
            dt.GetRasterBand(1).SetNoDataValue(0)
            tempfiles.append("'{}'".format(outfile))
            dt = None

        ## composite to RGB
        outfile = '{}_{}{}'.format(out, '{}_RGB'.format(base), '.tif')
        if os.path.exists(outfile):
            os.remove(outfile)
        cmd = ['gdal_merge.py',  "-o '{}'".format(outfile), '-separate'] + tempfiles
        ret = os.popen(' '.join(cmd)).read()
        #sp = subprocess.run(' '.join(cmd),shell=True, check=True, stdout=subprocess.PIPE)
        if os.path.exists(outfile):
            print('Wrote {}'.format(outfile))

        ## remove tempfiles
        if remove_temp_files:
            for file in tempfiles:
                os.remove(file.strip("'"))
