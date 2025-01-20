## def nc_to_geotiff_rgb
## makes a RGB GeoTIFF from an ACOLITE NetCDF
## written by Quinten Vanhellemont, RBINS
## 2022-01-02
## modifications: 2022-02-17 (QV) fix for paths with spaces
##                2022-12-29 (QV) use output directory if given
##                2022-12-30 (QV) use gdal_merge import to avoid gdal_merge.py not being recognised
##                2024-02-27 (QV) added COG options
##                2024-03-14 (QV) update settings handling
##                                removed some keywords
##                2024-04-16 (QV) use gem NetCDF handling
##                2025-01-20 (QV) added use of rgb_datasets and unified output naming with png outputs

def nc_to_geotiff_rgb(f, settings = None, min_wave = 2,
                        rgb_datasets = ['rhot', 'rhos', 'rhotc', 'rhosu', 'rhorc', 'rhow'],
                        use_gdal_merge_import = True,
                        remove_temp_files = True):

    import os, sys, subprocess
    import acolite as ac
    import numpy as np
    from osgeo import gdal
    gdal.DontUseExceptions()
    from packaging import version
    if version.parse(gdal.__version__) < version.parse('3.3'):
        from osgeo.utils import gdal_merge
    else:
        from osgeo_utils import gdal_merge

    ## combine default and user defined settings
    setu = ac.acolite.settings.parse(None, settings = settings)
    for k in ac.settings['user']: setu[k] = ac.settings['user'][k]

    creationOptions = None
    oformat = 'GTiff'
    if setu['export_cloud_optimized_geotiff']:
        oformat = 'COG'
        creationOptions = setu['export_cloud_optimized_geotiff_options']

    ## open file
    gem = ac.gem.gem(f)
    tags = ['xrange', 'yrange', 'pixel_size', 'proj4_string']
    if ('projection_key' not in gem.gatts) and (~all([t in gem.gatts for t in tags])):
        print('Unprojected data {}. Not outputting GeoTIFF files'.format(f))
        return()

    if 'ofile' in gem.gatts:
        out = gem.gatts['ofile'].replace('.nc', '')
    else:
        out = f.replace('.nc', '')

    ## use output directory if given
    if 'output' in setu:
        if setu['output'] is not None:
            out = '{}/{}'.format(setu['output'], os.path.basename(out))

    ## run through rho datasets
    for base in rgb_datasets:
        key = 'rgb_{}'.format(base)
        if key not in setu: continue
        if setu[key] is False: continue

        cds = [ds for ds in gem.datasets if ds.startswith('{}_'.format(base))]
        if len(cds) < min_wave: continue
        cwaves = [int(ds.split('_')[-1]) for ds in cds]

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
            data = None
            if not setu['rgb_autoscale']:
                scale_cur = [scale_in[0],scale_in[1]]
            else:
                data = gem.data(ds)
                scale_cur = np.nanpercentile(data, setu['rgb_autoscale_percentiles'])
                data = None

            ## Translate options
            ## output format
            options_list = ['-of {}'.format(oformat)]
            options_list += [' -a_nodata 0']
            if creationOptions is not None:
                for co in creationOptions: options_list += ['-co {}'.format(co)]
            options_string = ' '.join(options_list)

            ## Byte scaling
            byte_scaling_string = ' -ot Byte -scale {:.2f} {:.2f} {} {}'.format(scale_cur[0], scale_cur[1], scale_out[0], scale_out[1])
            if oformat == 'COG':
                options_scaling_string = options_string[1] + byte_scaling_string
            else:
                options_scaling_string = options_string + byte_scaling_string

            ## TranslateOptions
            options_creation = gdal.TranslateOptions(options=options_list[0])
            options_scaling = gdal.TranslateOptions(options=options_scaling_string)

            ## warp
            if os.path.exists(ifile): ## use existing tif file
                dt = gdal.Translate(outfile, ifile, options=options_scaling)
            else:
                dt = gdal.Translate(outfile, 'NETCDF:"{}":{}'.format(f, ds), options=options_scaling)

            ## do custom RGB stretch
            if True:
                if data is None: data = gem.data(ds)
                bsc = np.asarray(scale_cur)
                gamma = setu['rgb_gamma'][ii]
                tmp = ac.shared.rgb_stretch(data, gamma = gamma, bsc = bsc, stretch=setu['rgb_stretch'])
                tmp = (tmp*255).astype(np.uint8)
                dt.GetRasterBand(1).WriteArray(tmp)
                dt = None

                ## move to another temporary file to resave as COG
                # temp file does not need to be COG, just convert the RGB to COG
                #if oformat == 'COG':
                #    outfile_temp = '{}_{}{}'.format(out, ds, '_temp2.tif')
                #    os.rename(outfile, outfile_temp)
                #    dt = gdal.Translate(outfile, outfile_temp, options=options_creation)
                #    dt = None
                #    if remove_temp_files: os.remove(outfile_temp)

            tempfiles.append("{}".format(outfile))

        ## composite to RGB
        outfile = '{}_{}{}'.format(out, 'rgb_{}'.format(base), '.tif')
        outfile_temp = '{}_{}{}'.format(out, 'rgb_{}_temp'.format(base), '.tif')
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(outfile_temp): os.remove(outfile_temp)

        ## use gdal_merge to merge RGB bands
        if setu['use_gdal_merge_import']:
            ## run using imported gdal_merge
            ## gives file name errors when running batch processing
            ## so generate default out.tif and rename
            cwd = os.getcwd()
            odir = os.path.dirname(outfile)
            os.chdir(odir)
            gdal_merge.main(['', '-quiet', '-separate'] + tempfiles) ## calling without -o will generate out.tif in cwd
            if os.path.exists('out.tif'):
                if oformat == 'COG': ## COG needs to be converted again
                    os.rename('out.tif', os.path.basename(outfile_temp))
                else:
                    os.rename('out.tif', os.path.basename(outfile))
            os.chdir(cwd)
            ## end imported gdal_merge
        else:
            ## use os to run gdal_merge.py
            tempfiles = ["'{}'".format(tf) for tf in tempfiles] ## add quotes to allow spaces in file path
            if oformat == 'COG': ## COG needs to be converted again
                cmd = ['gdal_merge.py', "-o '{}'".format(outfile_temp), '-separate'] + tempfiles
            else:
                cmd = ['gdal_merge.py', "-o '{}'".format(outfile), '-separate'] + tempfiles
            ret = os.popen(' '.join(cmd)).read()
            ##sp = subprocess.run(' '.join(cmd),shell=True, check=True, stdout=subprocess.PIPE)

        ## convert the RGB to COG
        if os.path.exists(outfile_temp) and oformat == 'COG':
            dt = gdal.Translate(outfile, outfile_temp, options = options_creation)
            if remove_temp_files: os.remove(outfile_temp.strip("'"))
        ## end RGB COG

        if os.path.exists(outfile): print('Wrote {}'.format(outfile))

        ## remove tempfiles
        if remove_temp_files:
            for file in tempfiles:
                os.remove(file.strip("'"))
    ## close file
    gem.close()
