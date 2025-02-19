## def l2_convert
## converts Sentinel-2 L2A Sen2Cor .SAFE bundle data to s2c_l2a NetCDF
## written by Quinten Vanhellemont, RBINS
## 2024-06-04
## modifications: 2024-07-30 (QV) added skip_bands to not include band 9 (or others if so desired)
##                2024-10-16 (QV) added RSR versioning support
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l2_convert(inputfile, output = None, settings = None, skip_bands = ['9']):
    import numpy as np
    import acolite as ac
    import os, dateutil.parser, time

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ofile = None
    ofiles = []
    sub = None

    for bundle in inputfile:
        t0 = time.time()
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        try:
            safe_files = ac.sentinel2.safe_test(bundle)
        except:
            print('File not recognised: {}'.format(bundle))
            continue


        if 'granules' not in safe_files:
            print('File not recognised: {}'.format(bundle))
            continue

        if len(safe_files['granules']) > 1:
            print('Multi granule files are no longer supported.')
            continue

        granule = safe_files['granules'][0]

        if verbosity > 1: print('Importing metadata from {}'.format(granule))
        grmeta = ac.sentinel2.metadata_granule(safe_files[granule]['metadata']['path'])
        meta, band_data= ac.sentinel2.metadata_scene(safe_files['metadata']['path'])

        ## get relevant data from meta
        if meta['SPACECRAFT_NAME'] == 'Sentinel-2A':
            sensor = 'S2A_MSI'
        elif meta['SPACECRAFT_NAME'] == 'Sentinel-2B':
            sensor = 'S2B_MSI'
        else:
            print('{} not supported'.format(meta['SPACECRAFT_NAME']))
            continue
        if meta['PROCESSING_LEVEL'] != 'Level-2A':
            print('Processing of {} Sentinel-2 {} data not supported by convert_l2'.format(bundle, meta['PROCESSING_LEVEL']))
            continue

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults
        verbosity = setu['verbosity']
        if output is None: output = setu['output']

        ## output metadata
        dtime = dateutil.parser.parse(grmeta['SENSING_TIME'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()
        dct = ac.sentinel2.projection(grmeta, s2_target_res=int(setu['s2_target_res']))
        global_dims = dct['dimensions']
        mgrs_tile = grmeta['TILE_ID'].split('_')[-2]

        ## read rsr for band names
        if setu['rsr_version'] is None:
            rsrf = ac.path+'/data/RSR/{}.txt'.format(sensor)
        else:
            rsrf = ac.path+'/data/RSR/{}_{}.txt'.format(sensor, setu['rsr_version'])
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims':global_dims,
                 'granule': granule, 'mgrs_tile': mgrs_tile, 'acolite_file_type': 'L2A'}
        gatts['tile_code'] = '{}'.format(gatts['mgrs_tile'])
        stime = dateutil.parser.parse(gatts['isodate'])

        ## output filename
        oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), gatts['tile_code'])
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]

        ## get scene projection and extent
        dct = ac.sentinel2.projection(grmeta, s2_target_res=int(setu['s2_target_res']))

        ## full scene
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        ## check crop
        if (sub is None) & (setu['limit'] is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
        if sub is None: ## full tile processing
             dct_prj = {k:dct[k] for k in dct}
        else:
            gatts['sub'] = sub
            gatts['limit'] = setu['limit']
            dct_prj = {k:dct_sub[k] for k in dct_sub}

        ## get projection info for netcdf
        if setu['netcdf_projection']:
            nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=True)
        else:
            nc_projection = None

        ## save projection keys in gatts
        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]

        ## with subsetting fix the offsets should not be required 2021-10-28
        xyr = [min(dct_prj['xrange']),
               min(dct_prj['yrange']),
               max(dct_prj['xrange']),
               max(dct_prj['yrange']),
               dct_prj['proj4_string']]

        ## warp settings for read_band
        res_method = 'average'
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']

        ## if we are clipping to a given polygon get the clip_mask here
        if setu['polygon_clip']:
            clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub=False)
            clip_mask = clip_mask.astype(bool) == False
            print('clip mask', clip_mask.shape)

        new = True
        if new:
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.nc_projection = nc_projection
            new = False

        ## write lat/lon
        if (setu['output_geolocation']):
            if verbosity > 1: print('Writing geolocation lon/lat')
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
            gemo.write('lon', lon)
            if verbosity > 1: print('Wrote lon {}'.format(lon.shape))
            lon = None
            gemo.write('lat', lat)
            if verbosity > 1: print('Wrote lat {}'.format(lat.shape))
            lat = None

        ##
        nodata = int(meta['NODATA'])
        quant = float(meta['BOA_QUANTIFICATION_VALUE'])
        dilate = setu['s2_dilate_blackfill']
        dilate_iterations = setu['s2_dilate_blackfill_iterations']

        if verbosity > 1: print('Converting bands')
        for bi, b in enumerate(rsr_bands):
            if skip_bands is not None:
                if b in skip_bands: continue
            Bn = 'B{}'.format(b)
            if Bn not in safe_files[granule]: continue
            if os.path.exists(safe_files[granule][Bn]['path']):
                if b in waves_names:
                    data = ac.shared.read_band(safe_files[granule][Bn]['path'], sub=sub, warp_to=warp_to)
                    data_mask = data == nodata
                    if dilate: data_mask = scipy.ndimage.binary_dilation(data_mask, iterations=dilate_iterations)
                    data = data.astype(np.float32)
                    if 'BOA_ADD_OFFSET' in band_data: data += band_data['BOA_ADD_OFFSET'][Bn]

                    data /= quant
                    data[data_mask] = np.nan
                    if (setu['polygon_clip']): data[clip_mask] = np.nan
                    ds = 'rhos_l2a_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}

                    ## write to ms file
                    gemo.write(ds, data, replace_nan = True, ds_att = ds_att)
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
            else:
                continue

        ## update attributes
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
