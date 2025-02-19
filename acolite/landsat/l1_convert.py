## def l1_convert
## converts landsat bundle data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-07 (QV) added support for tile merging and extending the size of the output file to the requested limit
##                2021-02-09 (QV) added cross zone warping support (test)
##                2021-02-10 (QV) fixed cross zone warping and added support for full tile warping
##                2021-02-11 (QV) added checks for merging tiles of the same sensor and close in time
##                2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2023-02-04 (QV) added QA band output
##                2023-04-20 (QV) fix for changed extension case
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-03-28 (QV) added multi-tile merging
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit
##                2025-02-02 (QV) removed percentiles
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None,

                output_pan = True,
                output_pan_ms = True,
                output_thermal = True,

                usgs_reflectance = True,

                check_sensor = True,
                check_time = True,
                max_merge_time = 600, # seconds
                ):

    import os, glob, datetime, dateutil.parser, time, copy
    import acolite as ac
    import scipy.ndimage
    import numpy as np
    t0 = time.time()

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

    new = True
    new_pan = True
    warp_to = None
    warp_to_pan = None

    ofile = None
    ofiles = []

    for bundle in inputfile:
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        mtl = glob.glob('{}/{}'.format(bundle, '*MTL.txt'))
        ## add ALI MTL files
        mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1T.TXT'))
        mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1GST.TXT'))

        if len(mtl) == 0:
            if verbosity > 0: print('No metadata file found for {}'.format(bundle))
            continue
        else:
            mtl = mtl[0]

        ## read landsat metadata and check files
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        meta = ac.landsat.metadata_read(mtl)
        fmeta = ac.landsat.metadata_bands(bundle, meta)

        ## get relevant data from meta
        if 'PRODUCT_CONTENTS' in meta: ## COLL2
            pk = 'IMAGE_ATTRIBUTES'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PROJECTION_ATTRIBUTES'
            level = meta['PRODUCT_CONTENTS']['PROCESSING_LEVEL']
        elif 'PRODUCT_METADATA' in meta: ## COLL1
            pk = 'PRODUCT_METADATA'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PRODUCT_METADATA'
            level = meta['PRODUCT_METADATA']['DATA_TYPE']

        if level[0:2].upper() != 'L1':
            print('Only Level 1 Landsat data can be processed {} is {}'.format(bundle, level))
            continue

        spacecraft_id = meta[pk]['SPACECRAFT_ID']
        sensor_id = meta[pk]['SENSOR_ID']
        path = meta[pk]['WRS_PATH']
        row = meta[pk]['WRS_ROW']
        isodate = meta[pk]['DATE_ACQUIRED']+'T'+meta[pk]['SCENE_CENTER_TIME']
        global_dims = int(meta[rk]['REFLECTIVE_LINES']), int(meta[rk]['REFLECTIVE_SAMPLES'])

        ## some hard coded info
        azi_use_band = '5'
        sat = 'L{}'.format(spacecraft_id[-1])
        pan_scale = 2
        satellite_sensor = None
        if sensor_id in ['TM']:# Landsat 5
            sen = 'TM'
            pan_bands = []
            thermal_bands = ['6']
            if sat == 'L4':
                satellite_sensor = '{}_{}'.format(sat,sen)
                sat = 'L5' ## to use the same LUT as L5
        elif sensor_id in ['ETM']:# Landsat 7
            sen = 'ETM'
            pan_bands = ['8']
            thermal_bands = ['6_VCID_1', '6_VCID_2']
        elif sensor_id in ['OLI', 'OLI_TIRS']:# Landsat 8
            sen = 'OLI'
            pan_bands = ['8']
            thermal_bands = ['10', '11']
        elif sensor_id in ['ALI']:# EO1/ALI
            pan_scale = 3
            sat = 'EO1'
            sen = 'ALI'
            pan_bands = ['1']
            thermal_bands = []
        elif sensor_id in ['MSS']:# Landsat 1--5
            sen = 'MSS'
            pan_bands = []
            # 3 has one
            if sat == 'L3':
                thermal_bands = ['8']
            else:
                thermal_bands = []
            if sat in ['L4', 'L5']:
                azi_use_band = '4'
        else:
            print(spacecraft_id, sensor_id)
            print('Not configured')
            continue

        ## sensor name
        sensor = '{}_{}'.format(sat,sen)

        ## merge sensor specific settings
        if new:
            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults
            verbosity = setu['verbosity']
            if output is None: output = setu['output']

            ## check if merging settings make sense
            extend_region = setu['extend_region']
            if (setu['merge_tiles']):
                if (setu['limit'] is None):
                    if not setu['merge_full_tiles']:
                        if verbosity > 0: print("Merging tiles without ROI limit, merging to first tile extent")
                    else:
                        if verbosity > 0: print("Merging tiles without ROI limit, merging to all tiles extent")
                        dct_tiles = ac.landsat.multi_tile_extent(inputfile, dct = None)
                else:
                    extend_region = True

        ## sub is set to None
        sub = None

        ## scene average geometry
        vza = 0
        sza = 90-float(meta[ik]['SUN_ELEVATION'])
        saa = float(meta[ik]['SUN_AZIMUTH'])
        ## compute view zenith angle
        r,l,t,b, nadir_top, nadir_bottom,nadir_middle = ac.landsat.image_corners(bundle, meta, use_band = azi_use_band)
        vaa = ac.shared.azimuth_two_points(nadir_top[0],nadir_top[1],nadir_bottom[0],nadir_bottom[1])
        raa = np.abs(saa-vaa)
        while raa > 180: raa = np.abs(360 - raa)
        se_distance = float(meta[ik]['EARTH_SUN_DISTANCE'])

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## gains
        gains_dict = None
        if setu['gains'] & (setu['gains_toa'] is not None):
            if len(setu['gains_toa']) == len(rsr_bands):
                gains_dict = {b: float(setu['gains_toa'][ib]) for ib, b in enumerate(rsr_bands)}

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims':global_dims,
                 'sza':sza, 'vza':vza, 'vaa':vaa, 'saa':saa, 'raa':raa, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'wrs_path': path, 'wrs_row': row,
                 'acolite_file_type': 'L1R'}
        ## track L4 satellite sensor to use same LUT as L5, but not duplicate LUT space
        if satellite_sensor is not None: gatts['satellite_sensor'] = satellite_sensor
        if setu['merge_tiles']:
            gatts['tile_code'] = 'merged'
        else:
            gatts['tile_code'] = '{}{}'.format(gatts['wrs_path'].zfill(3),gatts['wrs_row'].zfill(3))

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), gatts['tile_code'])
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ## output file information
        if (setu['merge_tiles'] is False) | (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile
        elif (setu['merge_tiles']) & (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

        ## check if we should merge these tiles
        if (setu['merge_tiles']) & (not new) & (os.path.exists(ofile)):
                fgatts = ac.shared.nc_gatts(ofile)
                if (check_sensor) & (fgatts['sensor'] != gatts['sensor']):
                    print('Sensors do not match, skipping {}'.format(bundle))
                    continue
                if check_time:
                    tdiff = dateutil.parser.parse(fgatts['isodate'])-dateutil.parser.parse(gatts['isodate'])
                    tdiff = abs(tdiff.days*86400 + tdiff.seconds)
                    if (tdiff > max_merge_time):
                        print('Time difference too large, skipping {}'.format(bundle))
                        continue

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]
            if b in fmeta:
                fmeta[b]['f0'] = f0_b[b]/10
                fmeta[b]['se_distance'] = gatts['se_distance']

        ## get scene projection and extent
        dct = ac.landsat.projection(meta)

        ## full scene
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        if (sub is None) & (setu['limit'] is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
        else:
            if extend_region:
                print("Can't extend region if no ROI limits given")
                extend_region = False

        ## remove warp_to from previous run if merge_tiles is not set
        if (setu['merge_tiles'] is False): warp_to = None
        if sub is None:
            sub_pan = None
            ## determine warping target
            if (warp_to is None):
                if (setu['merge_tiles'] & setu['merge_full_tiles']): ## warp to all tile extent
                    dct_prj = {k:dct_tiles[k] for k in dct_tiles}
                else: ## warp to current/first tile
                    dct_prj = {k:dct[k] for k in dct}
        else:
            pan_dims = sub[3]*pan_scale, sub[2]*pan_scale
            sub_pan = [s*pan_scale for s in sub]
            gatts['sub'] = sub
            gatts['pan_sub'] = sub_pan
            gatts['limit'] = setu['limit']

            ## get the target NetCDF dimensions and dataset offset
            if (warp_to is None):
                if (extend_region): ## include part of the roi not covered by the scene
                    dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                else: ## just include roi that is covered by the scene
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
        ## end cropped

        ## get projection info for netcdf
        if setu['netcdf_projection']:
            nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=False)
            ## PAN band projection - not used but why not compute it
            dct_prj_pan = {k: dct_prj[k] for k in dct_prj}
            dct_prj_pan['pixel_size'] = dct_prj_pan['pixel_size'][0]/pan_scale, dct_prj_pan['pixel_size'][1]/pan_scale
            dct_prj_pan['xdim'] *= pan_scale
            dct_prj_pan['ydim'] *= pan_scale
            nc_projection_pan = ac.shared.projection_netcdf(dct_prj_pan, add_half_pixel=False)
        else:
            nc_projection = None
            nc_projection_pan = None

        ## save projection keys in gatts
        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = copy.copy(dct_prj[k])

        ## new after S2 changes 2021-10-28
        xyr = [min(dct_prj['xrange']),
               min(dct_prj['yrange']),
               max(dct_prj['xrange']),
               max(dct_prj['yrange']),
               dct_prj['proj4_string']]
        xyr_pan = [min(dct_prj['xrange']),
                   min(dct_prj['yrange']),
                   max(dct_prj['xrange']),
                   max(dct_prj['yrange']),
                   dct_prj['proj4_string']]

        ## warp settings for read_band
        res_method = 'near'
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)
        warp_to_pan = (dct_prj['proj4_string'], xyr_pan, dct_prj['pixel_size'][0]/pan_scale,dct_prj['pixel_size'][1]/pan_scale, res_method)

        ## store scene and output dimensions
        gatts['ydim'], gatts['xdim'] = dct['ydim'], dct['xdim']
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']
        gatts['pan_dims'] =  dct_prj['dimensions'][0]*pan_scale, dct_prj['dimensions'][1]*pan_scale

        ## new file for every bundle if not merging
        if (setu['merge_tiles'] is False):
            new = True
            new_pan = True

        if new: ## half pixel offset for writing geotiff
            gatts['xrange'][0]-=gatts['pixel_size'][0]/2
            gatts['yrange'][0]-=gatts['pixel_size'][1]/2
            gatts['xrange'][1]-=gatts['pixel_size'][0]/2
            gatts['yrange'][1]-=gatts['pixel_size'][1]/2

        ## copy thermal constants to metadata
        mts = ['LEVEL1_THERMAL_CONSTANTS', 'TIRS_THERMAL_CONSTANTS', 'THERMAL_CONSTANTS'] ## Coll2, Coll1 L8, Coll1 L5/7
        for mt in mts:
            if mt in meta:
                for k in meta[mt]:
                    if k not in gatts: gatts[k] = float(meta[mt][k])


        ## if we are clipping to a given polygon get the clip_mask here
        if setu['polygon_clip']:
            clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub=False)
            clip_mask = clip_mask.astype(bool) == False

        ## create new file
        if new:
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.nc_projection = nc_projection
            new = False

        ## start the conversion
        ## write geometry
        if ('VAA' in fmeta) & ('SAA' in fmeta) & ('VZA' in fmeta) & ('SZA' in fmeta):
            if verbosity > 1: print('Reading per pixel geometry')
            sza = ac.shared.read_band(fmeta['SZA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
            mus = np.cos(sza*(np.pi/180.)) ## per pixel cos sun zenith
            if (setu['output_geometry']):
                saa = ac.shared.read_band(fmeta['SAA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                vza = ac.shared.read_band(fmeta['VZA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                vaa = ac.shared.read_band(fmeta['VAA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                mask = (vaa == 0) * (vza == 0) * (saa == 0) * (sza == 0)
                vza[mask] = np.nan
                sza[mask] = np.nan
                ## clip geometry data
                if (setu['polygon_clip']):
                    sza[clip_mask] = np.nan
                    saa[clip_mask] = np.nan
                    vza[clip_mask] = np.nan
                    vaa[clip_mask] = np.nan

                raa = np.abs(saa-vaa)
                tmp = np.where(raa>180)
                raa[tmp]=np.abs(360 - raa[tmp])
                raa[mask] = np.nan
                vaa = None
                saa = None
                mask = None
                gemo.write('raa', raa, replace_nan=True,)
                if verbosity > 1: print('Wrote raa')
                gemo.write('vza', vza, replace_nan=True,)
                if verbosity > 1: print('Wrote vza')
                gemo.write('sza', sza, replace_nan=True,)
                if verbosity > 1: print('Wrote sza')
                sza = None
                vza = None
        else:
            mus = np.asarray(gatts['mus'])  ## average cos sun zenith
            #mus.shape+=(1,1)

        if os.path.exists(ofile) & (not new):
            gemo.datasets_read()
            datasets = gemo.datasets
            #datasets = ac.shared.nc_datasets(ofile)
        else:
            datasets = []

        ## write lat/lon
        if (setu['output_geolocation']):
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=False)
                gemo.write('lon', lon)
                if verbosity > 1: print('Wrote lon')
                gemo.write('lat', lat)
                if verbosity > 1: print('Wrote lat')

        ## write x/y
        if (setu['output_xy']):
            if ('xm' not in datasets) or ('ym' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=False)
                gemo.write('xm', x)
                x = None
                if verbosity > 1: print('Wrote xm')
                gemo.write('ym', y)
                y = None
                if verbosity > 1: print('Wrote ym')

        ## write TOA bands
        if verbosity > 1: print('Converting bands')
        for b in fmeta:
            if '.TIF' not in fmeta[b]['FILE'].upper(): continue
            if b in setu['landsat_qa_bands']: continue
            if os.path.exists(fmeta[b]['FILE']):
                if b in waves_names:
                    pan = False
                    if b in pan_bands: ## pan band
                        if (not output_pan) & (not output_pan_ms): continue
                        pan = True
                        mus_pan = scipy.ndimage.zoom(mus, zoom=pan_scale, order=1) if len(np.atleast_1d(mus))>1 else mus * 1
                        data = ac.landsat.read_toa(fmeta[b], sub=sub_pan, mus=mus_pan, warp_to=warp_to_pan)
                        mus_pan = None
                    else: ## not a pan band
                        data = ac.landsat.read_toa(fmeta[b], sub=sub, mus=mus, warp_to=warp_to)
                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}
                    for k in fmeta[b]: ds_att[k] = fmeta[b][k]

                    if setu['gains'] & (gains_dict is not None):
                        ds_att['toa_gain'] = gains_dict[b]
                        data *= ds_att['toa_gain']
                        if verbosity > 1: print('Converting bands: Applied TOA gain {} to {}'.format(ds_att['toa_gain'], ds))

                    if output_pan & pan:
                        ## write output
                        ofile_pan = ofile.replace('_L1R.nc', '_L1R_pan.nc')
                        if new_pan:
                            gemop = ac.gem.gem(ofile_pan, new = True)
                            gemop.gatts = {k: gatts[k] for k in gatts}
                            gemop.nc_projection = nc_projection_pan
                            new_pan = False
                        gemop.write(ds, data, ds_att = ds_att, replace_nan = True)
                        if verbosity > 1: print('Converting bands: Wrote {} to separate L1R_pan'.format(ds))

                    ## prepare for low res output
                    if output_pan_ms & pan: data = scipy.ndimage.zoom(data, zoom=1/pan_scale, order=1)

                    ## clip data
                    if (setu['polygon_clip']): data[clip_mask] = np.nan

                    ## write to ms file
                    gemo.write(ds, data, ds_att = ds_att, replace_nan = True)
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
                else:
                    if b in thermal_bands:
                        if output_thermal:
                            ds = 'bt{}'.format(b).lower()
                            ds_att = {'band':b}
                            for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                            data = ac.landsat.read_toa(fmeta[b], sub=sub, warp_to=warp_to)

                            ## clip data
                            if (setu['polygon_clip']): data[clip_mask] = np.nan
                            gemo.write(ds, data, ds_att = ds_att, replace_nan = True)
                            if verbosity > 1: print('Converting bands: Wrote {}'.format(ds))
                    else:
                        continue

        ## output quality assesment bands
        if setu['landsat_qa_output']:
            if verbosity > 1: print('Writing QA bands')
            for b in fmeta:
                if '.TIF' not in fmeta[b]['FILE']: continue
                if b not in setu['landsat_qa_bands']: continue
                ds = 'QA_{}'.format(b)
                ds_att = {'band':b}
                for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                data = ac.shared.read_band(fmeta[b]['FILE'], sub=sub, warp_to=warp_to)
                gemo.write(ds, data, ds_att = ds_att, replace_nan = True)
                if verbosity > 1: print('Writing QA bands: Wrote {}'.format(ds))

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        gemo.close()
        try:
            gemop.close()
        except:
            pass

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
