## def l1_convert
## convert PACE OCI NetCDF file to ACOLITE L1R NetCDF
##
## written by Quinten Vanhellemont, RBINS
## 2024-02-21
## modifications: 2024-04-15 (QV) updated for first light data, integrated as function
##                2024-04-16 (QV) use new gem NetCDF handling
##                2024-07-01 (QV) added L2 conversion
##                2024-07-03 (QV) store band irradiance
##                2024-07-11 (QV) changed attributes loading, added instrument_gain for SWIR
##                2024-07-22 (QV) include SWIR RSR
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming
##                2025-02-13 (QV) added tile merging, flip data
##                2025-03-03 (QV) added support for gains
##                2025-03-17 (QV) fixed application of gains
##                2025-03-18 (QV) fix for when no limit or sub is supplied

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import dateutil.parser, time
    import numpy as np
    import acolite as ac

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if setu['verbosity'] > 1: print('Starting conversion of {} scene{}'.format(nscenes, '' if nscenes==1 else 's'))

    ## list to store output files
    ofiles = []

    if setu['merge_tiles'] & (nscenes == 1):
        if setu['verbosity'] > 1: print('One scene provided, and merge_tiles=True. Setting merge_tiles=False.')
        setu['merge_tiles'] = False
        ac.settings['run']['merge_tiles'] = False

    ## test if we need to merge
    if setu['merge_tiles']:
        if setu['verbosity'] > 1: print('Testing whether {} scene{} can be merged'.format(nscenes, '' if nscenes==1 else 's'))
        ret = ac.pace.pace_merge_test(inputfile, limit = setu['limit'])
        if ret is None: ## return with no result
            return(ofiles, setu)
        else: ## unpack returns and sort inputfiles by time
            sub_merged, data_shape_merged, sort_bundles, crop_in, crop_out = ret
            inputfile = [inputfile[bi] for bi in sort_bundles]

    new = True
    ## run through inputfiles
    for bi, file in enumerate(inputfile):
        ## create new file if not merging
        if not setu['merge_tiles']: new = True

        ## read attributes
        igatts = ac.shared.nc_gatts(file)

        try:
            platform = igatts['platform'] ## simulated data only?
        except:
            if 'PACE' in igatts['title']: platform = 'PACE'
        instrument = igatts['instrument']
        sensor = '{}_{}'.format(platform, instrument)

        ## update settings
        setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(file)

        if setu['limit'] is None:
            print('Warning processing of PACE/OCI data recommended with small ROI')
            print('Supply limit or polygon for processing!')

        if igatts['processing_level'] == 'L1B':
            geo_group = 'geolocation_data'
            acolite_file_type = 'L1R'
            level1 = True
        elif igatts['processing_level'] == 'L2':
            geo_group = 'navigation_data'
            acolite_file_type = 'L2A'
            level1 = False
        else:
            print('Format of {} not supported'.format(file))
            continue

        if not setu['merge_tiles']:
            sub = None
            ## read lat and lon
            lon = ac.shared.nc_data(file, 'longitude', group=geo_group)
            lat = ac.shared.nc_data(file, 'latitude', group=geo_group)

            ## get subset
            if setu['sub'] is not None:
                sub = setu['sub']
            elif setu['limit'] is not None:
                sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if (setu['limit'] is not None) & (sub is None):
                print('Limit not in scene {}'.format(file))
                continue
        else:
            if setu['limit'] is None:
                sub = None
                data_shape = data_shape_merged[0],  data_shape_merged[1]
            else:
                sub = crop_in[bi][0], crop_in[bi][2], crop_in[bi][1]-crop_in[bi][0], crop_in[bi][3]-crop_in[bi][2]
                data_shape = sub_merged[3], sub_merged[2]

        ## make new output file
        if new:
            isodate = igatts['time_coverage_start']
            time = dateutil.parser.parse(isodate)

            ## output attributes
            gatts = {'sensor': sensor, 'isodate': time.isoformat()}
            gatts['acolite_file_type'] = acolite_file_type
            oname =  '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
            if setu['merge_tiles']: oname+='_merged'
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

            ## create output file
            gemo = ac.gem.gem(ofile, new = True)
            gemo.verbosity = setu['verbosity']
            new = False

        ## read write lat
        if setu['merge_tiles']: ## update merged version
            if 'lat' not in gemo.data_mem: gemo.data_mem['lat'] = np.zeros(data_shape) + np.nan
            gemo.data_mem['lat'][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = \
                 np.flipud(ac.shared.nc_data(file, 'latitude', group = geo_group, sub = sub))
        elif sub is not None: ## read cropped version
            gemo.data_mem['lat'] = np.flipud(ac.shared.nc_data(file, 'latitude', group = geo_group, sub = sub))
        else: ## store full version
            gemo.data_mem['lat'] = np.flipud(lat)
            lat = None

        ## read write lon
        if setu['merge_tiles']: ## update merged version
            if 'lon' not in gemo.data_mem: gemo.data_mem['lon'] = np.zeros(data_shape) + np.nan
            gemo.data_mem['lon'][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = \
                np.flipud(ac.shared.nc_data(file, 'longitude', group = geo_group, sub = sub))
        elif sub is not None: ## read cropped version
            gemo.data_mem['lon'] = np.flipud(ac.shared.nc_data(file, 'longitude', group = geo_group, sub = sub))
        else: ## store full version
            gemo.data_mem['lon'] = np.flipud(lon)
            lon = None
        ## end write lat/lon

        ## level1 data
        if level1:
            ## read write geometry
            geo_datasets = {'sza': 'solar_zenith','saa': 'solar_azimuth',
                            'vza': 'sensor_zenith', 'vaa': 'sensor_azimuth'}
            for ko in geo_datasets:
                if setu['merge_tiles']:
                    if ko not in gemo.data_mem: gemo.data_mem[ko] = np.zeros(data_shape) + np.nan
                    gemo.data_mem[ko][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = \
                        np.flipud(ac.shared.nc_data(file, geo_datasets[ko], group = geo_group, sub = sub))
                else:
                    gemo.data_mem[ko] = np.flipud(ac.shared.nc_data(file, geo_datasets[ko], group = geo_group, sub = sub))

            ## read SWIR RSR
            rsrd_swir = ac.shared.rsr_dict('PACE_OCI_SWIR')

            ## track gain
            band_gains = []
            gain_index = 0

            ## read band data
            band_waves = []
            band_widths = []
            band_irradiance = []
            for det in ['blue', 'red', 'SWIR']:
                print('Reading data from detector {}'.format(det))
                f0_det, f0_att = ac.shared.nc_data(file, '{}_solar_irradiance'.format(det), \
                                                   group = 'sensor_band_parameters', attributes = True)
                wv_det, wv_att = ac.shared.nc_data(file, '{}_wavelength'.format(det), \
                                                   group = 'sensor_band_parameters', attributes = True)

                if det == 'SWIR':
                    bp_det, bp_att = ac.shared.nc_data(file, '{}_bandpass'.format(det), \
                                                   group = 'sensor_band_parameters', attributes = True)
                else:
                    bp_det = np.zeros(len(wv_det))+5.0

                ## read detector data
                rhot, rhot_att = ac.shared.nc_data(file, 'rhot_{}'.format(det), \
                                                   group = 'observation_data', attributes = True, sub = sub)
                print(det, rhot.shape, len(wv_det))

                for wi, wave in enumerate(wv_det):
                    if not np.isfinite(wave): continue

                    att = {'f0': f0_det[wi], 'wave': wave, 'wave_name': '{:.0f}'.format(wave), 'width': bp_det[wi]}

                    ## track gains
                    if setu['gains']:
                        att['gain'] = setu['gains_toa'][gain_index]
                        band_gains.append(att['gain'])
                        gain_index +=1
                        #print(wave, gain_index, att['gain'])

                    ## SWIR instrument gain and update band name
                    ## presumed same order as PACE_OCI_L1B_LUT_RSR_baseline_1.1.1.nc
                    ## 0 - 939.71497
                    ## 1 - 1038.315
                    ## 2 - 1250.375  standard gain
                    ## 3 - 1248.5525 high gain
                    ## 4 - 1378.165
                    ## 5 - 1619.625  standard gain
                    ## 6 - 1618.0349 high gain
                    ## 7 - 2130.5923
                    ## 8 - 2258.43
                    if (det == 'SWIR'):
                        if (wi in [3, 6]):
                            att['instrument_gain'] = 'high'
                        else:
                            att['instrument_gain'] = 'standard'
                        swir_b = rsrd_swir['PACE_OCI_SWIR']['rsr_bands'][wi]
                        att['wave_name'] = rsrd_swir['PACE_OCI_SWIR']['wave_name'][swir_b]
                    ## end SWIR gain and band name

                    band_waves.append(att['wave'])
                    band_widths.append(att['width'])
                    band_irradiance.append(att['f0'])

                    ds_name = 'rhot_{}_{}'.format(det, att['wave_name'])
                    if setu['merge_tiles']:
                        if ds_name not in gemo.data_mem:
                            gemo.data_mem[ds_name] = np.zeros(data_shape) + np.nan
                            gemo.data_att[ds_name] = {k: att[k] for k in att}
                        gemo.data_mem[ds_name][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = np.flipud(rhot[wi, :,:])
                    else:
                        gemo.data_mem[ds_name] = np.flipud(rhot[wi, :,:])
                        gemo.data_att[ds_name] = att
                rhot = None

            ## update attributes
            gatts['band_waves'] = band_waves
            gatts['band_widths'] = band_widths
            gatts['band_irradiance'] = band_irradiance
            if setu['gains']: gatts['band_gains'] = band_gains

            ## compute relative azimuth
            if (not setu['merge_tiles']) | (bi == len(inputfile)-1):
                raa = np.abs(gemo.data('saa') - gemo.data('vaa'))
                tmp = np.where(raa>180)
                raa[tmp]=np.abs(360 - raa[tmp])
                gemo.data_mem['raa'] = raa
                raa = None
        ## end level1 data

        ## level2 data
        if not level1:
            data_group = 'geophysical_data'

            ## get sensor band parameters
            from netCDF4 import Dataset
            with Dataset(file) as nc:
                band_pars = list(nc.groups['sensor_band_parameters'].variables.keys())
                band_atts = {}
                for p in band_pars:
                    band_atts[p] = nc.groups['sensor_band_parameters'][p][:]

            ## get datasets
            for ds_name in ['l2_flags', 'avw', 'aot_865', 'angstrom', 'Rrs']:
                d, att = ac.shared.nc_data(file, ds_name, group=data_group, sub=sub, attributes = True, axis_3d = 2)
                if d.dtype in [np.float32]:
                    d[d.mask] = np.nan
                if ds_name in ['Rrs']:
                    print(d.shape)
                    for wi in range(d.shape[2]):
                        ds_att = {k: att[k] for k in att}
                        for k in band_atts: ds_att[k] = band_atts[k][wi]
                        ds_out = '{}_{}'.format(ds_name, ds_att['wavelength_3d'])
                        if setu['merge_tiles']:
                            if ds_out not in gemo.data_mem:
                                gemo.data_mem[ds_out] = np.zeros(data_shape) + np.nan
                                gemo.data_att[ds_out] = ds_att
                            gemo.data_mem[ds_out][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = np.flipud(d[:,:, wi])
                        else:
                            gemo.data_mem[ds_out] = np.flipud(d[:,:, wi])
                            gemo.data_att[ds_out] = ds_att
                else:
                    if setu['merge_tiles']:
                        if ds_name not in gemo.data_mem:
                            gemo.data_mem[ds_name] = np.zeros(data_shape) + np.nan
                            gemo.data_att[ds_name] = att
                        gemo.data_mem[ds_name][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = np.flipud(d)
                    else:
                        gemo.data_mem[ds_name] = np.flipud(d)
                        gemo.data_att[ds_name] = att
        ## end level2 data

        ## write data if not merging, or if last image in merging
        if (not setu['merge_tiles']) | (bi == len(inputfile)-1):
            write_ds = list(gemo.data_mem.keys())
            for ds in write_ds:
                if 'rhot_' in ds:
                    if 'gain' in gemo.data_att[ds]:
                        gemo.data_mem[ds] *= gemo.data_att[ds]['gain']
                        print('Applied gain {} to {}'.format(gemo.data_att[ds]['gain'], ds))
                gemo.write_ds(ds, clear = True)
            ## update attributes
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.gatts_update()
            ## close file
            gemo.close()
            ofiles.append(ofile)

    return(ofiles, setu)
