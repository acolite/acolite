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

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import dateutil.parser, time
    import numpy as np
    import acolite as ac

    ## get run verbosity
    verbosity = ac.settings['run']['verbosity']
    if (settings is None) & ('user' in ac.settings):
        settings = {k: ac.settings['user'][k] for k in ac.settings['user']}

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## run through inputfiles
    ofiles = []
    for fi, file in enumerate(inputfile):

        ## read attributes
        igatts = ac.shared.nc_gatts(file)

        try:
            platform = igatts['platform'] ## simulated data only?
        except:
            if 'PACE' in igatts['title']: platform = 'PACE'
        instrument = igatts['instrument']
        sensor = '{}_{}'.format(platform, instrument)

        ## parse settings
        setu = ac.acolite.settings.parse(sensor, settings=settings, merge=True)
        if output is None: output = setu['output']

        ## get ROI from user settings
        limit = setu['limit']

        ## check if ROI polygon is given
        poly = setu['polygon']
        clip, clip_mask = False, None
        if poly is not None:
            if os.path.exists(poly):
                try:
                    limit = ac.shared.polygon_limit(poly)
                    if verbosity > 1: print('Using limit from polygon envelope: {}'.format(limit))
                    clip = True
                except:
                    if verbosity > 1: print('Failed to import polygon {}'.format(poly))

        ## add limit buffer
        if (limit is not None) & (setu['limit_buffer'] is not None):
            print('Applying limit buffer {}'.format(setu['limit_buffer']))
            print('Old limit: {}'.format(limit))
            setu['limit_old'] = limit
            limit = limit[0] - setu['limit_buffer'], limit[1] - setu['limit_buffer'], \
                    limit[2] + setu['limit_buffer'], limit[3] + setu['limit_buffer']
            print('New limit: {}'.format(limit))

        sub = None
        if limit is None:
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

        ## read lat and lon
        lon = ac.shared.nc_data(file, 'longitude', group=geo_group)
        lat = ac.shared.nc_data(file, 'latitude', group=geo_group)

        ## get subset
        sub = ac.shared.geolocation_sub(lat, lon, limit)
        if (limit is not None) & (sub is None):
            print('Limit not in scene {}'.format(file))
            continue

        isodate = igatts['time_coverage_start']
        time = dateutil.parser.parse(isodate)

        ## output attributes
        gatts = {'sensor': sensor, 'isodate': time.isoformat()}
        gatts['acolite_file_type'] = acolite_file_type
        gatts['obase'] = '{}_{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'), gatts['acolite_file_type'])
        ofile = '{}/{}.nc'.format(output, gatts['obase'])

        ## read write lat/lon
        if sub is not None: ## read cropped version
            lat = ac.shared.nc_data(file, 'latitude', group=geo_group, sub=sub)

        ## output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.write('lat', lat)
        lat = None

        if sub is not None: ## read cropped version
            lon = ac.shared.nc_data(file, 'longitude', group=geo_group, sub=sub)
        gemo.write('lon', lon)
        lon = None

        ## level1 data
        if level1:
            ## read write geometry
            sza = ac.shared.nc_data(file, 'solar_zenith', group=geo_group, sub=sub)
            gemo.write('sza', sza)
            sza = None
            vza = ac.shared.nc_data(file, 'sensor_zenith', group=geo_group, sub=sub)
            gemo.write('vza', vza)
            vza = None

            saa = ac.shared.nc_data(file, 'solar_azimuth', group=geo_group, sub=sub)
            vaa = ac.shared.nc_data(file, 'sensor_azimuth', group=geo_group, sub=sub)

            raa = np.abs(saa - vaa)
            tmp = np.where(raa>180)
            raa[tmp]=np.abs(360 - raa[tmp])
            gemo.write('saa', saa)
            gemo.write('vaa', vaa)
            saa, vaa = None, None

            gemo.write('raa', raa)
            raa = None

            ## read SWIR RSR
            rsrd_swir = ac.shared.rsr_dict('PACE_OCI_SWIR')

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

                for bi, wave in enumerate(wv_det):
                    if not np.isfinite(wave): continue

                    att = {'f0': f0_det[bi], 'wave': wave, 'wave_name': '{:.0f}'.format(wave), 'width': bp_det[bi]}

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
                        if (bi in [3, 6]):
                            att['instrument_gain'] = 'high'
                        else:
                            att['instrument_gain'] = 'standard'
                        swir_b = rsrd_swir['PACE_OCI_SWIR']['rsr_bands'][bi]
                        att['wave_name'] = rsrd_swir['PACE_OCI_SWIR']['wave_name'][swir_b]
                    ## end SWIR gain and band name

                    band_waves.append(att['wave'])
                    band_widths.append(att['width'])
                    band_irradiance.append(att['f0'])

                    ds_name = 'rhot_{}_{}'.format(det, att['wave_name'])
                    gemo.write(ds_name, rhot[bi, :,:], ds_att = att)
                    print('Wrote {}'.format(ds_name))
                rhot = None

            ## update attributes
            gatts['band_waves'] = band_waves
            gatts['band_widths'] = band_widths
            gatts['band_irradiance'] = band_irradiance
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
                    for bi in range(d.shape[2]):

                        ds_att = {k: att[k] for k in att}
                        for k in band_atts: ds_att[k] = band_atts[k][bi]
                        ds_out = '{}_{}'.format(ds_name, ds_att['wavelength_3d'])
                        gemo.write(ds_out, d[:,:, bi], ds_att = ds_att)
                        print('Wrote {}'.format(ds_out))
                else:
                    gemo.write(ds_name, d, ds_att = att)
                    print('Wrote {}'.format(ds_name))
        ## end level2 data

        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()
        ofiles.append(ofile)

    return(ofiles, setu)
