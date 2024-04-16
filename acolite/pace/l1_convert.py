## def l1_convert
## convert PACE OCI NetCDF file to ACOLITE L1R NetCDF
##
## written by Quinten Vanhellemont, RBINS
## 2024-02-21
## modifications: 2024-04-15 (QV) updated for first light data, integrated as function
##                2024-04-16 (QV) use new gem NetCDF handling

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import dateutil.parser, time
    import numpy as np
    import acolite as ac

    ## get run verbosity
    verbosity = ac.settings['run']['verbosity']

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

        ## read lat and lon
        lon = ac.shared.nc_data(file, 'longitude', group='geolocation_data')
        lat = ac.shared.nc_data(file, 'latitude', group='geolocation_data')

        ## get subset
        sub = ac.shared.geolocation_sub(lat, lon, limit)
        if (limit is not None) & (sub is None):
            print('Limit not in scene {}'.format(file))
            continue

        isodate = igatts['time_coverage_start']
        time = dateutil.parser.parse(isodate)

        ## output attributes
        gatts = {'sensor': sensor, 'isodate': time.isoformat()}
        gatts['acolite_file_type'] = 'L1R'
        gatts['obase'] = '{}_{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'), gatts['acolite_file_type'])
        ofile = '{}/{}.nc'.format(output, gatts['obase'])

        ## read write lat/lon
        if sub is not None: ## read cropped version
            lat = ac.shared.nc_data(file, 'latitude', group='geolocation_data', sub=sub)

        ## output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.write('lat', lat)
        lat = None

        if sub is not None: ## read cropped version
            lon = ac.shared.nc_data(file, 'longitude', group='geolocation_data', sub=sub)
        gemo.write('lon', lon)
        lon = None

        ## read write geometry
        sza = ac.shared.nc_data(file, 'solar_zenith', group='geolocation_data', sub=sub)
        gemo.write('sza', sza)
        sza = None
        vza = ac.shared.nc_data(file, 'sensor_zenith', group='geolocation_data', sub=sub)
        gemo.write('vza', vza)
        vza = None

        saa = ac.shared.nc_data(file, 'solar_azimuth', group='geolocation_data', sub=sub)
        vaa = ac.shared.nc_data(file, 'sensor_azimuth', group='geolocation_data', sub=sub)

        raa = np.abs(saa - vaa)
        tmp = np.where(raa>180)
        raa[tmp]=np.abs(360 - raa[tmp])
        gemo.write('saa', saa)
        gemo.write('vaa', vaa)
        saa, vaa = None, None

        gemo.write('raa', raa)
        raa = None

        ## read band data
        band_waves = []
        band_widths = []
        for det in ['blue', 'red', 'SWIR']:
            print('Reading data from detector {}'.format(det))
            f0_det = ac.shared.nc_data(file, '{}_solar_irradiance'.format(det), group='sensor_band_parameters')
            wv_det = ac.shared.nc_data(file, '{}_wavelength'.format(det), group='sensor_band_parameters')

            if det == 'SWIR':
                bp_det = ac.shared.nc_data(file, '{}_bandpass'.format(det), group='sensor_band_parameters')
            else:
                bp_det = np.zeros(len(wv_det))+5.0

            ## read detector data
            tmp = ac.shared.nc_data(file, 'rhot_{}'.format(det), group='observation_data', sub=sub)
            print(det, tmp.shape, len(wv_det))

            for bi, wave in enumerate(wv_det):
                if not np.isfinite(wave): continue

                att = {'f0': f0_det[bi], 'wave': wave, 'wave_name': '{:.0f}'.format(wave), 'width': bp_det[bi]}

                band_waves.append(att['wave'])
                band_widths.append(att['width'])

                ds_name = 'rhot_{}_{}'.format(det, att['wave_name'])
                gemo.write(ds_name, tmp[bi, :,:], ds_att = att)
                print('Wrote {}'.format(ds_name))
            tmp = None

        ## update attributes
        gatts['band_waves'] = band_waves
        gatts['band_widths'] = band_widths

        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()
        ofiles.append(ofile)

    return(ofiles, setu)
