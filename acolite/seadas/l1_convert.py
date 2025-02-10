# def l1_convert
# converts SeaDAS L1B file to l1r NetCDF for acolite
# written by Quinten Vanhellemont, RBINS
# 2024-10-16
# initially for SeaHawk1_HawkEye, but can presumably be adapted to other L1B data
#
# modifications:  2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac

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

    ## tested platforms
    seadas_platforms = ['Seahawk1']

    ## run through inputfiles
    ofiles = []
    for fi, file in enumerate(inputfile):

        ## read attributes
        igatts = ac.shared.nc_gatts(file)

        if igatts['platform'] not in seadas_platforms: continue
        if igatts['processing_level'] not in ['L1B', 'L2']: continue

        ## set up sensor name based on platform and instrument
        platform = igatts['platform']
        instrument = igatts['instrument']
        if platform.lower() in ['seahawk1', 'seahawk2']:
            platform = 'SeaHawk{}'.format(platform[-1])
        if instrument.lower() in ['hawkeye']:
            instrument = 'HawkEye'
        sensor = '{}_{}'.format(platform, instrument)

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        ## get F0 for radiance -> reflectance computation
        ## F0 is present in 'sensor_band_parameters' nc group but empty
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

        ## read sensor rsr and convolve F0
        rsrd = ac.shared.rsr_dict(sensor)
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd[sensor]['rsr'])

        verbosity = setu['verbosity']
        if output is None: output = setu['output']

        sub = None
        geo_group = 'navigation_data'
        sensor_group = 'sensor_band_parameters'
        data_group = 'geophysical_data'
        if igatts['processing_level'] == 'L1B':
            acolite_file_type = 'L1R'
            level1 = True
        elif igatts['processing_level'] == 'L2':
            acolite_file_type = 'L2A'
            level1 = False
        else:
            print('Format of {} not supported'.format(file))
            continue

        ## get band waves
        waves = ac.shared.nc_data(file, 'wavelength', group = sensor_group).data

        ## read lat and lon
        lon = ac.shared.nc_data(file, 'longitude', group=geo_group)
        lat = ac.shared.nc_data(file, 'latitude', group=geo_group)

        ## get subset
        if setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])

        if (setu['limit'] is not None) & (sub is None):
            print('Limit not in scene {}'.format(file))
            continue

        isodate = igatts['time_coverage_start']
        time = dateutil.parser.parse(isodate)

        ## output attributes
        gatts = {'sensor': sensor, 'isodate': time.isoformat()}
        gatts['acolite_file_type'] = acolite_file_type
        oname  = '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## read cropped version
        if sub is not None:
            lat = ac.shared.nc_data(file, 'latitude', group=geo_group, sub=sub)
            lon = ac.shared.nc_data(file, 'longitude', group=geo_group, sub=sub)

        t0 = igatts['time_coverage_start']
        t1 = igatts['time_coverage_end']
        start_time = dateutil.parser.parse(t0)
        stop_time = dateutil.parser.parse(t1)
        half_time_diff = (stop_time-start_time)/2
        dt = start_time + half_time_diff
        isodate = dt.isoformat()

        ## compute scene center sun position
        clon = np.nanmedian(lon)
        clat = np.nanmedian(lat)
        cspos = ac.shared.sun_position(dt, clon, clat)
        spos = ac.shared.sun_position(dt, lon, lat)
        d = spos['distance']

        ## output attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['acolite_file_type'] = acolite_file_type
        gatts['isodate'] = dt.isoformat()

        ## centre sun position
        gatts['doy'] = dt.strftime('%j')
        gatts['se_distance'] = spos['distance']

        ## output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.write('lat', lat)
        lat = None

        gemo.write('lon', lon)
        lon = None

        ## level1 data
        if level1:
            ## read write geometry
            sza = ac.shared.nc_data(file, 'solz', group=geo_group, sub=sub)
            gemo.write('sza', sza)
            gatts['sza'] = np.nanmean(sza)

            ## for radiance conversion
            cossza = np.cos(np.radians(sza))

            sza = None
            vza = ac.shared.nc_data(file, 'senz', group=geo_group, sub=sub)
            gemo.write('vza', vza)
            gatts['vza'] = np.nanmean(vza)
            vza = None

            saa = ac.shared.nc_data(file, 'sola', group=geo_group, sub=sub)
            vaa = ac.shared.nc_data(file, 'sena', group=geo_group, sub=sub)

            raa = np.abs(saa - vaa)
            tmp = np.where(raa>180)
            raa[tmp]=np.abs(360 - raa[tmp])
            gemo.write('saa', saa)
            gemo.write('vaa', vaa)
            gatts['saa'] = np.nanmean(saa)
            gatts['vaa'] = np.nanmean(vaa)
            saa, vaa = None, None

            gemo.write('raa', raa)
            gatts['raa'] = np.nanmean(raa)
            raa = None

            ## read TOA data
            for wi, wave in enumerate(waves):
                band = rsrd[sensor]['rsr_bands'][wi]
                ds_att = {'wave_nm': rsrd[sensor]['wave_nm'][band],
                          'wave_name':rsrd[sensor]['wave_name'][band],
                          'f0': f0d[band]}

                data = ac.shared.nc_data(file, 'Lt_{}'.format(wave), group=data_group, sub=sub)

                ## write toa radiance
                if setu['output_lt']:
                    gemo.write('Lt_{}'.format(ds_att['wave_name']), data, ds_att = ds_att)
                    print('Wrote Lt_{}'.format(ds_att['wave_name']))

                ## compute reflectance
                data = data * (np.pi * d * d) / (ds_att['f0'] * cossza)

                ## write toa reflectance
                gemo.write('rhot_{}'.format(ds_att['wave_name']), data, ds_att = ds_att)
                data = None
                print('Wrote rhot_{}'.format(ds_att['wave_name']))

        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()
        ofiles.append(ofile)
    return(ofiles, setu)
