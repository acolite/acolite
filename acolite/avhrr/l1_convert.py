## def l1_convert
## converts AVHRR L1B/1C European Data Set bundle to ACOLITE L1R
## written by Quinten Vanhellemont, RBINS
## 2024-11-14
## modifications: 2024-11-16 (QV) added data flip for ascending orbit only
##                2024-11-17 (QV) added 6 channel support
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output=None, settings = None):
    import numpy as np
    import dateutil.parser, os
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

    ## get F0 for radiance -> reflectance computation
    #f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

    ofiles = []
    for file in inputfile:
        image_file, meta_file = ac.avhrr.bundle_test(file)
        datasets = ac.shared.nc_datasets(image_file)
        print(datasets)

        metadata = ac.avhrr.metadata(meta_file)
        sdate = dateutil.parser.parse(metadata['start_date'])
        edate = dateutil.parser.parse(metadata['end_date'])
        sensor = metadata['sensor']
        sensor_tir = '{}_TIR'.format(sensor)

        print(metadata['eop:orbitDirection'])
        if metadata['eop:orbitDirection'] == 'DESCENDING':
            flip = False
        elif metadata['eop:orbitDirection'] == 'ASCENDING':
            flip = True
        else:
            print('eop:orbitDirection {} not configured'.format(metadata['eop:orbitDirection']))
            flip = True

        ## read RSR
        rsrd = ac.shared.rsr_dict(sensor)[sensor]
        #rsrd_tir = ac.shared.rsr_dict(sensor_tir)[sensor_tir]

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.basename(file)

        gatts = {}
        time = dateutil.parser.parse(metadata['start_date'])
        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)

        ## lon and lat
        lat = ac.shared.nc_data(image_file, 'latitude')
        lat[lat.mask] = np.nan
        lat = ac.shared.fillnan(lat)
        lon = ac.shared.nc_data(image_file, 'longitude')
        lon[lon.mask] = np.nan
        lon = ac.shared.fillnan(lon)

        sub = None
        if setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit outside of scene {}'.format(file))
                continue
            ## crop to sub
            lat = lat[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
            lon = lon[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
        ## end read geolocation

        ## read geometry
        vaa, saa = None, None
        #vaa = ac.avhrr.read_fill(image_file,'sensor_azimuth_angle', sub = sub)
        vza = ac.avhrr.read_fill(image_file,'sensor_zenith_angle', sub = sub)
        #saa = ac.avhrr.read_fill(image_file,'solar_azimuth_angle', sub = sub)
        sza = ac.avhrr.read_fill(image_file,'solar_zenith_angle', sub = sub)
        raa = ac.avhrr.read_fill(image_file,'sun_sensor_azimuth_difference_angle', sub = sub)

        gatts['vza'] = np.nanmean(np.abs(vza))
        gatts['raa'] = np.nanmean(np.abs(raa))
        gatts['sza'] = np.nanmean(np.abs(sza))

        #mus = np.cos(np.radians(sza))
        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        gatts['sensor'] = sensor
        gatts['isodate'] = time.isoformat()
        gatts['acolite_file_type'] = 'L1R'

        ## output file name
        oname  = '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## set up output file
        gemo = ac.gem.gem(ofile, new=True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        if flip:
            lon = np.flip(lon)
            lat = np.flip(lat)
            sza = np.flip(sza)
            raa = np.flip(raa)

        if (setu['output_geolocation']):
            if verbosity > 1: print('Writing geolocation lon/lat')
            gemo.write('lon', lon)
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
            gemo.write('lat', lat)
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))

        ## write geometry
        if (setu['output_geometry']):
            if verbosity > 1: print('Writing geometry')
            gemo.write('vza', vza)
            if verbosity > 1: print('Wrote vza ({})'.format(vza.shape))
            vza = None
            if vaa is not None:
                if flip: vaa = np.flip(vaa)
                gemo.write('vaa', vaa)
                if verbosity > 1: print('Wrote vaa ({})'.format(vaa.shape))
                vaa = None

            gemo.write('sza', sza)
            if verbosity > 1: print('Wrote sza ({})'.format(sza.shape))
            if saa is not None:
                if flip: saa = np.flip(saa)
                gemo.write('saa', saa)
                if verbosity > 1: print('Wrote saa ({})'.format(saa.shape))
                saa = None
            gemo.write('raa', raa)
            if verbosity > 1: print('Wrote raa ({})'.format(raa.shape))
            raa = None

        ## write TOA data
        bands_datasets = ['reflectance_channel_1',
                          'reflectance_channel_2',
                          'reflectance_channel_3a',
                         'brightness_temperature_channel_3',
                         'brightness_temperature_channel_3b',
                         'brightness_temperature_channel_4',
                         'brightness_temperature_channel_5']
        for bi, b in enumerate(bands_datasets):
            if b not in datasets: continue
            bid = b.split('_')[-1]

            ds_att = {}
            if 'reflectance' in b:
                ds = 'rhot_{}'.format(rsrd['wave_name'][bid])
                ds_att['wavelength'] = rsrd['wave_nm'][bid]
                factor = (1./100)
            elif 'brightness_temperature' in b:
                if not setu['output_bt']: continue
                ds = 'bt_{}'.format(bid)
                factor = 1.

            print('Reading {}'.format(b))

            data = ac.avhrr.read_fill(image_file, b, sub = sub)
            data *= factor

            ## write toa data
            if flip: data = np.flip(data)
            gemo.write(ds, data, ds_att = ds_att)
            print('Wrote {}'.format(ds))

        ## close files
        gemo.close()
        gemo = None

        ofiles.append(ofile)
    return(ofiles, setu)
