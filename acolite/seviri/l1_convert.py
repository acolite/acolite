## def l1_convert
## convert SEVIRI nat file to ACOLITE L1R
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-12
## modifications: 2024-04-17 (QV) use new gem NetCDF handling
##                2024-04-18 (QV) added geolocation/geometry functions , added masks for negative Lt and sza > 90
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import dateutil.parser, time
    import numpy as np
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

    ## get F0 from EUM documentation
    with open(ac.config['data_dir'] + '/GEO/SEVIRI/f0.json', 'r', encoding='utf8') as f:
        f0s = json.load(f)
    band_names = {1: 'VIS06', 2: 'VIS08', 3: 'NIR16', 12: 'HRV'} ## band index in nat file

    ## determine subset once
    sub, lat, lon, vaa, vza = [None] * 5

    ## run through inputfiles
    ofiles = []
    for fi, bundle in enumerate(inputfile):
        file = ac.seviri.bundle_test(bundle)

        ## get nat  metadata
        meta = ac.seviri.metadata(file)

        ## platform
        platform = meta['ASTI']
        instrument = None
        if meta['AIID'] == 'SEVI': instrument = 'SEVIRI'
        sensor = '{}_{}'.format(platform, instrument)

        if ('SEVIRI' in sensor) & ('MSG' in sensor):
            band_idx = [1,2,3,12]
        else:
            print('Sensor {} not supported'.format(sensor))
            continue

        ## sub satellite longitude
        lon_0 = float(meta['LLOS'])

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        if output is None: output = setu['output']

        ## some assumptions here
        pt = meta['SNIT']
        isotime = '{}-{}-{}T{}:{}:{}'.format(pt[0:4], pt[4:6], pt[6:8], pt[8:10],pt[10:12],pt[12:],)
        dt = dateutil.parser.parse(isotime)
        #start_time = meta['SSBT']
        #end_time = meta['SSST']

        ## global attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L1R'

        ## output name
        oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        ofile_hrv = '{}/{}_{}_HRV.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile
        gatts['ofile_hrv'] = ofile_hrv

        ## read rsr
        rsrd = ac.shared.rsr_dict(sensor)[sensor]

        ## figure out location
        if sub is None: ## track subset across runs
            if setu['sub'] is not None:
                sub = setu['sub']
            elif setu['limit'] is not None:
                lon, lat = ac.seviri.geom(lon_0 = lon_0, geometry = False)
                sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
                if sub is not None:
                    print('Using sub = {}'.format(sub))
                    lat = lat[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    lon = lon[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                else:
                    print('Limit {} not in full disk.'.format(limit))
                    continue
            else:
                print('Warning: running without subsetting full disk image.')
                print('It is recommended using a subset sub=x, y, nx, ny or limit=S, W, N, E.')

            ## HRV sub
            sub_hrv = None if sub is None else [v * 3 for v in sub]

            ## read geolocation
            if (lon is None) | (lat is None):
                lon, lat = ac.seviri.geom(lon_0 = lon_0, geometry = False, sub = sub)
            ## read geometry
            vaa, vza = ac.seviri.geom(lon_0 = lon_0, geolocation = False, sub = sub)

            ## geom function already provides mask
            mask = np.isnan(lon)
        ## end read location

        ## sun position
        spos = ac.shared.sun_position(dt, lon, lat)
        sza = spos['zenith']
        saa = spos['azimuth']
        d = spos['distance']
        spos = None

        ## relative azimuth
        raa = np.abs(saa-vaa)
        raa[raa>180] = 360 - raa[raa>180]

        ## new gem handling
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## write position and angles
        gemo.write('lon', lon)
        gemo.write('lat', lat)
        gemo.write('vza', vza)
        gemo.write('vaa', vaa)
        gemo.write('sza', sza)
        gemo.write('saa', saa)
        gemo.write('raa', raa)

        ## additional mask for sun below horizon
        ## can be made an option?
        mask = mask | (sza > 90)

        ## run through bands
        for b in band_idx:
            band = band_names[b]
            ds_att = {'wavelength':rsrd['wave_nm'][band]}

            ## output HRV data
            if (band == 'HRV'):
                if (setu['seviri_hrv']):
                    if verbosity > 1: print('Converting bands: Reading {} from {}'.format(band, file))
                    data = ac.seviri.read_nat(file, b, sub=sub_hrv) ## HRV only as DN
                    ac.output.nc_write(ofile_hrv, 'hrv', data, new = True)
                continue

            ## output nominal resolution bands
            if verbosity > 1: print('Converting bands: Reading {} from {}'.format(band, file))

            ## read radiance data
            data = ac.seviri.read_nat(file, b, sub=sub, radiance = True)
            ## add mask for negatives (out of disk?)
            # data[data<0] = np.nan

            if setu['output_lt']:
                ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                gemo.write(ds, data, ds_att = ds_att)

            ## convert to reflectance
            data = (np.pi * data * d ** 2) / (f0s[sensor][band] * np.cos(np.radians(sza)))
            data[mask] = np.nan
            ds = 'rhot_{}'.format(rsrd['wave_name'][band])
            gemo.write(ds, data, ds_att = ds_att)
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
            data = None

        gemo.close()
        gemo = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
