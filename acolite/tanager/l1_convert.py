## def l1_convert
## converts Tanager L1B RAD data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2025-07-22
## modifications:

def l1_convert(inputfile, output = None, settings = None):
    import netCDF4, os
    import datetime
    import numpy as np
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'Tanager'

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    if output is None: output = setu['output']

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

    if type(inputfile) is not list: inputfile = [inputfile]


    ofiles = []
    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)

        ### find obs data for geometry
        #bd = os.path.dirname(bundle)
        #bn = os.path.basename(bundle)

        ## open metadata
        meta = ac.tanager.metadata(bundle)

        ## read location data
        with netCDF4.Dataset(bundle) as nc:
            lat = nc.groups['HDFEOS']['SWATHS']['HYP']['Geolocation Fields']['Latitude'][:]
            lon = nc.groups['HDFEOS']['SWATHS']['HYP']['Geolocation Fields']['Longitude'][:]

        ## find image subset
        sub = None
        if setu['sub'] is not None:
            sub = setu['sub']
        elif setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit {} not in current product {}.'.format(setu['limit'], bn))
                continue

        ## subset geolocation
        column_range, row_range = None, None
        if sub is not None:
            column_range = sub[0], sub[0]+sub[2]
            row_range = sub[1], sub[1]+sub[3]
            lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
            lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]

        ## read geometry
        with netCDF4.Dataset(bundle) as nc:
            if sub is None:
                vaa = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sensor_azimuth'][:]
                vza = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sensor_zenith'][:]
                saa = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sun_azimuth'][:]
                sza = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sun_zenith'][:]
            else:
                vaa = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sensor_azimuth'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                vza = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sensor_zenith'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                saa = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sun_azimuth'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                sza = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['sun_zenith'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        ## relative azimuth
        raa = np.abs(saa-vaa)
        raa[raa>180] = 360 - raa[raa>180]

        ## cosine sun angle
        mus = np.cos(np.radians(sza))

        ## use strip id date/time
        sp = meta['strip_id'].split('_')
        date_start, time_start = sp[0], sp[1]
        dts = [int(v) for v in [date_start[0:4], date_start[4:6], date_start[6:8], time_start[0:2], time_start[2:4], time_start[4:6]]]
        dt = datetime.datetime(dts[0], dts[1], dts[2], dts[3], dts[4], dts[5], tzinfo = datetime.timezone.utc)
        isodate = dt.isoformat()

        spos = ac.shared.sun_position(isodate, np.nanmean(lon), np.nanmean(lat))
        d = spos['distance']

        ## output attributes
        gatts = {}
        gatts['sensor'] = 'Tanager'
        gatts['acolite_file_type'] = 'L1R'
        gatts['isodate'] = dt.isoformat()

        ## centre sun position
        gatts['doy'] = dt.strftime('%j')
        gatts['se_distance'] = spos['distance']

        gatts['sza'] = np.nanmean(sza)
        gatts['saa'] = np.nanmean(saa)
        gatts['mus'] = np.cos(np.radians(gatts['sza']))
        gatts['vza'] = np.nanmean(vza)
        gatts['vaa'] = np.nanmean(vaa)

        gatts['raa'] = np.abs(gatts['saa']-gatts['vaa'])
        while gatts['raa'] > 180: gatts['raa'] = np.abs(gatts['raa']-360)

        ## create band dataset
        band_names = ['{}'.format(b) for b in range(0, len(meta['wavelengths']))]
        rsr = {'{}'.format(b): ac.shared.gauss_response(meta['wavelengths'][bi], meta['fwhm'][bi], step=0.1) \
               for bi, b in enumerate(band_names)}
        band_rsr= {b: {'wave': rsr[b][0]/1000, 'response': rsr[b][1]} for b in rsr}
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)
        bands = {}
        for bi, b in enumerate(band_names):
            cwave = meta['wavelengths'][bi]
            swave = '{:.0f}'.format(cwave)
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':swave, 'width': meta['fwhm'][bi],
                           'rsr': band_rsr[b],'f0': f0d[b]}

        ## store band info in gatts
        gatts['band_waves'] = [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = [bands[w]['width'] for w in bands]

        ## output file
        oname  = '{}_{}'.format(gatts['sensor'],  dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## outputfile
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## write lat/lon
        if (setu['output_geolocation']):
            if setu['verbosity']  > 1: print('Writing geolocation lon/lat')
            gemo.write('lon', lon)
            if setu['verbosity']  > 1: print('Wrote lon ({})'.format(lon.shape))
            del lon
            gemo.write('lat', lat)
            if setu['verbosity']  > 1: print('Wrote lat ({})'.format(lat.shape))
            del lat

        ## write geometry
        if setu['output_geometry']:
            if setu['verbosity'] > 1: print('Writing geometry')
            gemo.write('sza', sza)
            if setu['verbosity']  > 1: print('Wrote sza ({})'.format(sza.shape))
            del sza
            gemo.write('saa', saa)
            if setu['verbosity']  > 1: print('Wrote saa ({})'.format(saa.shape))
            del saa
            gemo.write('vza', vza)
            if setu['verbosity']  > 1: print('Wrote vza ({})'.format(vza.shape))
            del vza
            gemo.write('vaa',vaa)
            if setu['verbosity']  > 1: print('Wrote vaa ({})'.format(vaa.shape))
            del vaa
            gemo.write('raa', raa)
            if setu['verbosity']  > 1: print('Wrote raa ({})'.format(raa.shape))
            del raa

        ## read radiance data
        setu['hyper_read_cube'] = True
        if setu['hyper_read_cube']:
            with netCDF4.Dataset(bundle) as nc:
                if setu['verbosity'] > 1: print('Reading radiance cube')
                if sub is None:
                    rad = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['toa_radiance'][:]
                else:
                    rad = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['toa_radiance'][:, sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        ## run through bands and store rhot
        for bi, b in enumerate(bands):
            ## get dataset attributes
            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## copy radiance
            if setu['hyper_read_cube']:
                cdata_radiance = rad[bi, :,:]

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)
                print('Wrote Lt_{}'.format(bands[b]['wave_name']))

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * mus)
            cdata_radiance = None

            ## write toa reflectance
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))
        gemo.close()
        print('Wrote {}'.format(ofile))
        ofiles.append(ofile)
    return(ofiles, setu)
