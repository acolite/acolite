## def l1_convert
## converts EMIT L1B RAD data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-02-09
## modifications: 2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) check for rad/obs files, cleaned up settings use, output naming
##                2025-03-27 (QV) added subsetting support

def l1_convert(inputfile, output=None, settings = None):
    import netCDF4, os
    import dateutil.parser
    import numpy as np
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'ISS_EMIT'

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

        ## find obs data for geometry
        bd = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        if 'EMIT_L1B_OBS' in bn:
            obs_file = '{}/{}'.format(bd, bn)
            rad_file = '{}/{}'.format(bd, bn.replace('EMIT_L1B_OBS', 'EMIT_L1B_RAD'))
        if 'EMIT_L1B_RAD' in bn:
            obs_file = '{}/{}'.format(bd, bn.replace('EMIT_L1B_RAD', 'EMIT_L1B_OBS'))
            rad_file = '{}/{}'.format(bd, bn)
        if setu['verbosity'] > 1:
            print('Selected RAD file: {}'.format(rad_file))
            print('Selected OBS file: {}'.format(obs_file))

        if not os.path.exists(obs_file):
            print('OBS file missing: {}'.format(obs_file))
            continue
        if not os.path.exists(rad_file):
            print('RAD file missing: {}'.format(rad_file))
            continue

        ## open metadata
        meta = ac.shared.nc_gatts(rad_file)

        ## read location and band data
        sub = None
        with netCDF4.Dataset(rad_file) as ds:
            ## get band data
            dsb = ds['/sensor_band_parameters/']
            band_data = {k:dsb[k][:].data for k in dsb.variables.keys()}

            ## get location data
            dsl = ds['/location/']
            loc = {k:dsl[k][:].data for k in dsl.variables.keys() if k not in ['elev', 'glt_x', 'glt_y']}

            ## determine subset
            if setu['limit'] is not None:
                sub = ac.shared.geolocation_sub(loc['lat'], loc['lon'], setu['limit'])
                if sub is None:
                    print('Limit outside of scene {}'.format(file))
                    continue
                ## crop to sub
                for k in loc:
                    loc[k] = loc[k][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    print(k, loc[k].shape)
            #band_data = {k:dsb[k][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]].data for k in dsb.variables.keys()}

        ## read geometry
        with netCDF4.Dataset(obs_file) as ds:
            obs_datasets = ds['sensor_band_parameters']['observation_bands'][:]
            geom = {}
            if sub is None:
                geom['vaa'] = ds['obs'][:,:, 1]
                geom['vza'] = ds['obs'][:,:, 2]
                geom['saa'] = ds['obs'][:,:, 3]
                geom['sza'] = ds['obs'][:,:, 4]
            else:
                geom['vaa'] = ds['obs'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2], 1]
                geom['vza'] = ds['obs'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2], 2]
                geom['saa'] = ds['obs'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2], 3]
                geom['sza'] = ds['obs'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2], 4]

        ## use centre date/time
        dt0 = dateutil.parser.parse(meta['time_coverage_start'])
        dt1 = dateutil.parser.parse(meta['time_coverage_end'])
        dt = dt0 + (dt1-dt0)/2

        ## compute scene center sun position
        clon = np.nanmedian(loc['lon'])
        clat = np.nanmedian(loc['lat'])
        spos = ac.shared.sun_position(dt, clon, clat)
        #cossza = np.cos(np.radians(spos['zenith']))
        cossza = np.cos(np.radians(geom['sza']))
        d = spos['distance']

        ## output attributes
        gatts = {}
        gatts['satellite'] = meta['platform']
        gatts['instrument'] = meta['instrument']
        gatts['sensor'] = '{}_{}'.format(gatts['satellite'], gatts['instrument'])
        gatts['acolite_file_type'] = 'L1R'
        gatts['isodate'] = dt.isoformat()

        ## centre sun position
        gatts['doy'] = dt.strftime('%j')
        gatts['se_distance'] = spos['distance']

        #gatts['sza'] = spos['zenith'][0]
        #gatts['saa'] = spos['azimuth'][0]
        gatts['sza'] = np.nanmean(geom['sza'])
        gatts['saa'] = np.nanmean(geom['saa'])
        gatts['mus'] = np.cos(np.radians(gatts['sza']))
        gatts['vza'] = np.nanmean(geom['vza'])
        gatts['vaa'] = np.nanmean(geom['vaa'])

        gatts['raa'] = np.abs(gatts['saa']-gatts['vaa'])
        while gatts['raa'] > 180: gatts['raa'] = np.abs(gatts['raa']-360)

        ## create band dataset
        band_names = ['{}'.format(b) for b in range(0, len(band_data['wavelengths']))]
        rsr = {'{}'.format(b): ac.shared.gauss_response(band_data['wavelengths'][bi], band_data['fwhm'][bi], step=0.1) \
               for bi, b in enumerate(band_names)}
        band_rsr= {b: {'wave': rsr[b][0]/1000, 'response': rsr[b][1]} for b in rsr}
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)
        bands = {}
        for bi, b in enumerate(band_names):
            cwave = band_data['wavelengths'][bi]
            swave = '{:.0f}'.format(cwave)
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':swave, 'width': band_data['fwhm'][bi],
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
            gemo.write('lon', loc['lon'])
            if setu['verbosity']  > 1: print('Wrote lon ({})'.format(loc['lon'].shape))
            gemo.write('lat', loc['lat'])
            if setu['verbosity']  > 1: print('Wrote lat ({})'.format(loc['lat'].shape))

        if setu['output_geometry']:
            if setu['verbosity'] > 1: print('Writing geometry')
            for k in geom:
                gemo.write(k, geom[k])
                if setu['verbosity'] > 1: print('Wrote {} ({})'.format(k, geom[k].shape))

        ## read radiance data
        setu['hyper_read_cube'] = True
        if setu['hyper_read_cube']:
            if setu['verbosity'] > 1: print('Reading radiance cube')
            data, att = ac.shared.nc_data(rad_file, 'radiance', sub = sub, axis_3d = 2, attributes = True)
            print(data.shape)

        ## run through bands and store rhot
        for bi, b in enumerate(bands):
            ## get dataset attributes
            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## copy radiance
            if setu['hyper_read_cube']:
                cdata_radiance = data[:,:,bi]

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)
                print('Wrote Lt_{}'.format(bands[b]['wave_name']))

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] / 10 * cossza)
            cdata_radiance = None

            ## write toa reflectance
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))
        gemo.close()
        print('Wrote {}'.format(ofile))
        ofiles.append(ofile)
    return(ofiles, setu)
