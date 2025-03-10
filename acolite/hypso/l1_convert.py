## def l1_convert
## converts HYPSO file to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-09
## modifications: 2023-03-27 (QV) fixed issue with band naming
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming
##                2025-03-10 (QV) added HYPSO-2 support, added subsetting

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac
    import h5py

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']
    if output is None: output = setu['output']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ofiles = []
    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)

        gatts = ac.shared.nc_gatts(bundle)
        if gatts['instrument'] == 'HYPSO-1 Hyperspectral Imager':
            sensor = 'HYPSO1'
        elif gatts['instrument'] == 'HYPSO-2 Hyperspectral Imager':
            sensor = 'HYPSO2'
        else:
            print('HYPSO processing of sensor {} not configured'.format(gatts['instrument']))
            continue

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

        f = h5py.File(bundle, mode='r')

        ## get band information
        if 'Lt' in f['/products']:
            lt_pars = 'Lt'
            attributes = {k: f['/products']['Lt'].attrs[k] for k in f['/products']['Lt'].attrs.keys() if k != 'DIMENSION_LIST'}
            waves = attributes['wavelengths']
            fwhm = attributes['fwhm']
        else:
            lt_pars = []
            waves = []
            fwhm = []
            for ds in f['/products']:
                if ds.startswith('Lt_'):
                    lt_pars.append(ds)
                    att = {k: f['/products'][ds].attrs[k] for k in f['/products'][ds].attrs.keys() if k != 'DIMENSION_LIST'}
                    waves.append(att['wavelength'][0])
                    fwhm.append(att['fwhm'][0])

        rsr = ac.shared.rsr_hyper(waves, fwhm, step=0.1)
        rsrd = ac.shared.rsr_dict(rsrd={sensor:{'rsr':rsr}})

        ## create band dataset
        band_names = ['{}'.format(b) for b in range(0, len(waves))]
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsr)
        bands = {}

        for bi, b in enumerate(band_names):
            cwave = rsrd[sensor]['wave_nm'][bi]
            swave = rsrd[sensor]['wave_name'][bi]
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':swave, 'width': fwhm[bi],
                           'rsr': rsr[bi],'f0': f0d[bi]}

        ## time
        if 'unixtime' in f['/navigation/']:
            utime = f['/navigation/']['unixtime'][:]
            start_time = datetime.datetime.utcfromtimestamp(utime[0])
            stop_time = datetime.datetime.utcfromtimestamp(utime[-1])
            dt = start_time + (stop_time-start_time)
        else:
            dt = dateutil.parser.parse(gatts['timestamp_acquired'].strip('Z')) ## timestamp format is wrong, has both +00:00 and Z

        ## read lat and lon
        lat = f['/navigation/']['latitude'][:]
        lon = f['/navigation/']['longitude'][:]

        sub = None
        if setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit outside of scene {}'.format(bundle))
                continue
            lat = lat[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
            lon = lon[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]

        ## geometry
        if sub is None:
            vza = f['/navigation/']['sensor_zenith'][:]
            vaa = f['/navigation/']['sensor_azimuth'][:]
            sza = f['/navigation/']['solar_zenith'][:]
            saa = f['/navigation/']['solar_azimuth'][:]
        else:
            vza = f['/navigation/']['sensor_zenith'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
            vaa = f['/navigation/']['sensor_azimuth'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
            sza = f['/navigation/']['solar_zenith'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
            saa = f['/navigation/']['solar_azimuth'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]

        ## compute relative azimuth angle
        raa = np.abs(saa-vaa)
        ## raa along 180 degree symmetry
        tmp = np.where(raa>180)
        raa[tmp]=np.abs(raa[tmp] - 360)


        ## date time
        isodate = dt.isoformat()

        ## compute scene center sun position
        clon = np.nanmedian(lon)
        clat = np.nanmedian(lat)
        cspos = ac.shared.sun_position(dt, clon, clat)
        #cossza = np.cos(np.radians(spos['zenith']))

        spos = ac.shared.sun_position(dt, lon, lat)
        d = spos['distance']
        #sza = spos['zenith']
        #saa = spos['azimuth']
        spos = None

        ## output attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['acolite_file_type'] = 'L1R'
        gatts['isodate'] = dt.isoformat()

        ## centre sun position
        gatts['doy'] = dt.strftime('%j')
        gatts['se_distance'] = cspos['distance']

        gatts['sza'] = np.nanmean(sza)
        gatts['saa'] = np.nanmean(saa)
        gatts['mus'] = np.cos(np.radians(gatts['sza']))
        gatts['vza'] = np.nanmean(vza)
        gatts['vaa'] = np.nanmean(vaa)

        gatts['raa'] = np.abs(gatts['saa']-gatts['vaa'])
        while gatts['raa'] > 180: gatts['raa'] = np.abs(gatts['raa']-360)

        ## store band info in gatts
        gatts['band_waves'] = waves # [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = fwhm # [bands[w]['width'] for w in bands]

        ## output file
        oname  = '{}_{}'.format(gatts['sensor'],  dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## output file
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## write lat/lon
        if (setu['output_geolocation']):
            if verbosity > 1: print('Writing geolocation lon/lat')
            gemo.write('lon', lon)
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
            lon = None
            gemo.write('lat', lat)
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
            lat = None

        ## write geometry
        if (setu['output_geometry']):
            if verbosity > 1: print('Writing geometry')
            gemo.write('vza', vza)
            if verbosity > 1: print('Wrote vza ({})'.format(vza.shape))
            vza = None
            gemo.write('vaa', vaa)
            if verbosity > 1: print('Wrote vaa ({})'.format(vaa.shape))
            vaa = None
            gemo.write('sza', sza)
            if verbosity > 1: print('Wrote sza ({})'.format(sza.shape))
            gemo.write('saa', saa)
            if verbosity > 1: print('Wrote saa ({})'.format(saa.shape))
            saa = None
            gemo.write('raa', raa)
            if verbosity > 1: print('Wrote raa ({})'.format(raa.shape))
            raa = None

        ## cosine of sun angle
        cossza = np.cos(np.radians(sza))
        sza = None

        ## read data cube
        if (lt_pars == 'Lt') & (setu['hyper_read_cube']):
            if sub is None:
                data = f['/products']['Lt'][:]
            else:
                data= f['/products']['Lt'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2], :]
            nbands = data.shape[2]
            data_dimensions = data.shape[0], data.shape[1]

        ## run through bands and store rhot
        for bi, b in enumerate(bands):
            ## get dataset attributes
            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## copy radiance
            if (lt_pars == 'Lt'):
                if (setu['hyper_read_cube']):
                    cdata_radiance = data[:,:,bi]
                else:
                    if sub is None:
                        cdata_radiance = f['/products']['Lt'][:, :, bi]
                    else:
                        cdata_radiance = f['/products']['Lt'][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2], bi]
            else:
                if sub is None:
                    cdata_radiance = f['/products'][lt_pars[bi]][:]
                else:
                    cdata_radiance = f['/products'][lt_pars[bi]][sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)
                print('Wrote Lt_{}'.format(bands[b]['wave_name']))

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
            cdata_radiance = None

            ## write toa reflectance
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))
        data = None

        ## close file
        if f is not None:
            f.close()
            f = None

        gemo.close()
        ofiles.append(ofile)
    return(ofiles, setu)
