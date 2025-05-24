## def l1_convert
## convert EarthCare bundle to ACOLITE L1R
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-13
## modifications: 2025-05-22 (QV) converted to ACOLITE l1_convert function

def l1_convert(inputfile, output = None, settings = None):
    import os, datetime, dateutil.parser, time
    import numpy as np
    from netCDF4 import Dataset
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
            inputfile = inputfile.split(';') ## use semicolon as EUMETSAT uses commas in their file names
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## run through inputfiles
    ofiles = []
    for fi, bundle in enumerate(inputfile):
        ## test bundle and get igatts
        sensor, igatts = ac.earthcare.bundle_test(bundle)

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults
        if output is None: output = setu['output']

        ## read rsr
        if setu['rsr_version'] is not None:
            if setu['rsr_version'] == 'extracted':
                print('Warning: Using extracted RSR for {}'.format(sensor))
                sensor_version = '{}_{}'.format(sensor, setu['rsr_version'])
            else:
                sensor_version = '{}'.format(sensor)
        rsrd = ac.shared.rsr_dict(sensor_version)[sensor_version]

        ## get image time
        st = igatts['sensingStartTime']
        if st.startswith('UTC='): st = st[4:] + '+00:00'
        et = igatts['sensingStopTime']
        if et.startswith('UTC='): et = et[4:] + '+00:00'

        dt0 = dateutil.parser.parse(st)
        dt1 = dateutil.parser.parse(et)
        dt = dt0 + (dt1-dt0) ## should be weigthed depending on subset
        isodate = dt.isoformat()

        ## global attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L1R'

        ## output name
        oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## available bands
        bands = ac.shared.nc_data(bundle, 'band', group = 'ScienceData')

        ## read lat/lon
        lat = ac.shared.nc_data(bundle, 'latitude', group = 'ScienceData')
        lon = ac.shared.nc_data(bundle, 'longitude', group = 'ScienceData')

        ## find image subset
        sub = None
        if setu['sub'] is not None:
            sub = setu['sub']
        elif setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit {} not in bundle {}.'.format(limit, bundle))
                continue

        ## subset lat lon
        if sub is not None:
            column_range = sub[0], sub[0]+sub[2]
            row_range = sub[1], sub[1]+sub[3]
            ## subset geolocation
            lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
            lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]

        ## read geometry
        saa = ac.shared.nc_data(bundle, 'solar_azimuth_angle', group = 'ScienceData', sub = sub)
        sza = 90 - ac.shared.nc_data(bundle, 'solar_elevation_angle', group = 'ScienceData', sub = sub)
        vaa = ac.shared.nc_data(bundle, 'sensor_azimuth_angle', group = 'ScienceData', sub = sub)
        vza = 90 - ac.shared.nc_data(bundle, 'sensor_elevation_angle', group = 'ScienceData', sub = sub)

        ## sun position
        spos = ac.shared.sun_position(dt, np.nanmean(lon), np.nanmean(lat))
        se_distance = spos['distance']

        ## compute cosine solar zenith angle
        mus = np.cos(np.radians(sza))

        ## relative azimuth
        raa = np.abs(saa-vaa)
        raa[raa>180] = 360 - raa[raa>180]


        ## write geometry and geolocation data
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## write position and angles
        gemo.write('lon', lon)
        del lon
        gemo.write('lat', lat)
        del lat
        gemo.write('vza', vza)
        del vza
        gemo.write('vaa', vaa)
        del vaa
        gemo.write('sza', sza)
        del sza
        gemo.write('saa', saa)
        del saa
        gemo.write('raa', raa)
        del raa

        ## read image data
        nc = Dataset(bundle)
        for bi, band in enumerate(bands):
            ## convert to radiance
            radiance = False
            if band in rsrd['rsr_bands']:
                radiance = True
                if sub is not None:
                    irr = nc.groups['ScienceData']['solar_spectral_irradiance'][bi, column_range[0]:column_range[1]]
                else:
                    irr = nc.groups['ScienceData']['solar_spectral_irradiance'][bi, :]
            if (not radiance) & (not setu['output_bt']): continue

            ## read data
            if sub is not None:
                data = nc.groups['ScienceData']['pixel_values'][bi, row_range[0]:row_range[1], column_range[0]:column_range[1]]
            else:
                data = nc.groups['ScienceData']['pixel_values'][bi, :, :]


            ## convert to reflectance
            if radiance:
                ds_att = {'wavelength': rsrd['wave_nm'][band]}

                if setu['output_lt']:
                    ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                    gemo.write(ds, data, ds_att = ds_att)

                ## convert to reflectance
                data *= (np.pi * se_distance ** 2) / (mus)
                data /= irr

                ds = 'rhot_{}'.format(rsrd['wave_name'][band])
                gemo.write(ds, data, ds_att = ds_att)
            else:
                ds_att = {}
                ds = 'bt_{}'.format(band)
                gemo.write(ds, data, ds_att = ds_att)
            del data
        nc = None

        gemo.close()
        gemo = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
