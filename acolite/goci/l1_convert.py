## def l1_convert
## convert GOCI2 image to L1R ACOLITE NetCDF file
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-16
## modifications: 2025-07-17 (QV) added SVC, added segment to oname
##                2025-07-24 (QV) added sensor degradation coefficients

def l1_convert(inputfile, output = None, settings = None):
    import os, datetime
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
            inputfile = inputfile.split(';') ## use semicolon as EUMETSAT uses commas in their file names
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## run through inputfiles
    ofiles = []
    for bundle in inputfile:
        gatts = ac.shared.nc_gatts(bundle)
        if (gatts['platform'] != 'GK-2B') & (gatts['instrument'] != 'GOCI-II'):
            continue
        if gatts['title'] not in ['GK2B GOCI-II Level-1B Radiances']:
            continue

        bn = os.path.basename(bundle)
        segment = bn[bn.find('LA_')+3:-3]

        ## get sensor defaults
        if (gatts['platform'] == 'GK-2B') & (gatts['instrument'] == 'GOCI-II'):
            sensor = 'GK2_GOCI2'
        else:
            continue

        rsrd = ac.shared.rsr_dict(sensor)[sensor]
        bands = rsrd['rsr_bands']
        l_datasets = ['L_TOA_380', 'L_TOA_412', 'L_TOA_443', 'L_TOA_490', 'L_TOA_510', 'L_TOA_555',
                      'L_TOA_620', 'L_TOA_660', 'L_TOA_680', 'L_TOA_709', 'L_TOA_745', 'L_TOA_865']

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults
        if output is None: output = setu['output']

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

        ## compute date time
        date_start, time_start =  gatts['observation_start_time'].split('_')
        dts = [int(v) for v in [date_start[0:4], date_start[4:6], date_start[6:8], time_start[0:2], time_start[2:4], time_start[4:6]]]
        date_end, time_end = gatts['observation_end_time'].split('_')
        dte = [int(v) for v in [date_end[0:4], date_end[4:6], date_end[6:8], time_end[0:2], time_end[2:4], time_end[4:6]]]
        start_dt = datetime.datetime(dts[0], dts[1], dts[2], dts[3], dts[4], dts[5], tzinfo = datetime.timezone.utc)
        end_dt = datetime.datetime(dte[0], dte[1], dte[2], dte[3], dte[4], dte[5], tzinfo = datetime.timezone.utc)
        dt = start_dt + (end_dt - start_dt)
        isodate = dt.isoformat()

        ## read sensor degradation
        sensor_degradation = ac.goci.goci2_degradation(isodate)

        ## read lat/lon
        lat = ac.shared.nc_data(bundle, 'latitude', group = 'navigation_data')
        lon = ac.shared.nc_data(bundle, 'longitude', group = 'navigation_data')

        ## find image subset
        sub = None
        if setu['sub'] is not None:
            sub = setu['sub']
        elif setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit {} not in current product {}.'.format(setu['limit'], bundle))
                continue

        ## subset geolocation
        column_range, row_range = None, None
        if sub is not None:
            column_range = sub[0], sub[0]+sub[2]
            row_range = sub[1], sub[1]+sub[3]
            lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
            lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]

        ## compute viewing geometry
        vaa, vza = ac.goes.geometry(lon, lat, lon_0 = gatts['sub_longitude'],
                    r = gatts['nominal_satellite_height'],
                    a = gatts['semi_major_axis'],
                    b = gatts['semi_minor_axis'],
                    ealt = 0,)

        ## compute sun geometry
        print('Computing sun position for {}'.format(dt.isoformat()))
        spos = ac.shared.sun_position(isodate, lon, lat)
        se_distance = spos['distance']
        sza, saa = spos['zenith'], spos['azimuth']
        del spos

        sza_mask = sza > 90
        ## cosine sun zenith angle
        mus = np.cos(np.radians(sza))
        print('Cosine sun zenith angle shape: {}'.format(mus.shape))

        ## relative azimuth
        raa = np.abs(saa-vaa)
        raa[raa>180] = 360 - raa[raa>180]

        ## set up gatts
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L1R'

        ## output name
        oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if len(segment) > 0: oname += '_{}'.format(segment)
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## new gem handling
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

        ## run through bands
        for bi, band in enumerate(bands):
            ds = l_datasets[bi]
            d, a = ac.shared.nc_data(bundle, ds, group = 'geophysical_data', sub = sub, attributes = True)

            ## output dataset attributes
            ds_att = {'wavelength':rsrd['wave_nm'][band], 'f0': f0d[band]}

            ## apply sensor degradation factors
            if setu['goci2_sensor_degradation']:
                sd = sensor_degradation[bi]
                ds_att['sensor_degradation'] = sd
                if setu['verbosity'] > 2: print('Applying sensor degradation factor {:.5f} for {}'.format(sd, ds))
                d *= sd

            ## apply gains
            if (setu['gains']):
                if len(setu['gains_toa']) == len(bands):
                    cg = float(setu['gains_toa'][bi])
                    ds_att['gain'] = cg
                    if setu['verbosity'] > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                    d *= cg

            ## output radiance
            if setu['output_lt']:
                ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                gemo.write(ds, d, ds_att = a)

            ## convert to reflectance
            scale = (np.pi * se_distance ** 2) / (f0d[band] * mus)
            d *= scale

            ## mask negative TOA
            d[d <= 0] = np.nan

            ## write dataset
            ds = 'rhot_{}'.format(rsrd['wave_name'][band])
            gemo.write(ds, d, ds_att = ds_att)
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
            d = None

        ## close output file
        gemo.close()
        gemo = None

        if ofile is not None:
            if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
