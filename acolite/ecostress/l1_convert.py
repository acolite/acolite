## def l1_convert
## converts ECOSTRESS L1B data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2022-08-11
## modifications: 2022-08-19 (QV) added ECO1BRAD support
##                2024-01-17 (QV) added geofile keyword
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output=None, settings = None):
    import os, h5py, json
    import numpy as np
    import dateutil.parser, datetime
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'ISS_ECOSTRESS'

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    verbosity = setu['verbosity']
    if output is None: output = setu['output']

    ## read rsr
    rsrd = ac.shared.rsr_dict(sensor, wave_range=[7,14], wave_step=0.05)[sensor]

    sub = None
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

        dn = os.path.dirname(bundle)
        bn, ext = os.path.splitext(os.path.basename(bundle))
        ## read metadata
        meta = ac.ecostress.attributes(bundle)

        sdate = dateutil.parser.parse(meta['RangeBeginningDate']+'T'+meta['RangeBeginningTime'])
        edate = dateutil.parser.parse(meta['RangeEndingDate']+'T'+meta['RangeEndingTime'])
        tdiff = (edate - sdate)
        dt = sdate + datetime.timedelta(days=tdiff.days/2, seconds=tdiff.seconds/2)
        isodate = dt.isoformat()[0:19]
        if 'Z' not in isodate: isodate+='Z'

        use_bt_lut = False
        ## read BT LUT
        if use_bt_lut: btlut = ac.ecostress.bt_lut()

        ## read thermal coefficients
        tcfile = ac.config['data_dir']+'/ECOSTRESS/ECOSTRESS_thermal_coefficients.json'
        coeffs = json.load(open(tcfile, 'r'))

        ## set up global attributes
        gatts = {'sensor': sensor, 'isodate': isodate,
                 'acolite_file_type': 'L1R',
                 'thermal_sensor': sensor, 'thermal_bands': ['1', '2', '3', '4', '5']}
        for bk in coeffs['thermal_coefficients']:
            for k in coeffs['thermal_coefficients'][bk]:
                gatts['{}_CONSTANT_BAND_{}'.format(k, bk)] = coeffs['thermal_coefficients'][bk][k]

        ## output file name
        oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        f = h5py.File(bundle, mode='r')

        ## find out file type
        file_type = None
        geo_file, geo_meta = None, None
        if 'Mapped' in f:
            file_type = 'ECO1BMAPRAD'
            data_key = 'Mapped'
            geo_key = 'Mapped'
            fg = f
        elif 'Radiance' in f:
            file_type = 'ECO1BRAD'
            data_key = 'Radiance'
            geo_key = 'Geolocation'
            if ac.settings['run']['geofile'] is not None: ## use provided geo file
                geo_file = '{}'.format(ac.settings['run']['geofile'])
            else: ## find geo file in local directory
                geo_file = dn + os.path.sep + bn.replace('ECOSTRESS_L1B_RAD_', 'ECOSTRESS_L1B_GEO_') + ext
            if not os.path.exists(geo_file):
                print('ECO1BGEO file required for processing ECO1BRAD.')
                print('{} not found'.format(geo_file))
                continue
            if ac.settings['run']['verbosity'] > 2: print('Using geo file {}'.format(geo_file))
            geo_meta = ac.ecostress.attributes(geo_file)
            fg = h5py.File(geo_file, mode='r')

        datasets = {'latitude': 'lat', 'longitude':'lon',
                    'height': 'altitude', 'solar_azimuth':'saa', 'solar_zenith':'sza',
                    'view_azimuth':'vaa', 'view_zenith':'vza',
                   }

        ## find crop
        if setu['limit'] is not None:
            lat = fg[geo_key]['latitude'][()]
            lon = fg[geo_key]['longitude'][()]

            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            lat = None
            lon = None
            if sub is None:
                print('Image {} does not cover {}'.format(bundle, setu['limit']))
                continue

        ## output file
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = gatts

        for ds in datasets:
            if sub is None:
                data = fg[geo_key][ds][()]
            else:
                data = fg[geo_key][ds][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

            if len(np.where(np.isfinite(data))[0]) == 0:
                print('Limit {} in blackfill of {}'.format(setu['limit'], bundle))
                continue

            ds_att = {'name': ds}
            print('Writing {}'.format(datasets[ds]))
            gemo.write(datasets[ds], data, ds_att = ds_att)

        ## run through bands
        for b in range(1, 6):
            bk = '{}'.format(b)
            ds = 'Lt{}'.format(b)
            if sub is None:
                Lt = f[data_key]['radiance_{}'.format(b)][()]
            else:
                Lt = f[data_key]['radiance_{}'.format(b)][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

            ds_att = {'name': 'radiance_{}'.format(b)}
            for k in rsrd:
                if 'rsr' in k: continue
                if bk in rsrd[k].keys():
                    ds_att[k] = rsrd[k][bk]
            for k in coeffs['thermal_coefficients'][bk]:
                ds_att[k] = coeffs['thermal_coefficients'][bk][k]

            if len(np.where(Lt > 0)[0]) == 0:
                print('Skipping B{}'.format(bk))
                continue

            if use_bt_lut:
                flat = Lt.flatten()
                flat = np.where(np.logical_and(0 <= flat, flat <= 60), flat, 0)  # Filter fill values
                idx = np.int32(flat / 0.001)

                # interpolate LUT data and compute BT
                # from https://git.earthdata.nasa.gov/projects/LPDUR/repos/ecostress_swath2grid (accessed 2022-08-11)
                Lt_x0 = idx * 0.001
                Lt_x1 = Lt_x0 + 0.001
                factor0 = (Lt_x1 - flat) / 0.001
                factor1 = (flat - Lt_x0) / 0.001
                bt = (factor0 * btlut[bk][idx]) + (factor1 * btlut[bk][idx + 1])
                bt = np.where((flat != 0), bt, np.nan)
                bt = bt.reshape(Lt.shape)
            else:
                bt = ds_att['K2']/(np.log(ds_att['K1']/Lt)+1)

            Lt[Lt<0] = np.nan
            bt[bt<0] = np.nan

            print('Writing {}'.format('lt{}'.format(b)))
            gemo.write('lt{}'.format(b), Lt, ds_att = ds_att)
            Lt = None

            print('Writing {}'.format('bt{}'.format(b)))
            gemo.write('bt{}'.format(b), bt, ds_att = ds_att)
            bt = None
        f, fg = None, None
        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)
        gemo.close()

    return(ofiles, setu)
