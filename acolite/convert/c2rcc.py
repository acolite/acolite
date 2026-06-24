# def convert.c2rcc
# converts C2RCC nc file to L2A ACOLITE NetCDF
# written by Quinten Vanhellemont, RBINS
# 2026-06-23
#
# modifications:

def c2rcc(inputfile, output = None, settings = None, Rrs = True):
    import numpy as np
    import datetime, dateutil.parser, os
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

    ## run through inputfiles
    ofiles = []
    for fi, bundle in enumerate(inputfile):
        ## at the moment the global atts are empty for OCSMART V2.2
        igatts = ac.shared.nc_gatts(bundle)
        if not ('C2RCC' in igatts['product_type']): continue

        ## determine sensor and date from file name
        bn = os.path.basename(bundle)
        sp = bn.split('_')
        if 'MSIL1C' in bn:
            dt = dateutil.parser.parse(sp[2]+'Z')
            isodate = dt.isoformat()
            sensor = bn[0:7]
        else:
            print('{} not recognised'.format(bn))
            continue

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        ## do not convert L2 by default
        if (not setu['convert_l2']):
            print('To convert L2 files, set convert_l2 = True')
            continue

        if output is None: output = setu['output']

        sub = setu['sub']
        if setu['limit'] is not None:
            geo_group = None
            ## read lat and lon
            lon = ac.shared.nc_data(bundle, 'lon', group = geo_group)
            lat = ac.shared.nc_data(bundle, 'lat', group = geo_group)

            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])

            if (sub is None):
                print('Limit not in scene {}'.format(bundle))
                continue
            del lat, lon

        ## set up global attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L2A'

        ## output file name
        oname = '{}_{}'.format(sensor, dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## write data
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## dataset aliases
        odict = {}
        datasets_skip = [ 'metadata', 'crs',
                          'unc_apig', 'unc_adet', 'unc_agelb', 'unc_bpart', 'unc_bwit', 'unc_adg', 'unc_atot', 'unc_btot', 'unc_tsm', 'unc_chl',
                          'tcwv', 'tco3','msl','r','_0u','_0v','aux_latitude','aux_longitude',
                          'aod550','z','bcaod550', 'omaod550','ssaod550','suaod550','aod469','aod670','aod865','aod1240',
                          'Rtosa_OOS_mask', 'Rtosa_OOR_mask', 'Rhow_OOR_mask','Cloud_risk_mask', 'Iop_OOR_mask',
                          'Apig_at_max_mask', 'Adet_at_max_mask', 'Agelb_at_max_mask', 'Bpart_at_max_mask',
                          'Bwit_at_max_mask', 'Apig_at_min_mask', 'Adet_at_min_mask', 'Agelb_at_min_mask',
                          'Bpart_at_min_mask', 'Bwit_at_min_mask', 'Rhow_OOS_mask', 'Kd489_OOR_mask',
                          'Kdmin_OOR_mask', 'Kd489_at_max_mask', 'Kdmin_at_max_mask', 'Valid_PE_mask',
                           #'c2rcc_flags','lat','lon',
                        ]
                        
        ## read RSR with version
        sensor_version = '{}'.format(sensor)
        if setu['rsr_version'] is not None:
            sensor_version = '{}_{}'.format(sensor, setu['rsr_version'])
        rsrd = ac.shared.rsr_dict(sensor = sensor_version)[sensor_version]

        ## run through groups
        for group in [None]:
            datasets = ac.shared.nc_datasets(bundle, group = group)
            if datasets is None: continue

            if len(datasets) > 0:
                for ds in datasets:
                    if ds in datasets_skip: continue
                    d, att = ac.shared.nc_data(bundle, ds, group = group, sub = sub, attributes = True)
                    dso = '{}'.format(ds)
                    sp = ds.split('_')
                    if sp[0] == 'rrs':
                        band = sp[1][1:]
                        wave = rsrd['wave_nm'][band]
                        att['wavelength'] = wave
                        if not Rrs:
                            d *= np.pi
                            dso = '{}_{:.0f}'.format('rhow', wave)
                        else:
                            dso = '{}_{:.0f}'.format('Rrs', wave)

                    if ds in odict: dso = odict[ds]
                    print('Writing {} to {}'.format(ds, dso))
                    if att == {}: att = None
                    gemo.write(dso, d.data, ds_att = att)
            gemo.close()
            ofiles.append(ofile)
    return(ofiles, setu)
