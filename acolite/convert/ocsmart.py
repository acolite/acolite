# def convert.ocsmart
# converts OCSMART he5 file to L2A ACOLITE NetCDF
# written by Quinten Vanhellemont, RBINS
# 2026-06-18
#
# modifications: 2026-06-22 (QV) moved from ocsmart.l1_convert

def ocsmart(inputfile, output = None, settings = None):
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
        if igatts != {}: continue

        ## find groups
        groups_dict = ac.shared.nc_groups(bundle)
        if groups_dict == {}: continue

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
            lon = ac.shared.nc_data(bundle, 'Longitude', group = geo_group)
            lat = ac.shared.nc_data(bundle, 'Latitude', group = geo_group)

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
        odict = {'Latitude': 'lat', 'Longitude': 'lon', 'L2_flags': 'l2_flags',
                'Relative_azimuth': 'raa', 'Sensor_zenith': 'vza', 'Solar_zenith': 'sza',
                'chlor_a(oci)': 'chlor_a_oci', 'chlor_a(yoc)': 'chlor_a_yoc', }

        ## run through groups
        for group in [None, 'AOD', 'Rrs', 'Lt']:
            datasets = ac.shared.nc_datasets(bundle, group = group)
            if datasets is None: continue

            if len(datasets) > 0:
                #print(group, datasets)
                for ds in datasets:
                    sp = ds.split('_')
                    d, att = ac.shared.nc_data(bundle, ds, group = group, sub = sub, attributes = True)
                    dso = '{}'.format(ds)
                    #if dso in ['L2_flags']: continue
                    if sp[0] in ['AOD', 'Rrs', 'Lt']:
                        wave = float(sp[1][:sp[1].find('nm')])
                        att['wavelength'] = wave
                        dso = '{}_{:.0f}'.format(sp[0], wave)
                    if ds in odict: dso = odict[ds]
                    print('Writing {} to {}'.format(ds, dso))

                    if att == {}: att = None
                    gemo.write(dso, d.data, ds_att = att)
            gemo.close()
            ofiles.append(ofile)
    return(ofiles, setu)
