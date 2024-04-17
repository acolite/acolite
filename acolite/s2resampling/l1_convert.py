## def l1_convert
## converts S2Resampling NC file to ACOLITE L1R file
## written by Quinten Vanhellemont, RBINS
## 2023-05-02
## modifications: 2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling

def l1_convert(inputfile, output = None, settings = {}, verbosity = 5):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac

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
        ret = ac.s2resampling.bundle_test(bundle)
        if (ret is None): continue
        sensor, gatts, datasets = ret

        ## sensor rsrd
        rsrd = ac.shared.rsr_dict(sensor)[sensor]

        ## merge sensor specific settings
        setu = ac.acolite.settings.parse(sensor, settings=settings)
        output = setu['output']
        if output is None:
            odir = os.path.dirname(imagefile)
        else:
            odir = output

#     ## check if ROI polygon is given
#     clip, clip_mask = False, None
#     if poly is not None:
#         if os.path.exists(poly):
#             try:
#                 limit = ac.shared.polygon_limit(poly)
#                 print('Using limit from polygon envelope: {}'.format(limit))
#                 clip = True
#             except:
#                 print('Failed to import polygon {}'.format(poly))

        ## datetime
        dt = dateutil.parser.parse(gatts['start_date'])
        ## output file name
        obase = '{}_{}_S2R_L1R.nc'.format(sensor, dt.strftime('%Y_%m_%d_%H_%M_%S'))
        ofile = '{}/{}'.format(odir, obase)

        ## global attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L1R'
        gatts['obase'] = obase

        ## rename datasets to acolite L1R
        dsets = {'lat': 'lat', 'lon': 'lon',
                  'view_zenith_mean': 'vza', 'view_azimuth_mean': 'vaa',
                  'sun_zenith': 'sza', 'sun_azimuth': 'saa'}

        ## read data
        gemi = ac.gem.gem(bundle)
        data = {}
        for ds in dsets:
            print('Reading {}'.format(ds))
            d, att = gemi.data(ds, attributes = True)
            data[dsets[ds]] = d

        ## compute raa
        data['raa'] = np.abs(data['saa'] - data['vaa'])
        raasub = np.where(data['raa'] > 180)
        data['raa'][raasub] =np.abs(360-data['raa'][raasub])

        ## write data
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        for dso in data:
            print('Writing {}'.format(dso))
            gemo.write(dso, data[dso])
            data[dso] = None

        saa = None
        ## read and write band data
        for bi, b in enumerate(rsrd['rsr_bands']):
            ds = 'B{}'.format(b)
            dso = 'rhot_{}'.format(rsrd['wave_name'][b])
            print(ds, dso)
            d, att = gemi.data(ds, attributes = True)
            gemo.write(dso, d)

            ## add band specific geometry data
            if setu['geometry_per_band']:
                ## read sun azimuth angle
                if saa is None: saa = gemi.data('sun_azimuth')

                ## read band view zenith angle
                dso = 'vza_{}'.format(rsrd['wave_name'][b])
                print(dso)
                d, att = gemi.data('view_zenith_B{}'.format(b), attributes = True)
                gemo.write(dso, d)

                ## read band view azimuth angle
                dso = 'vaa_{}'.format(rsrd['wave_name'][b])
                print(dso)
                d, att = gemi.data('view_azimuth_B{}'.format(b), attributes = True)
                gemo.write(dso, d)

                dso = 'raa_{}'.format(rsrd['wave_name'][b])
                d = np.abs(saa - d)
                raasub = np.where(d > 180)
                d[raasub] = np.abs(360 - d[raasub])
                print(dso)
                gemo.write(dso, d)
        gemi.close()
        gemo.close()
        ofiles.append(ofile)
    return(ofiles, setu)
