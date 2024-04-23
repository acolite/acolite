## def l1_convert
## converts S2Resampling NC file to ACOLITE L1R file
## written by Quinten Vanhellemont, RBINS
## 2023-05-02
## modifications: 2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2024-04-23 (MB) read ancillary data of resampled input if requested

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
            odir = os.path.dirname(bundle)
        else:
            odir = output
        s2_auxiliary_default = setu['s2_auxiliary_default']

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

        ## auxiliary data, S2Resampling or msiresampling
        if s2_auxiliary_default:
            try:
                datasets = ac.shared.nc_datasets(bundle)
                # msiresampling can interpolate ancillary to each pixel
                if 'tcwv_interpolated' in datasets \
                        and 'msl_interpolated' in datasets \
                        and 'tco3_interpolated' in datasets \
                        and 'u10_interpolated' in datasets \
                        and 'v10_interpolated' in datasets:
                    for source,name,b in [('AUX_ECMWFT', 'tcwv', 'tcwv'),
                                          ('AUX_ECMWFT', 'msl', 'msl'),
                                          ('AUX_ECMWFT', 'tco3', 'tco3'),
                                          ('AUX_ECMWFT', 'u10', '_0u'),
                                          ('AUX_ECMWFT', 'v10', '_0v')]:
                        data = ac.shared.nc_data(bundle, f"{name}_interpolated")
                        height, width = data.shape
                        # ACOLITE so far uses a single ancillary value only and does not interpolate to the pixels.
                        # An array is expected in values, it is the central value that will be selected
                        gatts['{}_{}_{}'.format(source, b, 'dimensions')] = (3, 1)
                        gatts['{}_{}_{}'.format(source, b, 'values')] = [data[0, 0], data[height//2, width//2], data[-1, -1]]
                        gatts['{}_{}_{}'.format(source, b, 'longitudes')] = []
                        gatts['{}_{}_{}'.format(source, b, 'latitudes')] = []
                        del data
                    ac.shared.nc_gatts_update(ofile, gatts)
                    print(f"using embedded ECMWF ancillary data available on pixel resolution")
                else:
                    # S2Resampling in a patched version provides lat and lon for the ancillary data grid.
                    # A ticket is open to add this behaviour to SNAP.
                    # Until then s2_auxiliary_default must be used with s2_auxiliary_interpolate=False
                    try:
                        lat = ac.shared.nc_data(bundle, 'aux_latitude')
                        lon = ac.shared.nc_data(bundle, 'aux_longitude')
                    except Exception as e:
                        print("no aux_latitude and aux_longitude in tie-points. Use Calvalus S2 reader: {e}")
                        lat = []
                        lon = []
                    for source,b in [('AUX_ECMWFT', 'tcwv'),
                                     ('AUX_ECMWFT', 'msl'),
                                     ('AUX_ECMWFT', 'tco3'),
                                     ('AUX_ECMWFT', '_0u'),
                                     ('AUX_ECMWFT', '_0v')]:
                        data = ac.shared.nc_data(bundle, b)
                        gatts['{}_{}_{}'.format(source, b, 'dimensions')] = data.shape
                        gatts['{}_{}_{}'.format(source, b, 'values')] = data.flatten()
                        gatts['{}_{}_{}'.format(source, b, 'longitudes')] = lon.flatten()
                        gatts['{}_{}_{}'.format(source, b, 'latitudes')] = lat.flatten()
                    ac.shared.nc_gatts_update(ofile, gatts)
                    print(f"using embedded ECMWF ancillary data available on coarse geographic grid")
            except Exception as e:
                print(f"cannot use embedded ECMWF data: {e}")
                raise e

        ofiles.append(ofile)
    return(ofiles, setu)
