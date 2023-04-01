## def l1_convert
## converts VIIRS data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications:

def l1_convert(inputfile, output = None, settings = {}, verbosity = 0):
    import h5py
    import numpy as np
    import scipy.ndimage
    import dateutil.parser, time
    import glob, os
    import acolite as ac

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ofile = None
    ofiles = []

    for bundle in inputfile:
        new = True

        files = ac.viirs.bundle_test(bundle)
        l1b_meta = None
        for mkey in ['l1b_mod', 'l1b_img']:
            if mkey in files:
                l1b_meta = ac.shared.nc_gatts(files[mkey])
                break
        l1b_sensor = '{}_{}'.format(l1b_meta['platform'], l1b_meta['instrument']).upper()
        print(l1b_sensor, l1b_meta['title'], l1b_meta['processing_level'])

        ## parse settings
        sensor = '{}_{}'.format(l1b_meta['platform'], l1b_meta['instrument']).upper()
        setu = ac.acolite.settings.parse(sensor, settings=settings)
        opts = setu['viirs_option'].lower().split('+')
        ## configure scan lines
        if opts[0] == 'mod': slines = 16
        if opts[0] == 'img': slines = 32

        rotate = True

        skip = False
        for vm in ['mod', 'img']:
            if vm in opts:
                for ft in ['l1b', 'geo']:
                    if '{}_{}'.format(ft, vm) not in files:
                        print('{}_{} not available for VIIRS {} processing of {}'.format(ft, vm, setu['viirs_option'], bundle))
                        skip = True
        if skip: continue

        limit = setu['limit']
        poly = setu['polygon']
        if poly is not None:
            if os.path.exists(poly):
                try:
                    limit = ac.shared.polygon_limit(poly)
                    if setu['polygon_limit']:
                        print('Using limit from polygon envelope: {}'.format(limit))
                    else:
                        limit = setu['limit']
                    clip = True
                except:
                    print('Failed to import polygon {}'.format(poly))

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        output_lt = setu['output_lt']
        vname = setu['region_name']
        odir = '{}'.format(output)

        ## load rsr
        rsrd = ac.shared.rsr_dict(sensor=sensor)

        new = True
        ## run through viirs options
        for vi, viirs_res in enumerate(opts):
            print(vi, viirs_res)

            l1bt = 'l1b_{}'.format(viirs_res)
            geot = 'geo_{}'.format(viirs_res)

            if l1bt not in files:
                print('{} file not found for {}'.format(l1bt, bundle))
                continue
            if geot not in files:
                print('{} file not found for {}'.format(geot, bundle))
                continue

            l1b = files[l1bt]
            geo = files[geot]

            ## get metadata
            l1b_meta = ac.shared.nc_gatts(l1b)
            l1b_sensor = '{}_{}'.format(l1b_meta['platform'], l1b_meta['instrument']).upper()
            print(l1b_meta['title'], l1b_meta['processing_level'], l1b_sensor)
            geo_meta = ac.shared.nc_gatts(geo)
            geo_sensor = '{}_{}'.format(geo_meta['platform'], geo_meta['instrument']).upper()
            print(geo_meta['title'], geo_meta['processing_level'], geo_sensor)

            if new:
                ## set up attributes
                sisotime = l1b_meta['time_coverage_start']
                eisotime = l1b_meta['time_coverage_end']
                stime = dateutil.parser.parse(sisotime)
                etime = dateutil.parser.parse(eisotime)
                td = etime - stime
                dtime = stime + td/2


                gatts = {}
                gatts['sensor'] = sensor
                gatts['isodate'] = dtime.isoformat()
                gatts['obase']  = '{}_{}_{}_L1R'.format(gatts['sensor'],  dtime.strftime('%Y_%m_%d_%H_%M_%S'), setu['viirs_option'].upper())
                gatts['viirs_option'] = setu['viirs_option']
                gatts['viirs_slines'] = slines

                if not os.path.exists(odir): os.makedirs(odir)
                ofile = '{}/{}.nc'.format(odir, gatts['obase'])

                ## get image subset
                sub = None
                if limit is not None:
                    group = 'geolocation_data'
                    with h5py.File(geo, mode='r') as f:
                        lat = f[group]['latitude'][:]
                        lon = f[group]['longitude'][:]
                    full_shape = lat.shape

                    ## get sub
                    csub = ac.shared.geolocation_sub(lat, lon, limit)
                    ## make sure we are at even pixels
                    csub = [c - c%2 for c in csub]

                    ## make sure we crop at scan lines
                    csub[1] = csub[1] - csub[1]%slines if csub[1] >= slines else 0
                    csub[3] = csub[3] - csub[3]%slines + slines
                    if (csub[1] + csub[3] > full_shape[0]): csub[3] = (full_shape[0]-csub[1])

                    ## store crop subset
                    if opts[0] == 'mod':
                        sub = {}
                        sub['mod'] = [c for c in csub]
                        sub['img'] = [int(c*2) for c in csub]
                    elif opts[0] == 'img':
                        sub = {}
                        sub['mod'] = [int(c/2) for c in csub]
                        sub['img'] = [c for c in csub]

            ## load GEO data
            if vi == 0:
                with h5py.File(geo, mode='r') as f:
                    datasets = {'lat': {'name': 'latitude', 'group': 'geolocation_data'},
                                'lon': {'name': 'longitude', 'group': 'geolocation_data'},
                                'sza': {'name': 'solar_zenith', 'group': 'geolocation_data'},
                                'saa': {'name': 'solar_azimuth', 'group': 'geolocation_data'},
                                'vza': {'name': 'sensor_zenith', 'group': 'geolocation_data'},
                                'vaa': {'name': 'sensor_azimuth', 'group': 'geolocation_data'}
                               }

                    data = {}
                    for ds in datasets:
                        dds = f[datasets[ds]['group']][datasets[ds]['name']]
                        atts = {k: dds.attrs[k] for k in dds.attrs.keys()}
                        print(ds)

                        if sub is None:
                            data[ds] = dds[:]
                        else:
                            #data[ds] = dds[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
                            data[ds] = dds[sub[viirs_res][1]:sub[viirs_res][1]+sub[viirs_res][3],\
                                           sub[viirs_res][0]:sub[viirs_res][0]+sub[viirs_res][2]]

                        if data[ds].dtype in [np.int16, np.int32]:
                            data[ds] = (data[ds].astype(np.float32) *  atts['scale_factor']) + atts['add_offset']

                    ## compute relative azimuth
                    data['raa'] = np.abs(data['saa'] - data['vaa'])
                    data['raa'][data['raa']>180] = 360 - data['raa'][data['raa']>180]

                ## store means
                for ds in data:
                    gatts[ds] = np.nanmean(data[ds])

                ## store cos sun zenith angle
                mus = np.cos(np.radians(data['sza']))
                if rotate: mus = np.rot90(mus, k=2)

                ## track current shape
                data_shape = mus.shape

                for ds in datasets:
                    if ds not in data: continue
                    if (setu['output_geolocation']) & (ds in ['lat', 'lon']):
                        if verbosity > 1: print('Writing geolocation {}'.format(ds))
                    elif (setu['output_geometry']) & (ds in ['sza', 'saa', 'vza', 'vaa', 'raa']):
                        if verbosity > 1: print('Writing geometry {}'.format(ds))
                    else: continue

                    if rotate:
                        #data[ds] = np.flip(np.rot90(data[ds]))
                        data[ds] = np.rot90(data[ds], k=2)

                    ac.output.nc_write(ofile, ds, data[ds], new=new, attributes=gatts,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'])
                    if verbosity > 1: print('Wrote {} ({})'.format(ds, data[ds].shape))
                    new = False
                    del data[ds]
                data = None
            ## end store geometry

            ## read TOA data
            with h5py.File(l1b, mode='r') as f:
                group = 'observation_data'
                use_radiance = False

                for b in rsrd[sensor]['rsr_bands']:
                    ds = 'rhot_{}'.format(rsrd[sensor]['wave_name'][b])
                    if setu['add_band_name']: ds = 'rhot_{}_{}'.format(b, rsrd[sensor]['wave_name'][b])
                    #print(b, ds)

                    if b in f[group].keys():
                        dds = f[group][b]
                        atts = {k: dds.attrs[k] for k in dds.attrs.keys()}

                        if sub is None:
                            data = dds[:]
                        else:
                            #data = dds[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
                            data = dds[sub[viirs_res][1]:sub[viirs_res][1]+sub[viirs_res][3],\
                                           sub[viirs_res][0]:sub[viirs_res][0]+sub[viirs_res][2]]


                        ## get mask
                        mask = data >= 65533
                        data = data.astype(np.float32)
                        data[mask] = np.nan

                        ## convert to radiance or reflectance
                        if (use_radiance) & ('radiance_scale_factor' in atts):
                            data = (data * atts['radiance_scale_factor']) + atts['radiance_add_offset']
                            stop
                        elif (not use_radiance) & ('scale_factor' in atts):
                            data = (data * atts['scale_factor']) + atts['add_offset']
                            ## zoom data
                            if data.shape != data_shape:
                                data = scipy.ndimage.zoom(data, data_shape[0]/data.shape[0], order=1) ## order = 1 is linear interp
                            data /= mus
                        else:
                            print('Could not convert to radiance/reflectance')
                            continue

                        ds_att = {}
                        ds_att['band'] = b
                        ds_att['rhot_ds'] = ds
                        ds_att['wavelength'] = rsrd[sensor]['wave_nm'][b]

                        if rotate:
                            #data = np.flip(np.rot90(data))
                            data = np.rot90(data, k=2)

                        ## output dataset
                        ac.output.nc_write(ofile, ds, data, new=new, attributes=gatts,
                                            dataset_attributes = ds_att,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'])
                        if verbosity > 1: print('Wrote {} ({})'.format(ds, data.shape))
                        new = False
                        data = None

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
