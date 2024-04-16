## def l1_convert
## converts VIIRS data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications: 2023-04-18 (QV) check if outside limits
##                2023-04-19 (QV) added quality_flags check
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling, fixed raa writing, fix for new settings handling

def l1_convert(inputfile, output = None, settings = {}, verbosity = 0):
    import h5py
    import numpy as np
    import scipy.ndimage
    import dateutil.parser, time
    import glob, os
    import acolite as ac

    ## function to optimise BT lut
    from scipy import optimize
    def fun(l):
        k1, k2 = l[0], l[1]
        tmp = (k2/np.log((k1/ltlut[lutsub])+1))
        v = np.sqrt(np.nanmean(np.square((tmp-btlut[lutsub]))))
        return(v)

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
        setu = ac.acolite.settings.parse(sensor, settings = ac.settings['user'])
        opts = setu['viirs_option'].lower().split('+')
        ## configure scan lines
        if opts[0] == 'mod': slines = 16
        if opts[0] == 'img': slines = 32

        ## get quality flags
        viirs_quality_flags_sum = np.sum(setu['viirs_quality_flags'])

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
        if (limit is not None) & (setu['limit_buffer'] is not None):
            print('Applying limit buffer {}'.format(setu['limit_buffer']))
            print('Old limit: {}'.format(limit))
            setu['limit_old'] = limit
            limit = limit[0] - setu['limit_buffer'], limit[1] - setu['limit_buffer'], \
                    limit[2] + setu['limit_buffer'], limit[3] + setu['limit_buffer']
            print('New limit: {}'.format(limit))

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        output_lt = setu['output_lt']
        vname = setu['region_name']
        odir = '{}'.format(output)

        ## load rsr
        rsrd = ac.shared.rsr_dict(sensor=sensor)
        bands_vswir = [b for b in rsrd[sensor]['rsr_bands']]
        output_bands = [b for b in bands_vswir]
        if setu['viirs_output_tir']:
            sensor_tir = '{}_TIR'.format(sensor)
            rsrd_tir = ac.shared.rsr_dict(sensor=sensor_tir, wave_range = [3, 15])
            bands_tir = [b for b in rsrd_tir[sensor_tir]['rsr_bands']]
            output_bands += [b for b in bands_tir]

        new = True
        outside = False

        ## track gain index for RSB (0-2 are I bands 1-3, and 3-13 are M bands)
        gi = 0

        ## run through viirs options
        for vi, viirs_res in enumerate(opts):
            if outside: continue
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
                gatts['ofile'] = ofile

                ## new output gem
                gemo = ac.gem.gem(ofile, new = True)
                gemo.gatts = {k: gatts[k] for k in gatts}
                new = False

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
                    if csub is None:
                        print('Limit outside of scene {}'.format(bundle))
                        outside = True
                        continue
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

                    ## add subset to gatts
                    gatts['viirs_mod_sub'] = sub['mod']
                    gatts['viirs_img_sub'] = sub['img']

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
                #if rotate: mus = np.rot90(mus, k=2)

                ## track current shape
                data_shape = mus.shape
                datasets_ = list(data.keys())
                for ds in datasets_:
                    if ds not in data: continue
                    if (setu['output_geolocation']) & (ds in ['lat', 'lon']):
                        if verbosity > 1: print('Writing geolocation {}'.format(ds))
                    elif (setu['output_geometry']) & (ds in ['sza', 'saa', 'vza', 'vaa', 'raa']):
                        if verbosity > 1: print('Writing geometry {}'.format(ds))
                    else: continue
                    if rotate: data[ds] = np.rot90(data[ds], k=2)
                    gemo.write(ds, data[ds])
                    if verbosity > 1: print('Wrote {} ({})'.format(ds, data[ds].shape))
                    del data[ds]
                data = None
            ## end store geometry

            ## read TOA data
            with h5py.File(l1b, mode='r') as f:
                group = 'observation_data'
                use_radiance = False

                for b in output_bands:
                    if b in bands_vswir:
                        ds = 'rhot_{}'.format(rsrd[sensor]['wave_name'][b])
                        if setu['add_band_name']: ds = 'rhot_{}_{}'.format(b, rsrd[sensor]['wave_name'][b])
                    elif b in bands_tir:
                        ds = 'bt{}'.format(b)

                    if b in f[group].keys():
                        dds = f[group][b]
                        atts = {k: dds.attrs[k] for k in dds.attrs.keys()}

                        ## read band quality flags
                        mds = f[group]['{}_quality_flags'.format(b)]
                        matts = {k: mds.attrs[k] for k in mds.attrs.keys()}

                        if sub is None:
                            data = dds[:]
                            ql = mds[:]
                        else:
                            #data = dds[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
                            data = dds[sub[viirs_res][1]:sub[viirs_res][1]+sub[viirs_res][3],\
                                           sub[viirs_res][0]:sub[viirs_res][0]+sub[viirs_res][2]]
                            ql = mds[sub[viirs_res][1]:sub[viirs_res][1]+sub[viirs_res][3],\
                                           sub[viirs_res][0]:sub[viirs_res][0]+sub[viirs_res][2]]

                        ## get mask
                        mask = (data >= 65532) | ((ql & (viirs_quality_flags_sum)) != 0)
                        data = data.astype(np.float32)

                        ## track if we need to scale by mus
                        scale_mus = False
                        if b in bands_vswir:
                            ds_att = {}
                            ds_att['band'] = b
                            ds_att['rhot_ds'] = ds
                            ds_att['wavelength'] = rsrd[sensor]['wave_nm'][b]

                            ## convert to radiance or reflectance
                            if (use_radiance) & ('radiance_scale_factor' in atts):
                                data = (data * atts['radiance_scale_factor']) + atts['radiance_add_offset']
                                stop
                            elif (not use_radiance) & ('scale_factor' in atts):
                                data = (data * atts['scale_factor']) + atts['add_offset']
                                scale_mus = True
                            else:
                                print('Could not convert to radiance/reflectance')
                                continue
                        elif b in bands_tir:
                            ds_att = {}
                            ds_att['band'] = b
                            ds_att['bt_ds'] = ds
                            ds_att['wavelength'] = rsrd_tir[sensor_tir]['wave_nm'][b]

                            # viirs_output_tir_lt = False
                            # if viirs_output_tir_lt:
                            #     data_lt = data * 1
                            #     data_lt[mask] = np.nan
                            #
                            #     ## convert to radiance or reflectance
                            #     if ('scale_factor' in atts):
                            #         data = (data * atts['scale_factor']) + atts['add_offset']
                            #     else:
                            #         print('Could not convert to radiance/reflectance')
                            #         continue

                            ## read bt lut
                            btlut = f[group]['{}_brightness_temperature_lut'.format(b)][:]
                            ## bt lut dimension
                            xdlut = np.arange(0, 327681 if b == 'M13' else 2**16)

                            ## lut dimension in radiance
                            ltlut = (xdlut * atts['scale_factor']) + atts['add_offset']

                            ## fit btlut to get K1 and K2
                            lutsub = np.where(btlut > 0)
                            xro = optimize.minimize(fun, [1000,1000])
                            k1, k2 = xro.x[0], xro.x[1]
                            ds_att['K1_CONSTANT_BAND_{}'.format(b)] = k1
                            ds_att['K2_CONSTANT_BAND_{}'.format(b)] = k2

                            ## convert to BT
                            data = np.interp(data, xdlut, btlut)

                        ## mask data
                        data[mask] = np.nan

                        ## zoom data
                        if data.shape != data_shape:
                            data = scipy.ndimage.zoom(data, data_shape[0]/data.shape[0], order=1) ## order = 1 is linear interp

                        ## scale by mus
                        if scale_mus: data /= mus

                        ## rotate dataset
                        if rotate: data = np.rot90(data, k=2)

                        ## apply gains
                        if ('rhot_' in ds) & (setu['gains']):
                            ds_att['gain_toa'] = setu['gains_toa'][gi]
                            ds_att['gain_applied'] = 1
                            if verbosity > 0: print('Applying gain {} to {}'.format(setu['gains_toa'][gi], ds))
                            data *= setu['gains_toa'][gi]
                            gi+=1

                        ## output dataset
                        gemo.write(ds, data, ds_att = ds_att)
                        if verbosity > 1: print('Wrote {} ({})'.format(ds, data.shape))

                        ## output Lt
                        if (b in bands_tir) & (setu['viirs_output_tir_lt']):
                            ds = 'Lt{}'.format(b)
                            ds_att = {}
                            ds_att['band'] = b
                            ds_att['lt_ds'] = ds
                            ds_att['wavelength'] = rsrd_tir[sensor_tir]['wave_nm'][b]
                            ds_att['K1_CONSTANT_BAND_{}'.format(b)] = k1
                            ds_att['K2_CONSTANT_BAND_{}'.format(b)] = k2
                            ## convert from BTto Lt
                            data = ds_att['K1_CONSTANT_BAND_{}'.format(b)]/(np.exp(ds_att['K2_CONSTANT_BAND_{}'.format(b)]/data)-1)
                            gemo.write(ds, data, ds_att = ds_att)
                            if verbosity > 1: print('Wrote {} ({})'.format(ds, data.shape))

                        ## delete datasets
                        data = None
                        ql = None

        if limit is not None: sub = None
        if not outside:
            gemo.close()
            if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
