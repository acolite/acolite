## def l1_convert
## converts VIIRS data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications: 2023-04-18 (QV) check if outside limits
##                2023-04-19 (QV) added quality_flags check
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling, fixed raa writing, fix for new settings handling
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-04 (QV) added acolite_file_type, improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming
##                2025-02-13 (QV) changed to settings.merge, added tile merging

def l1_convert(inputfile, output = None, settings = None):
    import h5py
    import numpy as np
    import scipy.ndimage
    import dateutil.parser, time
    import glob, os
    import acolite as ac

    ## datasets in GEO file
    geo_datasets = {'lat': {'name': 'latitude', 'group': 'geolocation_data'},
                    'lon': {'name': 'longitude', 'group': 'geolocation_data'},
                    'sza': {'name': 'solar_zenith', 'group': 'geolocation_data'},
                    'saa': {'name': 'solar_azimuth', 'group': 'geolocation_data'},
                    'vza': {'name': 'sensor_zenith', 'group': 'geolocation_data'},
                    'vaa': {'name': 'sensor_azimuth', 'group': 'geolocation_data'}
                    }

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

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
    if setu['verbosity'] > 1: print('Starting conversion of {} scene{}'.format(nscenes, '' if nscenes==1 else 's'))

    ## list to store output files
    ofiles = []

    if setu['merge_tiles'] & (nscenes == 1):
        if setu['verbosity'] > 1: print('One scene provided, and merge_tiles=True. Setting merge_tiles=False.')
        setu['merge_tiles'] = False
        ac.settings['run']['merge_tiles'] = False

    ## test if we need to merge
    if setu['merge_tiles']:
        if setu['verbosity'] > 1: print('Testing whether {} scene{} can be merged'.format(nscenes, '' if nscenes==1 else 's'))
        ## get crop in mod data - faster
        ret = ac.viirs.viirs_merge_test(inputfile, limit = setu['limit'], opt = 'mod')
        if ret is None: ## return with no result
            return(ofiles, setu)
        else: ## unpack returns
            ## get mod merged data subsets
            sub_mod, data_shape_merged_mod, sort_bundles, crop_in_mod, crop_out_mod= ret
            ## sort inputfiles
            inputfile = [inputfile[bi] for bi in sort_bundles]
            ## get img data merged subsets
            sub_img = [s*2 for s in sub_mod]
            crop_in_img = [[s*2 for s in c] for c in crop_in_mod]
            crop_out_img = [[s*2 for s in c] for c in crop_out_mod]
            data_shape_merged_img = [s*2 for s in data_shape_merged_mod]
    ## end test merging

    ## run through bundles
    new = True
    outside = False
    for bi, bundle in enumerate(inputfile):
        ## make new file and do extent test if not merging
        if not setu['merge_tiles']:
            new = True
            outside = False

        files = ac.viirs.bundle_test(bundle)
        l1b_meta = None
        for mkey in ['l1b_mod', 'l1b_img']:
            if mkey in files:
                l1b_meta = ac.shared.nc_gatts(files[mkey])
                break
        l1b_sensor = '{}_{}'.format(l1b_meta['platform'], l1b_meta['instrument']).upper()
        print(l1b_sensor, l1b_meta['title'], l1b_meta['processing_level'])

        ## update settings
        sensor = '{}_{}'.format(l1b_meta['platform'], l1b_meta['instrument']).upper()
        setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

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

        verbosity = setu['verbosity']
        if output is None: output = setu['output']

        ## load rsr
        rsrd = ac.shared.rsr_dict(sensor=sensor)
        bands_vswir = [b for b in rsrd[sensor]['rsr_bands']]
        output_bands = [b for b in bands_vswir]
        if setu['viirs_output_tir']:
            sensor_tir = '{}_TIR'.format(sensor)
            rsrd_tir = ac.shared.rsr_dict(sensor=sensor_tir, wave_range = [3, 15])
            bands_tir = [b for b in rsrd_tir[sensor_tir]['rsr_bands']]
            output_bands += [b for b in bands_tir]

        ## track gain index for RSB (0-2 are I bands 1-3, and 3-13 are M bands)
        gi = 0

        ## run through viirs options
        for vi, viirs_res in enumerate(opts):
            if outside: continue

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
                ## empty dict for if new file tile
                data = {}
                data_att = {}

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
                gatts['viirs_option'] = setu['viirs_option']
                gatts['viirs_slines'] = slines
                gatts['acolite_file_type'] = 'L1R'

                ## output name
                oname  = '{}_{}_{}'.format(gatts['sensor'],  dtime.strftime('%Y_%m_%d_%H_%M_%S'), setu['viirs_option'].upper())
                if setu['merge_tiles']: oname+='_merged'
                if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
                ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
                gatts['oname'] = oname
                gatts['ofile'] = ofile

                ## new output gem
                gemo = ac.gem.gem(ofile, new = True)
                gemo.gatts = {k: gatts[k] for k in gatts}
                new = False

                ## get image subset if not merging
                sub = None
                if (setu['limit'] is not None) & (setu['merge_tiles'] is False):
                    group = 'geolocation_data'
                    with h5py.File(geo, mode='r') as f:
                        lat = f[group]['latitude'][:]
                        lon = f[group]['longitude'][:]
                    full_shape = lat.shape

                    ## get sub
                    csub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
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

            #### start of tile merging
            if setu['merge_tiles']:
                ## get correct shape and crops for output
                if opts[0] == 'mod':
                    co = [c for c in crop_out_mod[bi]]
                    cur_shape = crop_in_mod[bi][3]-crop_in_mod[bi][2], crop_in_mod[bi][1]-crop_in_mod[bi][0]
                    out_shape = sub_mod[3], sub_mod[2]
                elif opts[0] == 'img':
                    co = [c for c in crop_out_img[bi]]
                    cur_shape = crop_in_img[bi][3]-crop_in_img[bi][2], crop_in_img[bi][1]-crop_in_img[bi][0]
                    out_shape = sub_img[3], sub_img[2]
                ## get image subset crop for input
                if viirs_res == 'mod':
                    tsub = crop_in_mod[bi][0], crop_in_mod[bi][2], crop_in_mod[bi][1]-crop_in_mod[bi][0], crop_in_mod[bi][3]-crop_in_mod[bi][2]
                elif viirs_res == 'img':
                    tsub = crop_in_img[bi][0], crop_in_img[bi][2], crop_in_img[bi][1]-crop_in_img[bi][0], crop_in_img[bi][3]-crop_in_img[bi][2]
            else: ## non merged version here
                ## output dimensions in selected option
                co = [0, sub[opts[0]][2], 0, sub[opts[0]][3]]
                out_shape = sub[opts[0]][3], sub[opts[0]][2]
                cur_shape = sub[opts[0]][3], sub[opts[0]][2]
                ## sub in current image
                tsub = sub[viirs_res]

            ## read GEO data
            if vi == 0:
                with h5py.File(geo, mode='r') as f:
                    for ds in geo_datasets:
                        dds = f[geo_datasets[ds]['group']][geo_datasets[ds]['name']]
                        atts = {k: dds.attrs[k] for k in dds.attrs.keys()}

                        ## add array to data dict
                        if ds not in data: data[ds] = np.zeros(out_shape) + np.nan

                        ## read data
                        data_read = dds[tsub[1]:tsub[1]+tsub[3],tsub[0]:tsub[0]+tsub[2]]
                        if data_read.dtype in [np.int16, np.int32]:
                            data_read = (data_read.astype(np.float32) *  atts['scale_factor']) + atts['add_offset']
                        data[ds][co[2]:co[3], co[0]:co[1]] = data_read
                        del data_read
            ## end read GEO data

            ## read observation data
            with h5py.File(l1b, mode='r') as f:
                group = 'observation_data'

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

                        ## read data and quality
                        data_read = dds[tsub[1]:tsub[1]+tsub[3],tsub[0]:tsub[0]+tsub[2]]
                        ql = mds[tsub[1]:tsub[1]+tsub[3],tsub[0]:tsub[0]+tsub[2]]

                        ## get mask
                        mask = (data_read >= 65532) | ((ql & (viirs_quality_flags_sum)) != 0)
                        data_read = data_read.astype(np.float32)

                        ## VSWIR rhot data
                        if b in bands_vswir:
                            ds_att = {}
                            ds_att['band'] = b
                            ds_att['rhot_ds'] = ds
                            ds_att['wavelength'] = rsrd[sensor]['wave_nm'][b]

                            ## convert to reflectance
                            if ('scale_factor' in atts):
                                data_read = (data_read * atts['scale_factor']) + atts['add_offset']
                            else:
                                print('Could not convert to reflectance')
                                continue

                        ## TIR BT data
                        elif b in bands_tir:
                            ds_att = {}
                            ds_att['band'] = b
                            ds_att['bt_ds'] = ds
                            ds_att['wavelength'] = rsrd_tir[sensor_tir]['wave_nm'][b]

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
                            data_read = np.interp(data_read, xdlut, btlut)

                        ## mask data
                        data_read[mask] = np.nan

                        ## zoom data to output opt
                        if data_read.shape != cur_shape:
                            data_read = scipy.ndimage.zoom(data_read, cur_shape[0]/data_read.shape[0], order=1) ## order = 1 is linear interp

                        ## apply gains
                        if ('rhot_' in ds) & (setu['gains']):
                            ds_att['gain_toa'] = setu['gains_toa'][gi]
                            ds_att['gain_applied'] = 1
                            if setu['verbosity'] > 0: print('Applying gain {} to {}'.format(setu['gains_toa'][gi], ds))
                            data_read *= setu['gains_toa'][gi]
                            gi+=1

                        ## add array to data dict
                        if ds not in data: data[ds] = np.zeros(out_shape) + np.nan
                        data[ds][co[2]:co[3], co[0]:co[1]] = data_read
                        if ds not in data_att: data_att[ds] = ds_att

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
                            data_read = ds_att['K1_CONSTANT_BAND_{}'.format(b)]/(np.exp(ds_att['K2_CONSTANT_BAND_{}'.format(b)]/data_read)-1)
                            ## add array to data dict
                            if ds not in data: data[ds] = np.zeros(cur_shape) + np.nan
                            data[ds][co[2]:co[3], co[0]:co[1]] = data_read
                            if ds not in data_att: data_att[ds] = ds_att
                        del data_read
                        del ds_att
            ## end read observation data

            ## check if we need to write the data to NetCDF
            store = (vi == len(opts)-1)
            if setu['merge_tiles']:
                store = store & (bi == len(inputfile)-1)

            ## write data if last bundle and option
            if store:
                ## compute relative azimuth
                data['raa'] = np.abs(data['saa'] - data['vaa'])
                data['raa'][data['raa']>180] = 360 - data['raa'][data['raa']>180]

                ## store means
                for ds in geo_datasets: gatts[ds] = np.nanmean(data[ds])

                ## store cos sun zenith angle
                mus = np.cos(np.radians(data['sza']))

                ## track current shape
                data_shape = mus.shape

                datasets_ = list(data.keys())
                for ds in datasets_:
                    if ds not in data: continue
                    if (setu['output_geolocation']) & (ds in ['lat', 'lon']):
                        if setu['verbosity'] > 1: print('Writing geolocation {}'.format(ds))
                    elif (setu['output_geometry']) & (ds in ['sza', 'saa', 'vza', 'vaa', 'raa']):
                        if setu['verbosity'] > 1: print('Writing geometry {}'.format(ds))
                    else:
                        if setu['verbosity'] > 1: print('Writing dataset {}'.format(ds))
                    ## scale rhot with mus
                    if ds.startswith('rhot_'): data[ds] /= mus
                    ## rotate dataset
                    if rotate: data[ds] = np.rot90(data[ds], k=2)
                    ## write dataset
                    gemo.write(ds, data[ds])
                    if setu['verbosity'] > 1: print('Wrote {} ({})'.format(ds, data[ds].shape))
                    del data[ds]
                data = None

                if setu['limit'] is not None: sub = None
                if not outside:
                    gemo.close()
                    if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
