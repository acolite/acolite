## def zarr.l1_convert
## converts Sentinel-3 zarr dataset to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2025-11-04
## modifications: 2025-11-05 (QV) added scene merging

def l1_convert(inputfile, output = None, settings = None, write_l2_err = False):

    import os, glob, datetime, time, re
    import dateutil.parser

    import numpy as np
    import scipy.interpolate
    import acolite as ac
    import zarr

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

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
        ret = ac.sentinel3.zarr.olci_merge_test(inputfile, limit = setu['limit'], use_tpg = setu['use_tpg'])
        if ret is None: ## return with no result
            return(ofiles, setu)
        else: ## unpack returns
            sub_merged, data_shape_merged, sort_bundles, crop_in, crop_out = ret
            inputfile = [inputfile[bi] for bi in sort_bundles]

    ## start conversion
    new = True
    for bi, bundle in enumerate(inputfile):
        t0 = time.time() ## track time
        print(bi, bundle)

        ## open zarr bundle
        z = zarr.open(bundle, mode='r')

        ## get zarr metadata
        meta = ac.zarr.meta_parse(z)

        ## sensor settings
        sensor = meta['sensor']
        if sensor in ['S3A_OLCI', 'S3B_OLCI']:
            len_gains = 21
            bands_data = ac.sentinel3.olci_band_info()
            band_id = 'Oa'
        elif sensor in ['EN1_MERIS']:
            len_gains = 15
            bands_data = ac.sentinel3.meris_band_info()
            band_id = 'M'
        else:
            print('Sensor {} from file {} not configured.'.format(sensor, bundle))
            continue

        ## track if RR or FR
        s3_product_type = meta['product:type'][-2:]
        if 'proj:shape' in meta:
            s3_product_type = meta['proj:shape']['name']
        print(sensor, s3_product_type)

        ## load rsrd
        rsrd = ac.shared.rsr_dict(sensor)
        rsr_bands = rsrd[sensor]['rsr_bands']

        ## update settings
        setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

        if output is None: output = setu['output']

        ## for consistency with NetCDF header subsampling factors
        ac_subsampling_factor = meta['data_information']['sampling']['columns_per_tiepoint']
        al_subsampling_factor = meta['data_information']['sampling']['rows_per_tiepoint']

        ## get tpg and data shape
        tpg_shape = z['conditions']['geometry']['latitude'].shape
        data_shape = z['conditions']['image']['latitude'].shape

        ## if not merging tiles, create a new file
        if not setu['merge_tiles']:
            new = True
            sub = None
            ## find subset in scene
            if setu['limit'] is not None:
                if setu['use_tpg']:
                    print('Reading lat/lon TPG to determine subset')
                    lat_tpg = z['conditions']['geometry']['latitude'][:]
                    lon_tpg = z['conditions']['geometry']['longitude'][:]
                    sub_tpg = ac.shared.geolocation_sub(lat_tpg, lon_tpg, setu['limit'])
                    if sub_tpg is None:
                        sub = None
                    else:
                        sub = [s for s in sub_tpg]
                        ys = int((data_shape[1]-1)/(tpg_shape[1]-1))
                        y0s=max(0,sub[0]-1)*ys
                        y1s=min(tpg_shape[1],sub[0]+sub[2]+1)*ys
                        sub[0] = int(y0s)
                        sub[2] = int(y1s-y0s)
                else:
                    print('Reading lat/lon to determine subset')
                    lat = z['conditions']['image']['latitude'][:]
                    lon = z['conditions']['image']['longitude'][:]
                    sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])

                ## sub is None if out of scene
                if sub is None:
                    if setu['verbosity'] > 1: print('Limit {} out of scene {}'.format(setu['limit'], bundle))
                    continue

                ## limit to scene
                if sub is not None:
                    if sub[0]<0: sub[0] = 0
                    if sub[1]<0: sub[1] = 0
                    if sub[1]+sub[3] > data_shape[0]:
                        sub[3] = data_shape[0]-sub[1]
                    if sub[0]+sub[2] > data_shape[1]:
                        sub[2] = data_shape[1]-sub[0]
            ## user defined sub if limit or polygon are not set
            if (sub is None) & (setu['sub'] is not None): sub = setu['sub']
        else:
            ## sub in current bundle
            if setu['limit'] is None:
                sub = None
                data_shape = data_shape_merged[0],  data_shape_merged[1]
            else:
                sub = crop_in[bi][0], crop_in[bi][2], crop_in[bi][1]-crop_in[bi][0], crop_in[bi][3]-crop_in[bi][2]
                data_shape = sub_merged[3], sub_merged[2]

        ## new empty datasets
        if new:
            lfiles = {}
            data, metadata, tpg_data = {}, {}, {}
            instrument_data, image_data = {}, {}

        ## get instrument data
        instrument_datasets = [m[0] for m in z['conditions']['instrument'].members()]
        for k in instrument_datasets:
            if not setu['merge_tiles']:
                instrument_data[k] = z['conditions']['instrument'][k][:]
            else:
                ## stack the data if merging tiles
                if k not in instrument_data:
                    instrument_data[k] = z['conditions']['instrument'][k][:]
                else:
                    if len(instrument_data[k].shape) == 1:
                        ## use hstack if one dimension, in case the number of along track pixels differs
                        instrument_data[k] = np.hstack((instrument_data[k], z['conditions']['instrument'][k][:]))
                    else:
                        instrument_data[k] = np.vstack((instrument_data[k], z['conditions']['instrument'][k][:]))

        ## get image data
        #image_datasets = [m[0] for m in z['conditions']['image'].members()]
        image_datasets = ['detector_index']
        for k in image_datasets:
            if not setu['merge_tiles']:
                image_data[k] = z['conditions']['image'][k][:]
            else:
                ## stack the data if merging tiles
                if k not in image_data:
                    image_data[k] = z['conditions']['image'][k][:]
                else:
                    if len(image_data[k].shape) == 1:
                        ## use hstack if one dimension, in case the number of along track pixels differs
                        image_data[k] = np.hstack((image_data[k], z['conditions']['image'][k][:]))
                    else:
                        image_data[k] = np.vstack((image_data[k], z['conditions']['image'][k][:]))

        ## get tpg data
        geometry_datasets = [m[0] for m in z['conditions']['geometry'].members()]
        for k in geometry_datasets:
            if not setu['merge_tiles']:
                tpg_data[k] = z['conditions']['geometry'][k][:]
            else:
                ## stack the data if merging tiles
                if k not in tpg_data:
                    tpg_data[k] = z['conditions']['geometry'][k][:]
                else:
                    if len(tpg_data[k].shape) == 1:
                        ## use hstack if one dimension, in case the number of along track pixels differs
                        tpg_data[k] = np.hstack((tpg_data[k], z['conditions']['geometry'][k][:]))
                    else:
                        tpg_data[k] = np.vstack((tpg_data[k], z['conditions']['geometry'][k][:]))
                print(k, tpg_data[k].shape)

        meteorology_datasets = [m[0] for m in z['conditions']['meteorology'].members()]
        for k in meteorology_datasets:
            if k in ['latitude', 'longitude']: continue ## skip datasets loaded from geometry
            if not setu['merge_tiles']:
                tpg_data[k] = z['conditions']['meteorology'][k][:]
            else:
                ## stack the data if merging tiles
                if k not in tpg_data:
                    tpg_data[k] = z['conditions']['meteorology'][k][:]
                else:
                    if len(tpg_data[k].shape) == 1:
                        ## use hstack if one dimension, in case the number of along track pixels differs
                        tpg_data[k] = np.hstack((tpg_data[k], z['conditions']['meteorology'][k][:]))
                    else:
                        tpg_data[k] = np.vstack((tpg_data[k], z['conditions']['meteorology'][k][:]))

        ## read data
        datasets = [m[0] for m in z['measurements'].members()]
        for k in datasets:
            if k in ['orphans', 'time_stamp']: continue
            print('Reading {}'.format(k))
            if sub is None:
                d = z['measurements'][k][:]
            else:
                d = z['measurements'][k][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
            md = z['measurements'][k].metadata.to_dict()
            if k.endswith('_radiance'): k = 'O{}'.format(k[1:])

            if not setu['merge_tiles']:
                data[k] = d
            else:
                if k not in data: data[k] = np.zeros(data_shape) + np.nan
                data[k][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = d
            #if setu['verbosity'] > 2: print(k, data[k].shape)

            #data[k] = d
            #metadata[k] = md
            del d, md
        ## end read data

        ## if not final tile set new to False and continue
        if (setu['merge_tiles']) & (bi < len(inputfile)-1):
            new = False
            continue

        ## determine product level
        if meta['processing:level'] == 'L1':
            product_level = 'level1'
            acolite_file_type = 'L1R'

        ## global attributes
        dtime = dateutil.parser.parse(meta['start_datetime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## interpolate TPG
        if (setu['merge_tiles']):
            if (setu['limit'] is not None):
                subx = np.arange(sub_merged[0],sub_merged[0]+sub_merged[2]) + 0.5
                suby = np.arange(sub_merged[1],sub_merged[1]+sub_merged[3]) + 0.5
            else:
                subx = np.arange(0,data_shape_merged[1]) + 0.5
                suby = np.arange(0,data_shape_merged[0]) + 0.5
        elif (sub is None):
            subx = np.arange(0,data_shape[1]) + 0.5
            suby = np.arange(0,data_shape[0]) + 0.5
        else:
            subx = np.arange(sub[0],sub[0]+sub[2]) + 0.5
            suby = np.arange(sub[1],sub[1]+sub[3]) + 0.5

        ## interpolate tie point grids
        x_sub = ac_subsampling_factor * 1.0
        y_sub = al_subsampling_factor * 1.0
        x_off = 0
        y_off = 0
        tpg_shape = tpg_data['latitude'].shape
        tpx = (np.arange(tpg_shape[1])*x_sub) + x_off
        tpy = (np.arange(tpg_shape[0])*y_sub) + y_off

        ## 2d for RGI
        subx_ = np.tile(subx, (len(suby),1))
        suby_ = (np.tile(suby, len(subx)).reshape(subx_.shape[1], subx_.shape[0])).T

        tpg = {}
        for k in tpg_data.keys():
            if k in ['atmospheric_temperature_profile', 'horizontal_wind', 'reference_pressure_level',
                     'pressure_level', 'wind_vector']: continue
            if tpg_data[k].shape != tpg_shape:
                print('{} tpg shape {} not supported'.format(k,tpg_data[k].shape))
                continue
            print('Interpolating TPG {}'.format(k))
            dtype_in = tpg_data[k].dtype
            rgi = scipy.interpolate.RegularGridInterpolator([tpy, tpx], tpg_data[k].astype(np.float64), bounds_error = False, fill_value = None)
            tpg[k] = rgi((suby_,subx_)).astype(dtype_in)
        del tpg_data

        ## compute relative azimuth TPG
        tpg['raa'] = abs(tpg['saa']-tpg['oaa'])
        tpg['raa'][tpg['raa']>180]=np.abs(360-tpg['raa'][tpg['raa']>180])

        ## cosine of sun zenith angle
        mu = np.cos(tpg['sza']*(np.pi/180))

        ## average geometry
        sza = np.nanmean(tpg['sza'])
        vza = np.nanmean(tpg['oza'])
        raa = np.nanmean(tpg['raa'])
        ## center lat and lon
        clat = np.nanmean(tpg['latitude'])
        clon = np.nanmean(tpg['longitude'])

        ## compute gas transmittance
        uoz = None
        uwv = None
        pressure = None
        if (setu['use_supplied_ancillary']) & (not setu['ancillary_data']):
            ## convert ozone from kg.m-2 to cm.atm
            setu['uoz_default'] = np.nanmean(tpg['total_ozone'])/0.02141419
            ## convert water from kg.m-2 to g.cm-2
            setu['uwv_default'] = np.nanmean(tpg['total_columnar_water_vapour'])/10
            setu['pressure'] = np.nanmean(tpg['sea_level_pressure'])
        else:
            if setu['ancillary_data']:
                print('Getting ancillary data for {} {:.3f}E {:.3f}N'.format(isodate, clon, clat))
                anc = ac.ac.ancillary.get(dtime, clon, clat, verbosity=setu['verbosity'])
                ## overwrite the defaults
                if ('uoz' in anc): setu['uoz_default'] = anc['uoz']
                if ('uwv' in anc): setu['uwv_default'] = anc['uwv']
                if ('wind' in anc) & (setu['wind'] is None): setu['wind'] = anc['wind']
                if ('pressure' in anc) & (setu['pressure'] == setu['pressure_default']): setu['pressure'] = anc['pressure']

        if uoz is None: uoz = setu['uoz_default']
        if uwv is None: uwv = setu['uwv_default']
        if pressure is None: pressure = setu['pressure']
        print('current uoz: {:.2f} uwv: {:.2f} pressure: {:.2f}'.format(uoz, uwv, pressure))

        ## for smile correction
        ttg = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv, sensor=sensor)

        ## get per pixel detector index
        if (setu['merge_tiles']):
            if (setu['limit'] is not None):
                di = image_data['detector_index'][sub_merged[1]:sub_merged[1]+sub_merged[3], sub_merged[0]:sub_merged[0]+sub_merged[2]]
            else:
                di = image_data['detector_index']
        elif (sub is None):
            di = image_data['detector_index']
        else:
            di = image_data['detector_index'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        ## convert from float to int
        di = di.astype(int)
        ## track bad detector indices
        di_fill = np.where(di<0)
        di[di_fill] = 0

        ## smile correction - from l2gen smile.c
        if (setu['smile_correction']) & (product_level == 'level1'):
            if setu['verbosity'] > 1: print('Running smile correction')
            ## gas_correction
            if setu['smile_correction_tgas']:
                if setu['verbosity'] > 2: print('{} - Gas correction before smile correction'.format(datetime.datetime.now().isoformat()[0:19]), end='\n')
                for band in bands_data: data['{}_radiance'.format(band)]/=ttg['tt_gas'][band]

            smile = {}
            for band in bands_data:
                ## dataset name
                dname = '{}_radiance'.format(band)

                ## band index
                b_i = bands_data[band]['band']-1
                if setu['verbosity'] > 2: print('{} - Smile correction for band {} {} nm'.format(datetime.datetime.now().isoformat()[0:19], band, bands_data[band]['wavelength'] ), end='\n')

                ## bounding bands
                b1_i = bands_data[band]['lower_water']-1
                b2_i = bands_data[band]['upper_water']-1
                band1 = '{}{}'.format(band_id, str(bands_data[band]['lower_water']).zfill(2))
                band2 = '{}{}'.format(band_id, str(bands_data[band]['upper_water']).zfill(2))

                ## compute reflectance using per detector F0
                r_ = (data['{}_radiance'.format(band)]) / instrument_data['solar_flux'][b_i][di]

                ## based on that reflectance, compute radiance for target F0
                r_ *= bands_data[band]['E0']

                ## difference in radiance
                smile[band] = r_-data['{}_radiance'.format(band)]
                del r_ ## free memory

                ## do additional correction based on two bounding bands
                ## currently applying water everywhere
                if bands_data[band]['switch_water'] > 0:
                    if setu['verbosity'] > 2: print('{} - Smile correction - bounding bands {}/{}'.format(datetime.datetime.now().isoformat()[0:19], band1, band2), end='\n')

                    ## compute per pixel reflectance difference for bounding bands
                    r21_diff = (data['{}_radiance'.format(band2)]) / instrument_data['solar_flux'][b2_i][di]-\
                               (data['{}_radiance'.format(band1)]) / instrument_data['solar_flux'][b1_i][di]

                    ## wavelength difference ratio
                    wdiff_ratio = (bands_data[band]['wavelength'] - instrument_data['lambda0'][b_i][di])/\
                                  (instrument_data['lambda0'][b2_i][di] - instrument_data['lambda0'][b1_i][di])

                    ## additional smile
                    smile[band] += (r21_diff)*(wdiff_ratio)*(instrument_data['solar_flux'][b_i][di])
                    del r21_diff, wdiff_ratio

            ## add smile effect to radiance data
            for band in smile: data['{}_radiance'.format(band)]+=smile[band]
            del smile
            ## add back in gas transmittance
            if setu['smile_correction_tgas']:
                if setu['verbosity'] > 2: print('{} - Gas correction restored after smile correction'.format(datetime.datetime.now().isoformat()[0:19]), end='\n')
                for band in bands_data: data['{}_radiance'.format(band)]*=ttg['tt_gas'][band]

            ## mask bad detector indices
            for band in bands_data: data['{}_radiance'.format(band)][di_fill] = np.nan
        ## end smile correction

        ## read rsr
        waves_mu = rsrd[sensor]['wave_mu']
        waves_names = rsrd[sensor]['wave_name']

        ## take wavelengths and band names from external table
        ## the smile correction should bring things in line with these
        dnames = ['{}_radiance'.format(b) for b in bands_data]
        bnames = [b for b in bands_data]

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsrd[sensor]['rsr'])

        gatts = {'sensor':sensor, 'sza':sza, 'vza':vza, 'raa':raa,
                 'isodate':isodate, 'global_dims':data_shape,
                 'se_distance': se_distance, 'acolite_file_type': 'L1R',
                 's3_product_type': s3_product_type}

        gatts['pressure'] = pressure
        gatts['uoz'] = uoz
        gatts['uwv'] = uwv

        if setu['limit'] is not None: gatts['limit'] = setu['limit']
        if sub is not None: gatts['sub'] = sub

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        oname += '_{}'.format(gatts['s3_product_type'])
        if setu['merge_tiles']: oname+='_merged'
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])

        ofile = '{}/{}_{}.nc'.format(output, oname, acolite_file_type)
        if not os.path.exists(os.path.dirname(ofile)): os.makedirs(os.path.dirname(ofile))
        gatts['oname'] = oname
        gatts['ofile'] = ofile
        gatts['acolite_file_type'] = acolite_file_type

        ## add band info to gatts
        for bbi, b in enumerate(rsr_bands):
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## set up output file
        gemo = ac.gem.gem(ofile, new=True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## output lat/lon
        if setu['output_geolocation']:
            if setu['verbosity'] > 1: print('Writing geolocation')
            for ds in ['longitude', 'latitude']: gemo.write(ds[0:3], data[ds])

        ## output geometry
        if setu['output_geometry']:
            if setu['verbosity'] > 1: print('Writing geometry')
            for k in tpg:
                if k in ['sza', 'oza', 'raa', 'saa', 'oaa']:
                    if k == 'oza':
                        ko = 'vza'
                    elif k == 'oaa':
                        ko = 'vaa'
                    else:
                        ko = k.lower()
                    gemo.write(ko, tpg[k])
                elif k in ['sea_level_pressure']:
                    gemo.write('pressure', tpg[k])
                else:
                    continue

        ## is the F0 already corrected for Sun - Earth distance?
        ## see OLCI_L2_ATBD_Instrumental_Correction.pdf
        #se = ac.shared.distance_se(doy)
        se = 1.
        se2 = se**2

        ## read TOA
        if (product_level == 'level1'):
            if setu['verbosity'] > 1: print('Writing TOA reflectance')
            for iw, band in enumerate(rsr_bands):
                wave = waves_names[band]
                ds = 'rhot_{}'.format(wave)
                if setu['verbosity'] > 2: print('{} - Converting TOA data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')

                # per pixel wavelength
                l = instrument_data['lambda0'][iw][di]
                # per pixel f0
                f0 = instrument_data['solar_flux'][iw][di]
                # if smile corrected use nominal E0
                if setu['smile_correction']: f0 = bands_data[band]['E0']

                # per pixel fwhm
                fwhm = instrument_data['fwhm'][iw][di]
                dname = dnames[iw]

                ds_att  = {'wavelength':float(wave)}
                for key in ttg: ds_att[key]=ttg[key][bnames[iw]]

                ## write toa radiance
                if setu['output_lt']:
                    gemo.write('Lt_{}'.format(wave), data[dname], ds_att = ds_att)
                    if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format('Lt_{}'.format(wave), data[dname].shape))

                ## convert to reflectance
                print(mu.shape)
                d = (np.pi * data[dname] * se2) / (f0*mu)
                ## write dataset
                gemo.write(ds, d, ds_att = ds_att)
                if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
                d = None

        ## read BOA
        if (product_level == 'level2'):
            print('level2 conversion for ZARR not yet implemented')

        ## clear data
        del data, tpg
        del instrument_data, image_data

        ## close file
        gemo.close()
        gemo = None

        if setu['verbosity'] > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
