## def l1_convert
## converts Sentinel-3/OLCI or ENVISAT/MERIS data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-24
## modifications: 2021-12-22 (QV) added MERIS processing
##                2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-09-29 (QV) changed rsr loading
##                2022-10-06 (QV) added back in tgas after smile correction
##                2022-11-03 (QV) update sensor/metadata parsing
##                2022-11-09 (QV) change of F0 in radiance to reflectance correction, depending on whether smile correction was applied
##                2023-02-14 (QV) added Lt outputs, user defined subset,
##                                fixed tpg interpolation, update from scipy.interpolate.interp2d to scipy.interpolate.RegularGridInterpolator
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2023-11-25 (QV) added optional l2 conversion
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-02 (QV) removed percentiles
##                2025-02-04 (QV) improved settings handling
##                2025-02-07 (QV) added tile merging
##                2025-02-08 (QV) fixed for full tile merging
##                2025-03-19 (QV) added s3_product_type, added s3_product_type to oname

def l1_convert(inputfile, output = None, settings = None, convert_l2 = False, write_l2_err = False):

    import os, glob, datetime, time, re
    import dateutil.parser

    import numpy as np
    import scipy.interpolate
    import acolite as ac

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
        ret = ac.sentinel3.olci_merge_test(inputfile, limit = setu['limit'], use_tpg = setu['use_tpg'])
        if ret is None: ## return with no result
            return(ofiles, setu)
        else: ## unpack returns
            sub_merged, data_shape_merged, sort_bundles, crop_in, crop_out = ret
            inputfile = [inputfile[bi] for bi in sort_bundles]

    ## start conversion
    new = True
    for bi, bundle in enumerate(inputfile):
        print(bi, bundle)
        nbundle = glob.glob('{}/*.SEN3'.format(bundle))
        if len(nbundle) == 1:
            print('Found nested SEN3 bundle {}'.format(nbundle[0]))
            bundle = '{}'.format(nbundle[0])

        ## find data files
        dfiles = [os.path.basename(f) for f in glob.glob('{}/*.nc'.format(bundle))]
        dfiles.sort()

        ## find xml files
        mfile = [os.path.basename(f) for f in glob.glob('{}/xfdumanifest.xml'.format(bundle))]
        if len(mfile)==0: mfile = [os.path.basename(f) for f in glob.glob('{}/*.xml'.format(bundle))]
        if len(mfile)==1:
            mfile = mfile[0]
        else:
            print('No metadatafile found for {}'.format(bundle))
            continue

        t0 = time.time()

        ## sensor settings
        metafile = '{}/{}'.format(bundle, mfile)
        smeta = ac.sentinel3.metadata_parse(metafile)
        sensor = smeta['sensor']
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
        s3_product_type = None
        for k in ['sentinel3:productType','envisat:productType']:
            if k in smeta:
                if 'FR' in smeta[k]:
                    s3_product_type = 'FR'
                elif 'RR' in smeta[k]:
                    s3_product_type = 'RR'
                break

        ## load rsrd
        rsrd = ac.shared.rsr_dict(sensor)
        rsr_bands = rsrd[sensor]['rsr_bands']

        ## update settings
        setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

        if output is None: output = setu['output']

        ## if not merging tiles, create a new file
        if not setu['merge_tiles']:
            new = True
            sub = None
            data_shape = None
            if setu['limit'] is not None:
                sub, data_shape = ac.sentinel3.olci_sub(bundle, setu['limit'], use_tpg = setu['use_tpg'])
                if sub is None:
                    if setu['verbosity'] > 1: print('Limit {} out of scene {}'.format(setu['limit'], bundle))
                    continue
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
            data, meta = {}, {}

        ## read data
        for f in dfiles:
            fname = os.path.splitext(f)[0]
            file = '{}/{}'.format(bundle, f)
            gem = ac.gem.gem(file)
            gatts = gem.gatts
            start_time = gatts['start_time']
            stop_time = gatts['stop_time']

            ## tie point grids
            if ('tie' in fname) or (fname in ['removed_pixels',
                         'instrument_data', 'time_coordinates']):
                if fname not in meta: meta[fname]={}
                for ds in gem.datasets:
                    tmp = gem.data(ds)
                    if not setu['merge_tiles']:
                        meta[fname][ds] = gem.data(ds)
                    else:
                        ## stack the data if merging tiles
                        if ds not in meta[fname]:
                            meta[fname][ds] = gem.data(ds)
                        else:
                            if len(meta[fname][ds].shape) == 1:
                                ## use hstack if one dimension, in case the number of along track pixels differs
                                meta[fname][ds] = np.hstack((meta[fname][ds], gem.data(ds)))
                            else:
                                meta[fname][ds] = np.vstack((meta[fname][ds], gem.data(ds)))

            ## or full size data
            else:
                for ds in gem.datasets:
                    if ds[-9:] == '_radiance':
                        if not setu['merge_tiles']:
                            data[ds] = gem.data(ds, sub=sub)
                        else:
                            if ds not in data: data[ds] = np.zeros(data_shape) + np.nan
                            data[ds][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = gem.data(ds, sub=sub)
                        if setu['verbosity'] > 2: print(ds, data[ds].shape)

                        ## apply gains (if last tile)
                        if (setu['gains']) & ((not setu['merge_tiles']) | (bi == nscenes-1)):
                            cg = 1.0
                            if len(setu['gains_toa']) == len_gains:
                                gi = int(re.findall(r'\d+', ds)[0])-1
                                cg = float(setu['gains_toa'][gi])
                            if setu['verbosity'] > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                            data[ds]*=cg

                    elif setu['output_geolocation']:
                        data_ = gem.data(ds, sub=sub)
                        if data_shape is None: data_shape = data_.shape
                        ## save latitude and longitude already
                        dk = ds
                        if ds in ['latitude', 'longitude']: dk = ds[0:3]
                        if not setu['merge_tiles']:
                            data[dk] = data_
                        else:
                            if dk not in data: data[dk] = np.zeros(data_shape) + np.nan
                            data[dk][crop_out[bi][2]:crop_out[bi][3], crop_out[bi][0]:crop_out[bi][1]] = data_

        ## end read data

        ## if not final tile set new to False and continue
        if (setu['merge_tiles']) & (bi < len(inputfile)-1):
            new = False
            continue

        ## determine product level
        product_level = 'level1'
        acolite_file_type = 'L1R'
        if 'OLCI Level 2 WATER Product' in gatts['title']:
            product_level = 'level2'
            acolite_file_type = 'L2S'
        if (convert_l2 is False) & (product_level != 'level1'):
            print('File type = {}'.format(gatts['title']))
            print('Not converting to ACOLITE type, set convert_l2 = True.')
            continue

        ## global attributes
        dtime = dateutil.parser.parse(start_time)
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## create full scene tpgs
        ## added half a pixel offset - correspond to centre (and SNAP)
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
        x_sub = gatts['ac_subsampling_factor']
        y_sub = gatts['al_subsampling_factor']
        x_off = 0
        y_off = 0
        tpg_shape = meta['tie_geo_coordinates']['latitude'].shape
        tpx = (np.arange(tpg_shape[1])*x_sub) + x_off
        tpy = (np.arange(tpg_shape[0])*y_sub) + y_off

        ## 2d for RGI
        subx_ = np.tile(subx, (len(suby),1))
        suby_ = (np.tile(suby, len(subx)).reshape(subx_.shape[1], subx_.shape[0])).T

        tpg = {}
        for k in meta.keys():
            if 'tie_' in k:
                for l in meta[k].keys():
                    if l in ['atmospheric_temperature_profile',
                             'horizontal_wind', 'reference_pressure_level']: continue
                    if meta[k][l].shape != tpg_shape:
                        print('{}-{} tpg shape {} not supported'.format(k,l,meta[k][l].shape))
                        continue
                    dtype_in = meta[k][l].dtype
                    rgi = scipy.interpolate.RegularGridInterpolator([tpy, tpx], meta[k][l].astype(np.float64), bounds_error=False, fill_value=None)
                    tpg[l] = rgi((suby_,subx_)).astype(dtype_in)

        ## compute relative azimuth TPG
        tpg['RAA'] = abs(tpg['SAA']-tpg['OAA'])
        tpg['RAA'][tpg['RAA']>180]=np.abs(360-tpg['RAA'][tpg['RAA']>180])
        ## cosine of sun zenith angle
        mu = np.cos(tpg['SZA']*(np.pi/180))

        ## average geometry
        sza = np.nanmean(tpg['SZA'])
        vza = np.nanmean(tpg['OZA'])
        raa = np.nanmean(tpg['RAA'])
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
                di = meta['instrument_data']['detector_index'][sub_merged[1]:sub_merged[1]+sub_merged[3], sub_merged[0]:sub_merged[0]+sub_merged[2]]
            else:
                di = meta['instrument_data']['detector_index']
        elif (sub is None):
            di = meta['instrument_data']['detector_index']
        else:
            di = meta['instrument_data']['detector_index'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
        print(di.shape)
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
                r_ = (data['{}_radiance'.format(band)]) / meta['instrument_data']['solar_flux'][b_i][di]

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
                    r21_diff = (data['{}_radiance'.format(band2)]) / meta['instrument_data']['solar_flux'][b2_i][di]-\
                               (data['{}_radiance'.format(band1)]) / meta['instrument_data']['solar_flux'][b1_i][di]

                    ## wavelength difference ratio
                    wdiff_ratio = (bands_data[band]['wavelength'] - meta['instrument_data']['lambda0'][b_i][di])/\
                                  (meta['instrument_data']['lambda0'][b2_i][di] - meta['instrument_data']['lambda0'][b1_i][di])

                    ## additional smile
                    smile[band] += (r21_diff)*(wdiff_ratio)*(meta['instrument_data']['solar_flux'][b_i][di])
                    del r21_diff, wdiff_ratio

            ## add smile effect to radiance data
            for band in smile: data['{}_radiance'.format(band)]+=smile[band]
            del smile
            ## add back in gas transmittance
            if setu['smile_correction_tgas']:
                if setu['verbosity'] > 2: print('{} - Gas correction restored after smile correction'.format(datetime.datetime.now().isoformat()[0:19]), end='\n')
                for band in bands_data: data['{}_radiance'.format(band)]*=ttg['tt_gas'][band]
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
            for ds in ['lon', 'lat']: gemo.write( ds, data[ds])

        ## output geometry
        if setu['output_geometry']:
            if setu['verbosity'] > 1: print('Writing geometry')
            for k in tpg:
                if k in ['SZA', 'OZA', 'RAA', 'SAA', 'OAA']:
                    if k == 'OZA':
                        ko = 'vza'
                    elif k == 'OAA':
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
                if setu['verbosity'] > 2: print('{} - Reading TOA data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')

                # per pixel wavelength
                l = meta['instrument_data']['lambda0'][iw][di]
                # per pixel f0
                f0 = meta['instrument_data']['solar_flux'][iw][di]
                # if smile corrected use nominal E0
                if setu['smile_correction']: f0 = bands_data[band]['E0']

                # per pixel fwhm
                fwhm = meta['instrument_data']['FWHM'][iw][di]

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

        if (product_level == 'level2'):
            if setu['verbosity'] > 1: print('Writing water reflectance')
            for iw, band in enumerate(rsr_bands):
                wave = waves_names[band]
                ds = 'rhow_{}'.format(wave)
                dname = dnames[iw].replace('radiance', 'reflectance')
                print(dname)
                if dname not in data: continue
                if setu['verbosity'] > 2: print('{} - Reading data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')
                ds_att  = {'wavelength':float(wave)}

                ## write data
                gemo.write(ds, data[dname], ds_att = ds_att)
                if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format(ds, data[dname].shape))

                ## write error dataset
                if write_l2_err:
                    ds = '{}_err'.format(ds)
                    dname = '{}_err'.format(dname)
                    gemo.write(ds, data[dname], ds_att = ds_att)
                    if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format(ds, data[dname].shape))

            ## write other datasets
            for dname in data:
                if '_err' in dname: continue
                if '_reflectance' in dname: continue
                ds = '{}'.format(dname)
                ds_att = None
                if setu['verbosity'] > 2: print('{} - Reading dataset {}'.format(datetime.datetime.now().isoformat()[0:19], ds), end='\n')
                if data[dname].dtype not in (np.float32, np.float64): continue

                ## write data
                gemo.write(ds, data[dname], ds_att = ds_att)
                if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format(ds, data[dname].shape))

                ## write error dataset
                if write_l2_err:
                    ds = '{}_err'.format(ds)
                    dname = '{}_err'.format(dname)
                    if data[dname].dtype not in (np.float32, np.float64): continue
                    gemo.write(ds, data[dname], ds_att = ds_att)
                    if setu['verbosity'] > 2: print('Converting bands: Wrote {} ({})'.format(ds, data[dname].shape))

        ## clear data
        data = None
        gemo.close()
        gemo = None

        if setu['verbosity'] > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
