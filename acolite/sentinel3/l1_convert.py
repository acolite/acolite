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

def l1_convert(inputfile, output = None, settings = {},
                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),
                verbosity = 5):

    import os, glob, datetime, time, re
    import dateutil.parser

    import numpy as np
    import scipy.interpolate
    import acolite as ac

    if 'verbosity' in settings: verbosity = settings['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scene{}'.format(nscenes, 's' if nscenes==1 else ''))

    ## start conversion
    ofile = None
    ofiles = []
    for bundle in inputfile:
        nbundle = glob.glob('{}/*.SEN3'.format(bundle))
        if len(nbundle) == 1:
            print('Found nested SEN3 bundle {}'.format(nbundle[0]))
            bundle = '{}'.format(nbundle[0])

        ## find data files
        dfiles = [os.path.basename(f) for f in glob.glob('{}/*.nc'.format(bundle))]
        dfiles.sort()

        ## find xml files
        mfile = [os.path.basename(f) for f in glob.glob('{}/*.xml'.format(bundle))]
        if len(mfile)==1: mfile = mfile[0]

        t0 = time.time()
        new = True

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

        ## load rsrd
        rsrd = ac.shared.rsr_dict(sensor)
        rsr_bands = rsrd[sensor]['rsr_bands']

        ## merge sensor specific settings
        setu = ac.acolite.settings.parse(sensor, settings=settings)
        verbosity = setu['verbosity']

        ## extract sensor specific settings
        smile_correction = setu['smile_correction']
        use_tpg = setu['use_tpg']
        use_gains = setu['gains']
        gains = setu['gains_toa']

        ## get other settings
        limit = setu['limit']
        output_geolocation = setu['output_geolocation']
        output_geometry = setu['output_geometry']
        vname = setu['region_name']
        output = setu['output']

        ## check if ROI polygon is given
        poly = setu['polygon']
        clip, clip_mask = False, None
        if poly is not None:
            if os.path.exists(poly):
                try:
                    limit = ac.shared.polygon_limit(poly)
                    if verbosity > 1: print('Using limit from polygon envelope: {}'.format(limit))
                    clip = True
                except:
                    if verbosity > 1: print('Failed to import polygon {}'.format(poly))

        ## add limit buffer
        if (limit is not None) & (setu['limit_buffer'] is not None):
            print('Applying limit buffer {}'.format(setu['limit_buffer']))
            print('Old limit: {}'.format(limit))
            setu['limit_old'] = limit
            limit = limit[0] - setu['limit_buffer'], limit[1] - setu['limit_buffer'], \
                    limit[2] + setu['limit_buffer'], limit[3] + setu['limit_buffer']
            print('New limit: {}'.format(limit))

        sub = None
        data_shape = None
        if limit is not None:
            sub, data_shape = ac.sentinel3.olci_sub(bundle, limit, use_tpg=use_tpg)
            if sub is None:
                if verbosity > 1: print('Limit {} out of scene {}'.format(limit, bundle))
                continue

        ## user defined sub if limit or polygon are not set
        if sub is None:
            if 'sub' in setu: sub = setu['sub']

        ## read data
        dshape = None
        lfiles = {}
        data, meta = {}, {}
        for f in dfiles:
            fname = os.path.splitext(f)[0]
            file = '{}/{}'.format(bundle, f)
            datasets = ac.shared.nc_datasets(file)
            gatts = ac.shared.nc_gatts(file)
            start_time = gatts['start_time']
            stop_time = gatts['stop_time']

            ## tie point grids
            if ('tie' in fname) or (fname in ['removed_pixels',
                         'instrument_data', 'time_coordinates']):
                meta[fname]={}
                for ds in datasets:
                    meta[fname][ds]= ac.shared.nc_data(file, ds)
                    if ds == 'detector_index':
                        data[ds] = ac.shared.nc_data(file, ds, sub=sub)
            ## or full size data
            else:
                for ds in datasets:
                    if ds[-9:] == '_radiance':
                        data[ds] = ac.shared.nc_data(file, ds, sub=sub)
                        if verbosity > 2: print(ds, data[ds].shape)
                        if use_gains:
                            cg = 1.0
                            if len(gains) == len_gains:
                                gi = int(re.findall(r'\d+', ds)[0])-1
                                cg = float(gains[gi])
                            if verbosity > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                            data[ds]*=cg
                    elif output_geolocation:
                        data_ = ac.shared.nc_data(file, ds, sub=sub)
                        if dshape is None: dshape = data_.shape
                        if data_shape is None: data_shape = data_.shape
                        ## save latitude and longitude already
                        if ds in ['latitude', 'longitude']:
                            data[ds[0:3]] = data_
                        else:
                            data[ds] = data_
        ## end read data

        ## global attributes
        dtime = dateutil.parser.parse(start_time)
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## create full scene tpgs
        full_scene = False
        if sub is None:
            full_scene = True
            ## added half a pixel offset - correspond to centre (and SNAP)
            subx = np.arange(0,data_shape[1])+0.5
            suby = np.arange(0,data_shape[0])+0.5
        else:
            ## added half a pixel offset - correspond to centre (and SNAP)
            subx = np.arange(sub[0],sub[0]+sub[2])+0.5
            suby = np.arange(sub[1],sub[1]+sub[3])+0.5

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
                anc = ac.ac.ancillary.get(dtime, clon, clat)
                ## overwrite the defaults
                if ('ozone' in anc): setu['uoz_default'] = anc['ozone']['interp']/1000. ## convert from MET data
                if ('p_water' in anc): setu['uwv_default'] = anc['p_water']['interp']/10. ## convert from MET data
                if ('z_wind' in anc) & ('m_wind' in anc) & (setu['wind'] is None):
                    setu['wind'] = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
                if ('press' in anc) & (setu['pressure'] == setu['pressure_default']):
                    setu['pressure'] = anc['press']['interp']

        if uoz is None: uoz = setu['uoz_default']
        if uwv is None: uwv = setu['uwv_default']
        if pressure is None: pressure = setu['pressure']
        print('current uoz: {:.2f} uwv: {:.2f} pressure: {:.2f}'.format(uoz, uwv, pressure))

        ## for smile correction
        ttg = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv, sensor=sensor)

        ## get per pixel detector index
        if sub is None:
            di = meta['instrument_data']['detector_index']
        else:
            di = meta['instrument_data']['detector_index'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        ## smile correction - from l2gen smile.c
        if smile_correction:
            if verbosity > 1: print('Running smile correction')
            ## gas_correction
            if setu['smile_correction_tgas']:
                if verbosity > 2: print('{} - Gas correction before smile correction'.format(datetime.datetime.now().isoformat()[0:19]), end='\n')
                for band in bands_data: data['{}_radiance'.format(band)]/=ttg['tt_gas'][band]

            smile = {}
            for band in bands_data:
                ## dataset name
                dname = '{}_radiance'.format(band)
                ## band index
                b_i = bands_data[band]['band']-1
                if verbosity > 2: print('{} - Smile correction for band {} {} nm'.format(datetime.datetime.now().isoformat()[0:19], band, bands_data[band]['wavelength'] ), end='\n')

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
                    if verbosity > 2: print('{} - Smile correction - bounding bands {}/{}'.format(datetime.datetime.now().isoformat()[0:19], band1, band2), end='\n')

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
                if verbosity > 2: print('{} - Gas correction restored after smile correction'.format(datetime.datetime.now().isoformat()[0:19]), end='\n')
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
                     'se_distance': se_distance, 'acolite_file_type': 'L1R'}
        gatts['pressure'] = pressure
        gatts['uoz'] = uoz
        gatts['uwv'] = uwv

        if limit is not None: gatts['limit'] = limit
        if sub is not None: gatts['sub'] = sub

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if vname != '': oname+='_{}'.format(vname)

        ofile = '{}/{}_L1R.nc'.format(output, oname)
        if not os.path.exists(os.path.dirname(ofile)): os.makedirs(os.path.dirname(ofile))
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## add band info to gatts
        for bi, b in enumerate(rsr_bands):
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## output lat/lon
        if output_geolocation:
            if verbosity > 1: print('Writing geolocation')
            for ds in ['lon', 'lat']:
                ac.output.nc_write(ofile, ds, data[ds], new=new, attributes=gatts,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                new = False

        ## output geometry
        if output_geometry:
            if verbosity > 1: print('Writing geometry')
            for k in tpg:
                if k in ['SZA', 'OZA', 'RAA', 'SAA', 'OAA']:
                    if k == 'OZA':
                        ko = 'vza'
                    elif k == 'OAA':
                        ko = 'vaa'
                    else:
                        ko = k.lower()
                    ac.output.nc_write(ofile, ko, tpg[k], new=new, attributes=gatts,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                    new = False
                elif k in ['sea_level_pressure']:
                    ac.output.nc_write(ofile, 'pressure', tpg[k], new=new, attributes=gatts,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'],
                                    netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                    new = False
                else:
                    continue

        ## is the F0 already corrected for Sun - Earth distance?
        ## see OLCI_L2_ATBD_Instrumental_Correction.pdf
        #se = ac.shared.distance_se(doy)
        se = 1.
        se2 = se**2

        ## read TOA
        if verbosity > 1: print('Writing TOA reflectance')
        for iw, band in enumerate(rsr_bands):
            wave = waves_names[band]
            ds = 'rhot_{}'.format(wave)
            if verbosity > 2: print('{} - Reading TOA data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')

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
                ac.output.nc_write(ofile, 'Lt_{}'.format(wave), data[dname],
                              dataset_attributes = ds_att,
                              netcdf_compression=setu['netcdf_compression'],
                              netcdf_compression_level=setu['netcdf_compression_level'],
                              netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                if verbosity > 2: print('Converting bands: Wrote {} ({})'.format('Lt_{}'.format(wave), data[dname].shape))

            ## convert to reflectance
            d = (np.pi * data[dname] * se2) / (f0*mu)
            mask = d.mask
            d = d.data
            d[mask] = np.nan

            ac.output.nc_write(ofile, ds, d, dataset_attributes=ds_att, new=new, attributes=gatts,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'],
                                netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            if verbosity > 2: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
            new = False
            d = None

        ## clear data
        data = None

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
