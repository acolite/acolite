## def l1_convert
## converts Sentinel-3/OLCI or ENVISAT/MERIS data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-24
## modifications: 2021-12-22 (QV) added MERIS processing
##                2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression

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
        t0 = time.time()
        new = True
        ## identify sensor
        platform = os.path.basename(bundle)[0:3]
        if platform in ['S3A', 'S3B']:
            sensor = '{}_OLCI'.format(platform)
            len_gains = 21
            bands_data = ac.sentinel3.olci_band_info()
            band_id = 'Oa'
        if platform in ['EN1']:
            sensor = '{}_MERIS'.format(platform)
            len_gains = 15
            bands_data = ac.sentinel3.meris_band_info()
            band_id = 'M'
        rsr_file = ac.config['data_dir']+'/RSR/'+sensor+'.txt'
        rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

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

        ## find data files
        dfiles = [os.path.basename(f) for f in glob.glob('{}/*.nc'.format(bundle))]
        dfiles.sort()

        ## find xml files
        mfile = [os.path.basename(f) for f in glob.glob('{}/*.xml'.format(bundle))]
        if len(mfile)==1: mfile = mfile[0]

        sub = None
        data_shape = None
        if limit is not None:
            sub, data_shape = ac.sentinel3.olci_sub(bundle, limit, use_tpg=use_tpg)
            if sub is None:
                if verbosity > 1: print('Limit {} out of scene {}'.format(limit, bundle))
                continue

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
                    if ds[3:] == '_radiance':
                        data[ds] = ac.shared.nc_data(file, ds, sub=sub)
                        if verbosity > 2: print(ds, data[ds].shape)
                        if use_gains:
                            cg = 1.0
                            if len(gains) == len_gains:
                                gi = int(re.findall(r'\d+', ds)[0])-1
                                cg = float(gains[gi])
                            if verbosity > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                            data[ds]*=cg
                        continue

                    if output_geolocation:
                        data_ = ac.shared.nc_data(file, ds, sub=sub)
                        if dshape is None: dshape = data_.shape
                        if data_shape is None: data_shape = data_.shape
                        ## save latitude and longitude already
                        if ds in ['latitude', 'longitude']:
                            data[ds[0:3]] = data_
                        else:
                            data[ds] = data_
        ## end read data

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
        tpg_shape = meta['tie_geo_coordinates']['latitude'].shape
        tpx = np.linspace(0,data_shape[1]-1,tpg_shape[1])
        tpy = np.linspace(0,data_shape[0]-1,tpg_shape[0])
        tpg = {}
        for k in meta.keys():
            if 'tie_' in k:
                for l in meta[k].keys():
                    if l in ['atmospheric_temperature_profile',
                             'horizontal_wind', 'reference_pressure_level']: continue
                    if meta[k][l].shape != tpg_shape:
                        print('{}-{} tpg shape {} not supported'.format(k,l,meta[k][l].shape))
                        continue
                    z = scipy.interpolate.interp2d(tpx, tpy, meta[k][l])
                    tpg[l] = z(subx,suby)

        ## compute relative azimuth TPG
        tpg['RAA'] = abs(tpg['SAA']-tpg['OAA'])
        tpg['RAA'][tpg['RAA']>180]=np.abs(360-tpg['RAA'][tpg['RAA']>180])
        ## cosine of sun zenith angle
        mu = np.cos(tpg['SZA']*(np.pi/180))


        ## average geometry
        sza = np.nanmean(tpg['SZA'])
        vza = np.nanmean(tpg['OZA'])
        raa = np.nanmean(tpg['RAA'])

        ## compute gas transmittance
        use_supplied_ancillary = True
        uoz_default=0.3
        uwv_default=1.5
        if use_supplied_ancillary:
            ## convert ozone from kg.m-2 to cm.atm
            uoz = np.nanmean(tpg['total_ozone'])/0.02141419
            ## convert water from kg.m-2 to g.cm-2
            uwv = np.nanmean(tpg['total_columnar_water_vapour'])/10
        else:
            ## can get other ancillary here
            uoz = uoz_default
            uwv = uwv_default

        ## for smile correction
        ttg = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv, sensor=sensor)

        ## get per pixel detector index
        if sub is None:
            di = meta['instrument_data']['detector_index']
        else:
            di = meta['instrument_data']['detector_index'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        ## smile correction - from l2gen smile.c
        if verbosity > 1: print('Running smile correction')
        if smile_correction:
            smile = {}
            for band in bands_data:
                ## dataset name
                dname = '{}_radiance'.format(band)

                ## band index
                b_i = bands_data[band]['band']-1

                ## gas_correction:
                data['{}_radiance'.format(band)]/=ttg['tt_gas'][band]

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
        ## end smile correction

        ## global attributes
        dtime = dateutil.parser.parse(start_time)
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## read rsr
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## take wavelengths and band names from external table
        ## the smile correction should bring things in line with these
        #waves_names = ['{:.0f}'.format(bands_data[b]['wavelength']) for b in bands_data]
        dnames = ['{}_radiance'.format(b) for b in bands_data]
        bnames = [b for b in bands_data]

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        gatts = {'sensor':sensor, 'sza':sza, 'vza':vza, 'raa':raa,
                     'isodate':isodate, 'global_dims':data_shape,
                     'se_distance': se_distance, 'acolite_file_type': 'L1R'}

        if limit is not None: gatts['limit'] = limit

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
            # per pixel fwhm
            fwhm = meta['instrument_data']['FWHM'][iw][di]

            dname = dnames[iw]

            ## convert to reflectance
            d = (np.pi * data[dname] * se2) / (f0*mu)
            mask = d.mask
            d = d.data
            d[mask] = np.nan

            ds_att  = {'wavelength':float(wave)}
            for key in ttg: ds_att[key]=ttg[key][bnames[iw]]

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
