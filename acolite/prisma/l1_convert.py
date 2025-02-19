## def l1_convert
## converts PRISMA HDF file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-07-14
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-02-23 (QV) added option to output L2C reflectances
##                2023-05-09 (QV) added option to crop
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2024-05-11 (QV) update for L1G2 data
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output=None, settings = None):
    import numpy as np
    import h5py, dateutil.parser, os
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'PRISMA'

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    verbosity = setu['verbosity']
    if output is None: output = setu['output']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

    ofiles = []
    for file in inputfile:
        if output is None: output = os.path.dirname(file)

        #f = h5py.File(file, mode='r')
        h5_gatts = ac.prisma.attributes(file)

        waves_vnir = h5_gatts['List_Cw_Vnir']
        bands_vnir = ['{:.0f}'.format(w) for w in waves_vnir]
        fwhm_vnir = h5_gatts['List_Fwhm_Vnir']
        n_vnir = len(waves_vnir)

        waves_swir = h5_gatts['List_Cw_Swir']
        bands_swir = ['{:.0f}'.format(w) for w in waves_swir]
        fwhm_swir = h5_gatts['List_Fwhm_Swir']
        n_swir = len(waves_swir)

        waves = [w for w in waves_vnir] + [w for w in waves_swir]
        fwhm = [f for f in fwhm_vnir] + [f for f in fwhm_swir]
        waves_names = ['{:.0f}'.format(w) for w in waves]
        instrument = ['vnir']*n_vnir + ['swir']*n_swir
        band_index = [i for i in range(n_vnir)] + [i for i in range(n_swir)]

        band_names_vnir = ['vnir_{}'.format(b) for b in range(0, n_vnir)]
        band_names_swir = ['swir_{}'.format(b) for b in range(0, n_swir)]

        rsr_vnir = {'vnir_{}'.format(b): ac.shared.gauss_response(waves_vnir[b], fwhm_vnir[b], step=0.1) for b in range(0, n_vnir)}
        rsr_swir = {'swir_{}'.format(b): ac.shared.gauss_response(waves_swir[b], fwhm_swir[b], step=0.1) for b in range(0, n_swir)}

        band_names = band_names_vnir + band_names_swir
        band_rsr = {}
        for b in rsr_vnir: band_rsr[b] = {'wave': rsr_vnir[b][0]/1000, 'response': rsr_vnir[b][1]}
        for b in rsr_swir: band_rsr[b] = {'wave': rsr_swir[b][0]/1000, 'response': rsr_swir[b][1]}

        ## use same rsr as acolite_l2r
        #rsr = ac.shared.rsr_hyper(gatts['band_waves'], gatts['band_widths'], step=0.1)
        # rsrd = ac.shared.rsr_dict(rsrd={sensor:{'rsr':band_rsr}})
        # waves = [rsrd[sensor]['wave_nm'][b] for b in band_names]
        # waves_names = [rsrd[sensor]['wave_name'][b] for b in band_names]

        idx = np.argsort(waves)
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)

        bands = {}
        for i in idx:
            cwave = waves[i]
            if cwave == 0: continue
            swave = '{:.0f}'.format(cwave)
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':waves_names[i],
                           'width': fwhm[i],
                           'i':i, 'index':band_index[i],
                           'rsr': band_rsr[band_names[i]],
                           'f0': f0d[band_names[i]],
                           'instrument':instrument[i],}

        gatts = {}
        isotime = h5_gatts['Product_StartTime']
        time = dateutil.parser.parse(isotime)

        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)

        ## lon and lat keys
        bn = os.path.basename(file)
        l1g = ('PRS_L1G_STD_OFFL_' in bn) | ('PRS_L1G2_STD_OFFL_' in bn)
        if l1g:
            lat_key = 'Latitude'
            lon_key = 'Longitude'
        else:
            lat_key = 'Latitude_SWIR'
            lon_key = 'Longitude_SWIR'

        ## mask for L1G format
        mask_value = 65535
        dem = None

        ## reading settings
        src = 'HCO' ## coregistered radiance cube
        read_cube = True

        ## get geometry from l2 file if present
        if ac.settings['run']['l2cfile'] is not None:
            l2file = ac.settings['run']['l2cfile']
            print('Using user specified L2C file {}'.format(l2file))
        else:
            l2file = os.path.dirname(file) + os.path.sep + os.path.basename(file).replace('PRS_L1_STD_OFFL_', 'PRS_L2C_STD_')

        if not os.path.exists(l2file):
            print('PRISMA processing only supported when L2 geometry is present.')
            print('Please put {} in the same directory as {}'.format(os.path.basename(l2file), os.path.basename(file)))
            continue

        ## read geolocation
        with h5py.File(file, mode='r') as f:
            lat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields'][lat_key][:]
            lon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields'][lon_key][:]
            lat[lat>=mask_value] = np.nan
            lon[lon>=mask_value] = np.nan
        sub = None
        if setu['limit'] is not None:
            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
            if sub is None:
                print('Limit outside of scene {}'.format(file))
                continue
            ## crop to sub
            lat = lat[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
            lon = lon[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
        ## end read geolocation

        ## read geometry
        vza, vaa, sza, saa, raa = None, None, None, None, None
        with h5py.File(l2file, mode='r') as f:
            ## L1G format
            if l1g:
                if sub is None:
                    vza = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Sensor_Zenith_Angle'][:]
                    vaa = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Sensor_Azimuth_Angle'][:]
                    sza = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Solar_Zenith_Angle'][:]
                    saa = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Solar_Azimuth_Angle'][:]
                else:
                    vza = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Sensor_Zenith_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    vaa = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Sensor_Azimuth_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    sza = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Solar_Zenith_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    saa = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Geometric Fields']['Solar_Azimuth_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

                ## apply mask
                vza[vza>=mask_value] = np.nan
                sza[sza>=mask_value] = np.nan
                saa[saa>=mask_value] = np.nan
                vaa[vaa>=mask_value] = np.nan

                ## compute relative azimuth
                raa = np.abs(saa - vaa)
                raa[raa>180] = 360 - raa[raa>180]

                ## get DEM data
                if sub is None:
                    dem = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Terrain Fields']['DEM'][:]
                else:
                    dem = f['HDFEOS']['SWATHS']['PRS_L1_HCO']['Terrain Fields']['DEM'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                dem[dem>=mask_value] = np.nan
            else:
                if sub is None:
                    vza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Observing_Angle'][:]
                    raa = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Rel_Azimuth_Angle'][:]
                    sza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Solar_Zenith_Angle'][:]
                else:
                    vza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Observing_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    raa = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Rel_Azimuth_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                    sza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Solar_Zenith_Angle'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

        gatts['vza'] = np.nanmean(np.abs(vza))
        gatts['raa'] = np.nanmean(np.abs(raa))
        gatts['sza'] = np.nanmean(np.abs(sza))

        with h5py.File(file, mode='r') as f:
            ## read bands in spectral order
            if read_cube:
                if sub is None:
                    vnir_data = h5_gatts['Offset_Vnir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:]/h5_gatts['ScaleFactor_Vnir']
                    swir_data = h5_gatts['Offset_Swir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:]/h5_gatts['ScaleFactor_Swir']
                else:
                    vnir_data = h5_gatts['Offset_Vnir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][sub[1]:sub[1]+sub[3], :, sub[0]:sub[0]+sub[2]]/h5_gatts['ScaleFactor_Vnir']
                    swir_data = h5_gatts['Offset_Swir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][sub[1]:sub[1]+sub[3], :, sub[0]:sub[0]+sub[2]]/h5_gatts['ScaleFactor_Swir']

                vnir_data[vnir_data>=mask_value] = np.nan
                swir_data[swir_data>=mask_value] = np.nan

            ## read LOS vectors
            x_ = f['KDP_AUX']['LOS_Vnir'][:,0]
            y_ = f['KDP_AUX']['LOS_Vnir'][:,1]
            z_ = f['KDP_AUX']['LOS_Vnir'][:,2]

        ## get vza/vaa
        #dtor = np.pi/180
        #vza = np.arctan2(y_,x_)/dtor
        #vaa = np.arctan2(z_,np.sqrt(x_**2+y_**2))/dtor
        #vza_ave = np.nanmean(np.abs(vza))
        #vaa_ave = np.nanmean(np.abs(vaa))
        sza_ave = h5_gatts['Sun_zenith_angle']
        saa_ave = h5_gatts['Sun_azimuth_angle']

        if setu['prisma_rhot_per_pixel_sza']:
            cossza = np.cos(np.radians(sza))
        else:
            cossza = np.cos(np.radians(sza_ave))

        vza_ave = 0
        vaa_ave = 0

        if 'sza' not in gatts: gatts['sza'] = sza_ave
        if 'vza' not in gatts: gatts['vza'] = vza_ave
        if 'saa' not in gatts: gatts['saa'] = saa_ave
        if 'vaa' not in gatts: gatts['vaa'] = vaa_ave

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        gatts['sensor'] = sensor
        gatts['isodate'] = time.isoformat()
        gatts['acolite_file_type'] = 'L1R'
        gatts['band_waves'] = [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = [bands[w]['width'] for w in bands]

        ## output file name
        oname  = '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## set up output file
        gemo = ac.gem.gem(ofile, new=True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        if (setu['output_geolocation']):
            if verbosity > 1: print('Writing geolocation lon/lat')
            gemo.write('lon', np.flip(np.rot90(lon)))
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
            if not (setu['prisma_store_l2c'] & setu['prisma_store_l2c_separate_file']): lon = None
            gemo.write('lat', np.flip(np.rot90(lat)))
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
            if not (setu['prisma_store_l2c'] & setu['prisma_store_l2c_separate_file']): lat = None

        ## write geometry
        if os.path.exists(l2file):
            if (setu['output_geometry']):
                if verbosity > 1: print('Writing geometry')
                gemo.write('vza', np.flip(np.rot90(vza)))
                if verbosity > 1: print('Wrote vza ({})'.format(vza.shape))
                vza = None
                new = False
                if vaa is not None:
                    gemo.write('vaa', np.flip(np.rot90(vaa)))
                    if verbosity > 1: print('Wrote vaa ({})'.format(vaa.shape))
                    vaa = None
                gemo.write('sza', np.flip(np.rot90(sza)))
                if verbosity > 1: print('Wrote sza ({})'.format(sza.shape))

                if saa is not None:
                    gemo.write('saa', np.flip(np.rot90(saa)))
                    if verbosity > 1: print('Wrote saa ({})'.format(saa.shape))
                    saa = None
                gemo.write('raa', np.flip(np.rot90(raa)))
                if verbosity > 1: print('Wrote raa ({})'.format(raa.shape))
                raa = None

            if dem is not None:
                gemo.write('dem', np.flip(np.rot90(dem)))

        ## store l2c data
        if setu['prisma_store_l2c'] & read_cube:
            if setu['prisma_store_l2c_separate_file']:
                oname_l2c  = '{}_{}'.format('PRISMA',  time.strftime('%Y_%m_%d_%H_%M_%S'))
                if setu['region_name'] != '': oname_l2c+='_{}'.format(setu['region_name'])
                ofile_l2c = '{}/{}_L2C.nc'.format(output, oname_l2c)
                gemo_l2c = ac.gem.gem(ofile_l2c, new=True)
                gemo_l2c.gatts = {k: gatts[k] for k in gatts}
                gemo_l2c.gatts['oname'] = oname_l2c
                gemo_l2c.gatts['ofile'] = ofile_l2c
                gemo_l2c.gatts['acolite_file_type'] = 'L2C'
                gemo_l2c.write('lon', np.flip(np.rot90(lon)))
                gemo_l2c.write('lat', np.flip(np.rot90(lat)))
                lat = None
                lon = None
            else:
                ofile_l2c = '{}'.format(ofile)

            ## get l2c details for reflectance conversion
            h5_l2c_gatts = ac.prisma.attributes(l2file)
            scale_max = h5_l2c_gatts['L2ScaleVnirMax']
            scale_min = h5_l2c_gatts['L2ScaleVnirMin']

            ##  read in data cube
            with h5py.File(l2file, mode='r') as f:
                if sub is None:
                    vnir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['VNIR_Cube'][:]
                    swir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['SWIR_Cube'][:]
                else:
                    vnir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['VNIR_Cube'][sub[1]:sub[1]+sub[3], :, sub[0]:sub[0]+sub[2]]
                    swir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['SWIR_Cube'][sub[1]:sub[1]+sub[3], :, sub[0]:sub[0]+sub[2]]

        ## write TOA data
        for bi, b in enumerate(bands):
            wi = bands[b]['index']
            i = bands[b]['i']
            print('Reading rhot_{}'.format(bands[b]['wave_name']))

            if bands[b]['instrument'] == 'vnir':
                if read_cube:
                    cdata_radiance = vnir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                    if setu['prisma_store_l2c']:
                        cdata_l2c = scale_min + (vnir_l2c_data[:, wi, :] * (scale_max - scale_min)) / 65535
                else:
                    if sub is None:
                        cdata_radiance = h5_gatts['Offset_Vnir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Vnir']
                    else:
                        cdata_radiance = h5_gatts['Offset_Vnir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][sub[1]:sub[1]+sub[3], i, sub[0]:sub[0]+sub[2]]/h5_gatts['ScaleFactor_Vnir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            if bands[b]['instrument'] == 'swir':
                if read_cube:
                    cdata_radiance = swir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                    if setu['prisma_store_l2c']:
                        cdata_l2c = scale_min + (swir_l2c_data[:, wi, :] * (scale_max - scale_min)) / 65535
                else:
                    if sub is None:
                        cdata_radiance = h5_gatts['Offset_Swir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Swir']
                    else:
                        cdata_radiance = h5_gatts['Offset_Swir'] + \
                                f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][sub[1]:sub[1]+sub[3], i, sub[0]:sub[0]+sub[2]]/h5_gatts['ScaleFactor_Swir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']),
                                    np.flip(np.rot90(cdata_radiance)),ds_att = ds_att)
                cdata_radiance = None

            ## write toa reflectance
            gemo.write('rhot_{}'.format(bands[b]['wave_name']),
                                    np.flip(np.rot90(cdata)), ds_att = ds_att)
            cdata = None
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))

            ## store L2C data
            if setu['prisma_store_l2c'] & read_cube:
                if ofile_l2c != ofile:
                    gemo_l2c.write('rhos_l2c_{}'.format(bands[b]['wave_name']),
                                        np.flip(np.rot90(cdata_l2c)), ds_att = ds_att)
                else:
                    gemo.write('rhos_l2c_{}'.format(bands[b]['wave_name']),
                                        np.flip(np.rot90(cdata_l2c)), ds_att = ds_att)

                ofile_l2c_new = False
                cdata_l2c = None
                print('Wrote rhos_l2c_{}'.format(bands[b]['wave_name']))

        ## close files
        gemo.close()
        gemo = None
        if setu['prisma_store_l2c'] & setu['prisma_store_l2c_separate_file']:
            gemo_l2c.close()
            gemo_l2c = None

        ## output PAN
        if setu['prisma_output_pan']:
            psrc = src.replace('H', 'P')
            with h5py.File(file, mode='r') as f:
                if sub is None:
                    pan = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Data Fields']['Cube'][:]
                    plat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Geolocation Fields']['Latitude'][:]
                    plon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Geolocation Fields']['Longitude'][:]
                else:
                    psub = [s*6 for s in sub]
                    pan = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Data Fields']['Cube'][psub[1]:psub[1]+psub[3], psub[0]:psub[0]+psub[2]]
                    plat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Geolocation Fields']['Latitude'][psub[1]:psub[1]+psub[3], psub[0]:psub[0]+psub[2]]
                    plon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(psrc)]['Geolocation Fields']['Longitude'][psub[1]:psub[1]+psub[3], psub[0]:psub[0]+psub[2]]

            ## convert to radiance
            pan = h5_gatts['Offset_Pan'] + pan / h5_gatts['ScaleFactor_Pan']

            ## output netcdf
            ofile_pan = '{}/{}_pan.nc'.format(output, oname)
            gemo_pan = ac.gem.gem(ofile_pan, new=True)
            gemo_pan.write('lon', np.flip(np.rot90(plon)))
            plon = None
            gemo_pan.write('lat', np.flip(np.rot90(plat)))
            plat = None
            gemo_pan.write('pan', np.flip(np.rot90(pan)))
            pan = None
            gemo_pan.close()
            gemo_pan = None
         ## end PAN

        ofiles.append(ofile)
    return(ofiles, setu)
