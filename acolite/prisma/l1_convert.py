## def l1_convert
## converts PRISMA HDF file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-07-14
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-02-23 (QV) added option to output L2C reflectances

def l1_convert(inputfile, output=None, settings = {}, verbosity=0):
    import numpy as np
    import h5py, dateutil.parser, os
    import acolite as ac

    ## parse settings
    setu = ac.acolite.settings.parse('PRISMA', settings=settings)
    verbosity = setu['verbosity']
    if output is None: output = setu['output']
    output_lt = setu['output_lt']
    vname = setu['region_name']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get()

    ofiles = []
    for file in inputfile:
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

        idx = np.argsort(waves)
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)

        idx = np.argsort(waves)
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

        ## get geometry from l2 file if present
        l2file = os.path.dirname(file) + os.path.sep + os.path.basename(file).replace('PRS_L1_STD_OFFL_', 'PRS_L2C_STD_')
        if os.path.exists(l2file):
            with h5py.File(l2file, mode='r') as f:
                vza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Observing_Angle'][:]
                raa = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Rel_Azimuth_Angle'][:]
                sza = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Geometric Fields']['Solar_Zenith_Angle'][:]
                gatts['vza'] = np.nanmean(np.abs(vza))
                gatts['raa'] = np.nanmean(np.abs(raa))
                gatts['sza'] = np.nanmean(np.abs(sza))
        else:
            print('PRISMA processing only supported when L2 geometry is present.')
            print('Please put {} in the same directory as {}'.format(os.path.basename(l2file), os.path.basename(file)))
            continue

        src = 'HCO'
        read_cube = True
        with h5py.File(file, mode='r') as f:
            lat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Latitude_SWIR'][:]
            lon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Longitude_SWIR'][:]
            ## read bands in spectral order
            if read_cube:
                vnir_data = h5_gatts['Offset_Vnir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:]/h5_gatts['ScaleFactor_Vnir']
                swir_data = h5_gatts['Offset_Swir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:]/h5_gatts['ScaleFactor_Swir']
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

        if output is None:
            odir = os.path.dirname(file)
        else:
            odir = output

        obase  = '{}_{}_L1R'.format('PRISMA',  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)

        gatts['obase'] = obase
        gatts['sensor'] = 'PRISMA'
        gatts['isodate'] = time.isoformat()

        gatts['band_waves'] = [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = [bands[w]['width'] for w in bands]

        ac.output.nc_write(ofile, 'lat', np.flip(np.rot90(lat)), new=True, attributes=gatts,
                            netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
        ac.output.nc_write(ofile, 'lon', np.flip(np.rot90(lon)),
                            netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
        if os.path.exists(l2file):
            ac.output.nc_write(ofile, 'sza', np.flip(np.rot90(sza)),
                                netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
            ac.output.nc_write(ofile, 'vza', np.flip(np.rot90(vza)),
                                netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
            ac.output.nc_write(ofile, 'raa', np.flip(np.rot90(raa)),
                                netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])

        ## store l2c data
        store_l2c = setu['prisma_store_l2c']
        store_l2c_separate_file = setu['prisma_store_l2c_separate_file']
        if store_l2c & read_cube:
            if store_l2c_separate_file:
                obase_l2c  = '{}_{}_converted_L2C'.format('PRISMA',  time.strftime('%Y_%m_%d_%H_%M_%S'))
                ofile_l2c = '{}/{}.nc'.format(odir, obase_l2c)
                ac.output.nc_write(ofile_l2c, 'lat', np.flip(np.rot90(lat)), new=True, attributes=gatts,
                                    netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
                ac.output.nc_write(ofile_l2c, 'lon', np.flip(np.rot90(lon)),
                                    netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
            else:
                ofile_l2c = '{}'.format(ofile)

            ## get l2c details for reflectance conversion
            h5_l2c_gatts = ac.prisma.attributes(l2file)
            scale_max = h5_l2c_gatts['L2ScaleVnirMax']
            scale_min = h5_l2c_gatts['L2ScaleVnirMin']

            ##  read in data cube
            with h5py.File(l2file, mode='r') as f:
                vnir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['VNIR_Cube'][:]
                swir_l2c_data = f['HDFEOS']['SWATHS']['PRS_L2C_HCO']['Data Fields']['SWIR_Cube'][:]

        ## write TOA data
        for bi, b in enumerate(bands):
            wi = bands[b]['index']
            i = bands[b]['i']
            print('Reading rhot_{}'.format(bands[b]['wave_name']))

            if bands[b]['instrument'] == 'vnir':
                if read_cube:
                    cdata_radiance = vnir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                    if store_l2c:
                        cdata_l2c = scale_min + (vnir_l2c_data[:, wi, :] * (scale_max - scale_min)) / 65535
                else:
                    cdata_radiance = h5_gatts['Offset_Vnir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Vnir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            if bands[b]['instrument'] == 'swir':
                if read_cube:
                    cdata_radiance = swir_data[:,wi,:]
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
                    if store_l2c:
                        cdata_l2c = scale_min + (swir_l2c_data[:, wi, :] * (scale_max - scale_min)) / 65535
                else:
                    cdata_radiance = h5_gatts['Offset_Swir'] + \
                            f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Swir']
                    cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)

            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            if output_lt:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata_radiance)),
                              dataset_attributes = ds_att,
                              netcdf_compression=setu['netcdf_compression'],
                              netcdf_compression_level=setu['netcdf_compression_level'],
                              netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                cdata_radiance = None

            ## write toa reflectance
            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata)),
                              dataset_attributes = ds_att,
                              netcdf_compression=setu['netcdf_compression'],
                              netcdf_compression_level=setu['netcdf_compression_level'],
                              netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            cdata = None
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))

            ## store L2C data
            if store_l2c & read_cube:
                ac.output.nc_write(ofile_l2c, 'rhos_l2c_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata_l2c)),
                                  dataset_attributes = ds_att,
                                  netcdf_compression=setu['netcdf_compression'],
                                  netcdf_compression_level=setu['netcdf_compression_level'],
                                  netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                ofile_l2c_new = False
                cdata_l2c = None
                print('Wrote rhos_l2c_{}'.format(bands[b]['wave_name']))

        ofiles.append(ofile)
    return(ofiles, setu)
